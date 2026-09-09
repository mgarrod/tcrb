#!/usr/bin/env python3
"""
update_data.py — incremental refresh for the T CrB ELL O-C / light-curve pipeline.

What it does, each run:
  1. Reads the existing photometry CSV and finds the latest date already present.
  2. Downloads only AAVSO Johnson V observations newer than that (via the AAVSO
     v2 API), merges them into the CSV (deduped by observation id).
  3. Re-bins the full 2005-present light curve and re-fits every ELL primary-
     minimum epoch.
  4. Writes two JSON files that the HTML pages fetch at load time:
       - oc_data.json          -> O-C diagram
       - lightcurve_data.json  -> observed-vs-calculated light curve
  5. Prints a summary of what changed since the previous run.

Requires an AAVSO API key (Account Settings -> API Key on aavso.org), passed via
the AAVSO_API_KEY environment variable or --api-key. Without a key, the script
still re-bins/re-fits and regenerates the JSON from whatever is already in the
CSV — it just can't pull new observations.

Usage:
    AAVSO_API_KEY=xxxx python3 update_data.py
    python3 update_data.py --api-key xxxx --csv ../TCrB_V_1946-2026.csv

Everything the script needs — CSV I/O, nightly binning, the weighted quadratic
minimum-fit, epoch scanning — lives in this one file. There's nothing else that
imports this logic, so it isn't split into a separate library module.
"""
import argparse
import csv
import io
import json
import math
import os
import sys
import zipfile
from collections import defaultdict
from datetime import datetime, timedelta, timezone

import requests

# ======================================================================
# Config
# ======================================================================

DEFAULT_CSV = os.path.join(os.path.dirname(__file__), '..', 'TCrB_V_1946-2026.csv')
DEFAULT_OUT_DIR = os.path.dirname(__file__)

P_DEFAULT = 227.5528
T0_DEFAULT = 2455828.9  # HJD 2455828.9 = 2011-09-24, per VSX oid=10602

TARGET = 't crb'
BAND = 2  # Johnson V
HALF_WINDOW = 50.0
MIN_NIGHTS = 10
LIGHTCURVE_START = (2005, 1, 1)
API_BASE = 'https://apps.aavso.org/v2/api/observations/photometry/download/'

JD_UNIX_EPOCH = 2440587.5  # JD at 1970-01-01 00:00 UT

CSV_FIELDS = ['#', 'target', 'auid', 'jd', 'mag', 'uncertainty', 'fainterthan',
              'band', 'type', 'observer', 'airmass', 'transformed']


# ======================================================================
# Time utils
# ======================================================================

def jd_to_datetime(jd):
    return datetime(1970, 1, 1, tzinfo=timezone.utc) + timedelta(days=jd - JD_UNIX_EPOCH)


def jd_to_date_str(jd):
    return jd_to_datetime(jd).strftime('%Y-%m-%d')


def jd_to_decyear(jd):
    d = jd_to_datetime(jd)
    ys = datetime(d.year, 1, 1, tzinfo=timezone.utc)
    ye = datetime(d.year + 1, 1, 1, tzinfo=timezone.utc)
    return d.year + (d - ys).total_seconds() / (ye - ys).total_seconds()


def date_to_jd(y, m, d):
    dt = datetime(y, m, d, tzinfo=timezone.utc)
    return (dt - datetime(1970, 1, 1, tzinfo=timezone.utc)).total_seconds() / 86400 + JD_UNIX_EPOCH


# ======================================================================
# AAVSO download
# ======================================================================

def fetch_chunk(api_key, start_date, end_date, verify_tls):
    url = (f"{API_BASE}?target={TARGET.replace(' ', '+')}"
           f"&start_date={start_date}&end_date={end_date}&band={BAND}"
           f"&observer=&obs_campaign=&submit=Search&output_format=csv"
           f"&include_comparison_fields=false")
    resp = requests.get(url, headers={'Authorization': f'Token {api_key}'},
                         timeout=90, verify=verify_tls)
    resp.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
        csv_name = next(n for n in zf.namelist() if n.endswith('.csv'))
        with zf.open(csv_name) as f:
            text = io.TextIOWrapper(f, encoding='utf-8')
            return list(csv.DictReader(text))


def fetch_range(api_key, start_date, end_date, verify_tls, chunk_days=365):
    """Fetch start_date..end_date, splitting into <=chunk_days chunks
    (a single multi-year request can 502 on AAVSO's end for a busy target)."""
    start = datetime.strptime(start_date, '%Y-%m-%d')
    end = datetime.strptime(end_date, '%Y-%m-%d')
    if start > end:
        return []
    rows = []
    cur = start
    while cur <= end:
        chunk_end = min(cur + timedelta(days=chunk_days - 1), end)
        s, e = cur.strftime('%Y-%m-%d'), chunk_end.strftime('%Y-%m-%d')
        print(f"  fetching {s} .. {e} ...", end=' ', flush=True)
        try:
            chunk_rows = fetch_chunk(api_key, s, e, verify_tls)
            print(f"{len(chunk_rows)} rows")
            rows.extend(chunk_rows)
        except Exception as exc:
            print(f"FAILED ({exc})")
            raise
        cur = chunk_end + timedelta(days=1)
    return rows


# ======================================================================
# CSV I/O
# ======================================================================

def load_csv_rows(path):
    """Returns list of dict rows (raw, unfiltered) and the max JD present (or None)."""
    rows = []
    max_jd = None
    try:
        with open(path, newline='') as f:
            r = csv.DictReader(f)
            for row in r:
                rows.append(row)
                try:
                    jd = float(row['jd'])
                    if max_jd is None or jd > max_jd:
                        max_jd = jd
                except (ValueError, KeyError):
                    pass
    except FileNotFoundError:
        pass
    return rows, max_jd


def write_csv_rows(path, rows):
    """Writes rows sorted by jd. '#' is renumbered sequentially on every write —
    it's cosmetic only (see DEDUPE_FIELDS above), never used as an identity key."""
    rows_sorted = sorted(rows, key=lambda r: float(r['jd']))
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for i, row in enumerate(rows_sorted, 1):
            out = {k: row.get(k, '') for k in CSV_FIELDS}
            out['#'] = i
            w.writerow(out)


DEDUPE_FIELDS = [f for f in CSV_FIELDS if f != '#']


def _content_key(row):
    # NOTE: '#' is deliberately excluded — it's only a per-request row
    # counter (restarts at 1 for every separate download), not a stable
    # observation id, so it can't be used to detect duplicates across runs.
    return tuple(row.get(k, '') for k in DEDUPE_FIELDS)


def merge_rows(existing_rows, new_rows):
    """Dedupe by full row content (everything except '#')."""
    by_key = {_content_key(row): row for row in existing_rows}
    added = 0
    for row in new_rows:
        key = _content_key(row)
        if key not in by_key:
            added += 1
        by_key[key] = row
    return list(by_key.values()), added


def usable_points(rows):
    """(jd, mag) pairs, Johnson V only, excluding fainter-than (non-detection) rows."""
    pts = []
    for row in rows:
        if row.get('fainterthan') == 'True':
            continue
        try:
            jd = float(row['jd'])
            mag = float(row['mag'])
        except (ValueError, KeyError):
            continue
        pts.append((jd, mag))
    pts.sort()
    return pts


# ======================================================================
# Nightly binning
# ======================================================================

def nightly_bin(points, jd_start=None, roll_window=7.0):
    """points: list of (jd, mag). Returns list of dicts sorted by jd:
    {jd, mag, sd, n, roll, date} — mag/sd are the nightly mean/stdev.
    Observations are grouped into a night by round(jd), but the reported
    `jd` is the mean of the *actual* jds in that group (sub-day precision
    preserved — matters for the epoch quadratic fit), not the rounded key.
    `roll` is a `roll_window`-day centered rolling mean."""
    nightly = defaultdict(list)
    for jd, mag in points:
        if jd_start is not None and jd < jd_start:
            continue
        nightly[round(jd)].append((jd, mag))

    nights = []
    for key, obs in sorted(nightly.items()):
        jds = [o[0] for o in obs]
        mags = [o[1] for o in obs]
        n = len(mags)
        jd_mean = sum(jds) / n
        mean = sum(mags) / n
        sd = (sum((m - mean) ** 2 for m in mags) / n) ** 0.5 if n > 1 else 0.0
        nights.append({'jd': round(jd_mean, 5), 'mag': round(mean, 4), 'sd': round(sd, 4), 'n': n})

    for nt in nights:
        c = nt['jd']
        win = [x['mag'] for x in nights if abs(x['jd'] - c) <= roll_window]
        nt['roll'] = round(sum(win) / len(win), 4)
        nt['date'] = jd_to_date_str(c)

    return nights


# ======================================================================
# Weighted quadratic fit (with covariance, for error propagation)
# ======================================================================

def _det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def _inv3(m):
    D = _det3(m)
    if abs(D) < 1e-12:
        return None

    def minor(mm, i, j):
        rows_ = [r for k, r in enumerate(mm) if k != i]
        return [[v for k, v in enumerate(r) if k != j] for r in rows_]

    def det2(m2):
        return m2[0][0] * m2[1][1] - m2[0][1] * m2[1][0]

    cof = [[((-1) ** (i + j)) * det2(minor(m, i, j)) for j in range(3)] for i in range(3)]
    return [[cof[j][i] / D for j in range(3)] for i in range(3)]


def weighted_quadratic_fit(t, y, w):
    """Least squares y = a*t^2 + b*t + c, plus covariance of (a,b) for error propagation."""
    Sw = sum(w)
    Swt = sum(wi * ti for wi, ti in zip(w, t))
    Swt2 = sum(wi * ti * ti for wi, ti in zip(w, t))
    Swt3 = sum(wi * ti ** 3 for wi, ti in zip(w, t))
    Swt4 = sum(wi * ti ** 4 for wi, ti in zip(w, t))
    Swy = sum(wi * yi for wi, yi in zip(w, y))
    Swty = sum(wi * ti * yi for wi, ti, yi in zip(w, t, y))
    Swt2y = sum(wi * ti * ti * yi for wi, ti, yi in zip(w, t, y))
    A = [[Swt4, Swt3, Swt2], [Swt3, Swt2, Swt], [Swt2, Swt, Sw]]
    B = [Swt2y, Swty, Swy]
    D = _det3(A)
    if abs(D) < 1e-12:
        return None
    sols = []
    for col in range(3):
        Ai = [row[:] for row in A]
        for r_ in range(3):
            Ai[r_][col] = B[r_]
        sols.append(_det3(Ai) / D)
    a, b, c = sols

    Ainv = _inv3(A)
    resid = [yi - (a * ti * ti + b * ti + c) for ti, yi in zip(t, y)]
    n = len(t)
    if n <= 3 or Ainv is None:
        return a, b, c, None
    ssr = sum(wi * ri * ri for wi, ri in zip(w, resid))
    sigma2 = ssr / (n - 3)
    cov = {'var_a': Ainv[0][0] * sigma2, 'var_b': Ainv[1][1] * sigma2, 'cov_ab': Ainv[0][1] * sigma2}
    return a, b, c, cov


# ======================================================================
# Epoch fitting
# ======================================================================

def fit_epoch(nights, t_calc, half_window=50.0, min_nights=10):
    """Try to fit a single primary-minimum epoch. Returns (result_dict_or_None, diag_dict).
    diag_dict always has n/before/after so callers can explain a rejection."""
    window = [(nt['jd'], nt['mag']) for nt in nights if abs(nt['jd'] - t_calc) <= half_window]
    before = [w for w in window if w[0] < t_calc]
    after = [w for w in window if w[0] >= t_calc]
    diag = {'n': len(window), 'before': len(before), 'after': len(after)}

    if len(window) < min_nights or len(before) < 3 or len(after) < 3:
        return None, diag

    pts = window[:]
    fit = None
    for _ in range(3):
        t = [p[0] - t_calc for p in pts]
        y = [p[1] for p in pts]
        w = [1.0] * len(pts)
        fit = weighted_quadratic_fit(t, y, w)
        if fit is None:
            return None, diag
        a, b, c, cov = fit
        resid = [yi - (a * ti * ti + b * ti + c) for ti, yi in zip(t, y)]
        sd = (sum(r * r for r in resid) / len(resid)) ** 0.5
        if sd == 0:
            break
        new_pts = [p for p, r in zip(pts, resid) if abs(r) <= 3 * sd]
        if len(new_pts) == len(pts) or len(new_pts) < min_nights:
            pts = new_pts if len(new_pts) >= min_nights else pts
            break
        pts = new_pts

    if fit is None:
        return None, diag
    a, b, c, cov = fit
    if a >= 0:
        return None, diag  # concave-up -> not a real minimum in this window
    t_vertex_offset = -b / (2 * a)
    if abs(t_vertex_offset) > half_window:
        return None, diag

    t_obs = t_calc + t_vertex_offset
    oc = t_obs - t_calc
    mag_min = a * t_vertex_offset ** 2 + b * t_vertex_offset + c

    t_obs_err = None
    if cov is not None:
        dtv_da = b / (2 * a * a)
        dtv_db = -1.0 / (2 * a)
        var_tv = (dtv_da ** 2) * cov['var_a'] + (dtv_db ** 2) * cov['var_b'] + 2 * dtv_da * dtv_db * cov['cov_ab']
        if var_tv > 0:
            t_obs_err = var_tv ** 0.5

    result = {
        'E': None,  # filled in by caller
        't_calc': t_calc, 't_obs': t_obs, 'OC_days': oc, 'OC_err_days': t_obs_err,
        'n_nights': len(pts), 'year': jd_to_decyear(t_obs), 'mag_min': mag_min,
    }
    return result, diag


def scan_epochs(nights, P, T0, half_window=50.0, min_nights=10):
    """Scan every epoch spanned by the data. Returns (fitted, gaps) lists,
    each entry a dict; fitted entries have status='fit', gaps have status='gap'."""
    jd_min, jd_max = nights[0]['jd'], nights[-1]['jd']
    E_min = math.floor((jd_min - T0) / P) - 1
    E_max = math.ceil((jd_max - T0) / P) + 1

    fitted, gaps = [], []
    for E in range(E_min, E_max + 1):
        t_calc = T0 + E * P
        result, diag = fit_epoch(nights, t_calc, half_window, min_nights)
        if result is not None:
            result['E'] = E
            result['status'] = 'fit'
            fitted.append(result)
        else:
            gaps.append({'E': E, 't_calc': t_calc, 'year': jd_to_decyear(t_calc),
                          'status': 'gap', **diag})
    return fitted, gaps


def weighted_linear_fit(fitted):
    """OC_days = dT0 + dP*E, weighted by 1/err^2. Returns (dP, dT0) or (None, None)."""
    if len(fitted) < 2:
        return None, None
    we = [1.0 / (r['OC_err_days'] ** 2) if r['OC_err_days'] else 1.0 for r in fitted]
    Ee = [r['E'] for r in fitted]
    OCe = [r['OC_days'] for r in fitted]
    Sw = sum(we)
    SwE = sum(wi * Ei for wi, Ei in zip(we, Ee))
    SwEE = sum(wi * Ei * Ei for wi, Ei in zip(we, Ee))
    SwOC = sum(wi * oi for wi, oi in zip(we, OCe))
    SwEOC = sum(wi * Ei * oi for wi, Ei, oi in zip(we, Ee, OCe))
    denom = Sw * SwEE - SwE * SwE
    if abs(denom) < 1e-9:
        return None, None
    dP = (Sw * SwEOC - SwE * SwOC) / denom
    dT0 = (SwEE * SwOC - SwE * SwEOC) / denom
    return dP, dT0


def build_pending_and_gaps(gaps, jd_max_data, half_window):
    """Split 'gap' epochs into true historical gaps (window fully elapsed, just
    not enough data) vs. pending (the window extends past what we've downloaded
    yet, i.e. the minimum hasn't fully happened/been observed)."""
    out = []
    for g in gaps:
        g = dict(g)
        if g['t_calc'] + half_window > jd_max_data:
            g['status'] = 'pending'
        out.append(g)
    return out


# ======================================================================
# Main
# ======================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--csv', default=DEFAULT_CSV, help='Path to the master photometry CSV')
    ap.add_argument('--out-dir', default=DEFAULT_OUT_DIR, help='Where to write oc_data.json / lightcurve_data.json')
    ap.add_argument('--api-key', default=os.environ.get('AAVSO_API_KEY'),
                     help='AAVSO API key (or set AAVSO_API_KEY env var)')
    ap.add_argument('--end-date', default=None, help='Fetch through this date (default: today, UTC)')
    ap.add_argument('--period', type=float, default=P_DEFAULT)
    ap.add_argument('--t0', type=float, default=T0_DEFAULT)
    ap.add_argument('--skip-tls-verify', action='store_true',
                     help='Disable TLS verification (only for networks with TLS-intercepting proxies)')
    ap.add_argument('--no-fetch', action='store_true', help='Skip the AAVSO download; just re-fit existing CSV')
    args = ap.parse_args()

    csv_path = os.path.abspath(args.csv)
    out_dir = os.path.abspath(args.out_dir)
    verify_tls = not args.skip_tls_verify
    end_date = args.end_date or datetime.now(timezone.utc).strftime('%Y-%m-%d')

    print(f"CSV:     {csv_path}")
    print(f"Out dir: {out_dir}")

    existing_rows, max_jd = load_csv_rows(csv_path)
    print(f"Existing rows: {len(existing_rows)}"
          + (f" (through {jd_to_date_str(max_jd)})" if max_jd else " (no existing file)"))

    added = 0
    if args.no_fetch:
        print("--no-fetch: skipping AAVSO download.")
    elif not args.api_key:
        print("WARNING: no AAVSO API key (set AAVSO_API_KEY or pass --api-key) — "
              "skipping download, only re-fitting existing data.")
    else:
        if max_jd:
            start_dt = jd_to_datetime(max_jd) - timedelta(days=1)  # 1-day overlap, dedup handles it
        else:
            start_dt = datetime(1946, 1, 1, tzinfo=timezone.utc)
        start_date = start_dt.strftime('%Y-%m-%d')
        print(f"Fetching AAVSO observations {start_date} .. {end_date} ...")
        new_rows = fetch_range(args.api_key, start_date, end_date, verify_tls)
        existing_rows, added = merge_rows(existing_rows, new_rows)
        write_csv_rows(csv_path, existing_rows)
        print(f"Added {added} new observation(s). Total rows now: {len(existing_rows)}.")

    # ---- re-fit everything from the (possibly updated) CSV ----
    # Fitting uses the FULL history (unrestricted) so a window near the display
    # boundary can still use real pre-2005 nights; only the light-curve JSON's
    # displayed points get clipped to LIGHTCURVE_START.
    points = usable_points(existing_rows)
    all_nights = nightly_bin(points)
    if not all_nights:
        print("No usable nightly data found — aborting.")
        sys.exit(1)
    jd_lc_start = date_to_jd(*LIGHTCURVE_START)
    nights = [n for n in all_nights if n['jd'] >= jd_lc_start]

    fitted, gaps = scan_epochs(all_nights, args.period, args.t0, HALF_WINDOW, MIN_NIGHTS)
    dP, dT0 = weighted_linear_fit(fitted)
    jd_max_data = all_nights[-1]['jd']
    gaps = build_pending_and_gaps(gaps, jd_max_data, HALF_WINDOW)

    # Only keep gap/pending epochs whose window overlaps the display range
    # (E far before 2005 or many cycles past "today" aren't meaningful to show).
    gaps = [g for g in gaps if jd_lc_start - HALF_WINDOW <= g['t_calc'] <= jd_max_data + args.period]

    generated_at = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    data_through = jd_to_date_str(jd_max_data)
    params = {'P': args.period, 'T0': args.t0, 'dP': dP, 'dT0': dT0}

    # ---- diff against the previous oc_data.json, if any, for the summary ----
    oc_path = os.path.join(out_dir, 'oc_data.json')
    prev_fitted_epochs = set()
    if os.path.exists(oc_path):
        try:
            with open(oc_path) as f:
                prev = json.load(f)
            prev_fitted_epochs = {r['E'] for r in prev.get('results', [])}
        except (json.JSONDecodeError, KeyError):
            pass
    new_fitted_epochs = {r['E'] for r in fitted} - prev_fitted_epochs

    # ---- write oc_data.json ----
    with open(oc_path, 'w') as f:
        json.dump({
            'generated_at': generated_at,
            'data_through': data_through,
            'params': params,
            'results': fitted,
            'missing': gaps,
        }, f, indent=2)
    print(f"Wrote {oc_path}")

    # ---- write lightcurve_data.json ----
    minima = sorted(fitted + gaps, key=lambda r: r['E'])
    lc_path = os.path.join(out_dir, 'lightcurve_data.json')
    with open(lc_path, 'w') as f:
        json.dump({
            'generated_at': generated_at,
            'data_through': data_through,
            'params': {'P': args.period, 'T0': args.t0},
            'nights': nights,
            'minima': minima,
        }, f, separators=(',', ':'))
    print(f"Wrote {lc_path}")

    # ---- summary ----
    n_pending = sum(1 for g in gaps if g['status'] == 'pending')
    n_gap = sum(1 for g in gaps if g['status'] == 'gap')
    print("\n--- Summary ---")
    print(f"New observations added:   {added}")
    print(f"Nights in light curve:    {len(nights)} (from {LIGHTCURVE_START[0]}-01-01 to {data_through})")
    print(f"Fitted ELL minima:        {len(fitted)}")
    print(f"  of which newly fitted since last run: {sorted(new_fitted_epochs) if new_fitted_epochs else 'none'}")
    print(f"Historical coverage gaps: {n_gap}")
    print(f"Pending (not yet enough data): {n_pending}")
    if dP is not None:
        print(f"Weighted O-C trend: {dT0:+.4f} + ({dP:+.6f})*E days  ->  corrected P = {args.period+dP:.6f} d")


if __name__ == '__main__':
    main()
