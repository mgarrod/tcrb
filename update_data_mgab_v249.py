#!/usr/bin/env python3
"""
update_data_mgab_v249.py — eclipse-timing (O-C) refresh for MGAB-V249.

MGAB-V249 is an ultra-compact eclipsing binary (hot subdwarf + white dwarf,
VSX type E/WD+IBWD; Murawski / Kupfer et al.) with P = 0.0273195 d (~39.3 min)
and VSX epoch T0 = HJD 2456716.9169. That period is ~4 orders of magnitude
shorter than T CrB's 227.5 d ellipsoidal period, so the T CrB pipeline's
"one point per night" nightly-binning approach doesn't apply here — a single
night spans hundreds of orbital cycles, not a fraction of one.

What this script does instead, each run:
  1. Reads the existing photometry CSV and finds the latest date already
     present; downloads only AAVSO Johnson V observations newer than that
     (same incremental-merge approach as update_data.py).
  2. Groups observations into observing runs (gaps > RUN_GAP_DAYS apart).
  3. Refines a self-consistent eclipse-minimum reference epoch T0_ref by
     phase-wrapping ALL points onto the nearest predicted minimum and
     iterating a weighted quadratic fit (the VSX epoch's phase zero isn't
     at the deep minimum, so this re-centers on the actual eclipse).
  4. For each observing run, phase-wraps that run's own points onto its
     nearest predicted eclipse and quadratic-fits the minimum — one
     observed-minus-calculated timing residual per run (not per individual
     ~39-minute cycle; a single cycle only gets 1-2 points at AAVSO cadence).
  5. Writes oc_data_mgab_v249.json/.js (O-C results) consumed by
     mgab_v249_oc.html.

Requires an AAVSO API key (Account Settings -> API Key on aavso.org), passed
via AAVSO_API_KEY or --api-key. Without one, it just re-fits whatever is
already in the CSV.

Usage:
    AAVSO_API_KEY=xxxx python3 update_data_mgab_v249.py
    python3 update_data_mgab_v249.py --no-fetch
"""
import argparse
import csv
import io
import json
import math
import os
import sys
import zipfile
from datetime import datetime, timedelta, timezone

import requests

# ======================================================================
# Config
# ======================================================================

DEFAULT_CSV = os.path.join(os.path.dirname(__file__), '..', 'MGAB-V249_V.csv')
DEFAULT_OUT_DIR = os.path.dirname(__file__)

TARGET = 'MGAB-V249'
BAND = 2  # Johnson V
P_DEFAULT = 0.027319513001          # VSX, Kupfer et al.
T0_VSX_DEFAULT = 2456716.9169       # VSX epoch (phase reference, NOT at the deep minimum)
API_BASE = 'https://apps.aavso.org/v2/api/observations/photometry/download/'

RUN_GAP_DAYS = 0.3        # points more than this apart split into separate observing runs
HALF_WINDOW_FRAC = 0.27   # fit window = +/- this fraction of P around predicted minimum
MIN_POINTS = 6
MIN_BEFORE_AFTER = 3
REFINE_ITERS = 5
COARSE_BINS = 100         # for the initial phase-of-minimum guess

JD_UNIX_EPOCH = 2440587.5

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


# ======================================================================
# AAVSO download (same shape as update_data.py's fetch_chunk/fetch_range)
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
# CSV I/O (identical scheme to update_data.py)
# ======================================================================

def load_csv_rows(path):
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
    return tuple(row.get(k, '') for k in DEDUPE_FIELDS)


def merge_rows(existing_rows, new_rows):
    by_key = {_content_key(row): row for row in existing_rows}
    added = 0
    for row in new_rows:
        key = _content_key(row)
        if key not in by_key:
            added += 1
        by_key[key] = row
    return list(by_key.values()), added


def usable_points(rows):
    """(jd, mag, observer) triples, Johnson V only, excluding fainter-than rows."""
    pts = []
    for row in rows:
        if row.get('fainterthan') == 'True':
            continue
        try:
            jd = float(row['jd'])
            mag = float(row['mag'])
        except (ValueError, KeyError):
            continue
        pts.append((jd, mag, row.get('observer', '')))
    pts.sort()
    return pts


# ======================================================================
# Observing-run grouping
# ======================================================================

def group_runs(points, gap_days=RUN_GAP_DAYS):
    """points: sorted (jd, mag, observer). Splits into runs wherever the
    gap between consecutive points exceeds gap_days."""
    if not points:
        return []
    runs = [[points[0]]]
    for p in points[1:]:
        if p[0] - runs[-1][-1][0] > gap_days:
            runs.append([p])
        else:
            runs[-1].append(p)
    return runs


# ======================================================================
# Weighted quadratic fit (with covariance) — same algebra as update_data.py
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


def _vertex_with_err(a, b, cov):
    v = -b / (2 * a)
    err = None
    if cov is not None:
        dtv_da = b / (2 * a * a)
        dtv_db = -1.0 / (2 * a)
        var = (dtv_da ** 2) * cov['var_a'] + (dtv_db ** 2) * cov['var_b'] + 2 * dtv_da * dtv_db * cov['cov_ab']
        if var > 0:
            err = var ** 0.5
    return v, err


# ======================================================================
# Phase-wrapped minimum fitting (the MGAB-V249-specific part)
# ======================================================================

def wrap_dt(jd, center, P):
    """Signed time offset (days) from `center` to the nearest multiple of P."""
    return ((jd - center + P / 2) % P) - P / 2


def fit_minimum(points, t_calc, P, half_window_frac=HALF_WINDOW_FRAC,
                 min_points=MIN_POINTS, min_before_after=MIN_BEFORE_AFTER):
    """Phase-wraps `points` (jd, mag, ...) onto t_calc and quadratic-fits the
    minimum. `points` may span many cycles of P — all of them get folded onto
    the single nearest predicted minimum at t_calc. Returns a result dict or
    None, plus a diagnostic dict."""
    half = half_window_frac * P
    t, y = [], []
    for p in points:
        jd, mag = p[0], p[1]
        dt = wrap_dt(jd, t_calc, P)
        if abs(dt) <= half:
            t.append(dt)
            y.append(mag)
    before = sum(1 for x in t if x < 0)
    after = sum(1 for x in t if x >= 0)
    diag = {'n': len(t), 'before': before, 'after': after}

    if len(t) < min_points or before < min_before_after or after < min_before_after:
        return None, diag

    w = [1.0] * len(t)
    fit = weighted_quadratic_fit(t, y, w)
    if fit is None:
        return None, diag
    a, b, c, cov = fit
    if a >= 0:
        return None, diag  # concave-up: not a real minimum in this window
    v, v_err = _vertex_with_err(a, b, cov)
    if abs(v) > half:
        return None, diag

    t_obs = t_calc + v
    result = {
        't_calc': t_calc, 't_obs': t_obs,
        'OC_days': v, 'OC_err_days': v_err,
        'n_points': len(t), 'mag_min': a * v * v + b * v + c,
        'year': jd_to_decyear(t_obs),
    }
    return result, diag


def refine_reference_epoch(points, P, T0_vsx, half_window_frac=HALF_WINDOW_FRAC,
                            coarse_bins=COARSE_BINS, iters=REFINE_ITERS):
    """VSX's phase-zero isn't necessarily at the deep eclipse minimum. Finds
    a self-consistent T0_ref near the actual minimum, anchored to the cycle
    nearest the first data point, by: (1) coarse-binning phase vs mag to get
    an approximate minimum phase, then (2) iterating the pooled quadratic fit
    (all points, folded onto that one cycle) to convergence."""
    bins = {}
    for jd, mag, _obs in points:
        phase = ((jd - T0_vsx) / P) % 1.0
        b = min(int(phase * coarse_bins), coarse_bins - 1)
        bins.setdefault(b, []).append(mag)
    best_bin = max(bins, key=lambda b: sum(bins[b]) / len(bins[b]))
    phi0_approx = (best_bin + 0.5) / coarse_bins

    jd_first = points[0][0]
    E0 = round((jd_first - T0_vsx) / P)
    t0_ref = T0_vsx + E0 * P + phi0_approx * P

    for _ in range(iters):
        result, diag = fit_minimum(points, t0_ref, P, half_window_frac,
                                    min_points=MIN_POINTS, min_before_after=MIN_BEFORE_AFTER)
        if result is None:
            return None, diag
        new_ref = result['t_obs']
        if abs(new_ref - t0_ref) < 1e-7:
            t0_ref = new_ref
            break
        t0_ref = new_ref
    return t0_ref, None


def weighted_linear_fit(fitted, key_oc='OC_days', key_err='OC_err_days'):
    """OC = dT0 + dP*E, weighted by 1/err^2. Returns (dP, dT0) or (None, None)."""
    if len(fitted) < 2:
        return None, None
    we = [1.0 / (r[key_err] ** 2) if r.get(key_err) else 1.0 for r in fitted]
    Ee = [r['E'] for r in fitted]
    OCe = [r[key_oc] for r in fitted]
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


def write_json_and_js(path_no_ext, varname, obj, **json_kwargs):
    payload = json.dumps(obj, **json_kwargs)
    with open(path_no_ext + '.json', 'w') as f:
        f.write(payload)
    with open(path_no_ext + '.js', 'w') as f:
        f.write(f'const {varname} = {payload};\n')


# ======================================================================
# Main
# ======================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--csv', default=DEFAULT_CSV, help='Path to the master photometry CSV')
    ap.add_argument('--out-dir', default=DEFAULT_OUT_DIR, help='Where to write oc_data_mgab_v249.json/.js')
    ap.add_argument('--api-key', default=os.environ.get('AAVSO_API_KEY'),
                     help='AAVSO API key (or set AAVSO_API_KEY env var)')
    ap.add_argument('--end-date', default=None, help='Fetch through this date (default: today, UTC)')
    ap.add_argument('--period', type=float, default=P_DEFAULT)
    ap.add_argument('--t0-vsx', type=float, default=T0_VSX_DEFAULT,
                     help='VSX phase-reference epoch (not assumed to be at the minimum)')
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
            start_dt = jd_to_datetime(max_jd) - timedelta(days=1)
        else:
            start_dt = datetime(2026, 1, 1, tzinfo=timezone.utc)  # MGAB-V249 discovered 2026
        start_date = start_dt.strftime('%Y-%m-%d')
        print(f"Fetching AAVSO observations {start_date} .. {end_date} ...")
        new_rows = fetch_range(args.api_key, start_date, end_date, verify_tls)
        existing_rows, added = merge_rows(existing_rows, new_rows)
        write_csv_rows(csv_path, existing_rows)
        print(f"Added {added} new observation(s). Total rows now: {len(existing_rows)}.")

    points = usable_points(existing_rows)
    if not points:
        print("No usable V-band data found — aborting.")
        sys.exit(1)

    runs = group_runs(points)
    print(f"Observing runs: {len(runs)}")

    t0_ref, fail_diag = refine_reference_epoch(points, args.period, args.t0_vsx)
    if t0_ref is None:
        print(f"Could not refine a reference eclipse epoch (diag={fail_diag}) — aborting.")
        sys.exit(1)
    print(f"Refined eclipse-minimum reference epoch T0_ref = {t0_ref:.6f} HJD "
          f"(VSX phase-zero T0 = {args.t0_vsx})")

    fitted, missing = [], []
    for run in runs:
        run_mid = sum(p[0] for p in run) / len(run)
        E = round((run_mid - t0_ref) / args.period)
        t_calc = t0_ref + E * args.period
        result, diag = fit_minimum(run, t_calc, args.period)
        run_start, run_end = run[0][0], run[-1][0]
        observers = sorted({p[2] for p in run if p[2]})
        if result is not None:
            result.update({
                'E': E, 'run_start': jd_to_date_str(run_start), 'run_end': jd_to_date_str(run_end),
                'duration_hr': (run_end - run_start) * 24, 'observers': observers,
            })
            result['OC_sec'] = result['OC_days'] * 86400
            result['OC_err_sec'] = result['OC_err_days'] * 86400 if result['OC_err_days'] else None
            fitted.append(result)
        else:
            missing.append({
                'E': E, 't_calc': t_calc, 'year': jd_to_decyear(t_calc),
                'run_start': jd_to_date_str(run_start), 'run_end': jd_to_date_str(run_end),
                'observers': observers, **diag,
            })

    fitted.sort(key=lambda r: r['E'])
    dP_sec, dT0_sec = weighted_linear_fit(fitted, key_oc='OC_sec', key_err='OC_err_sec')

    # phase-folded light curve (full dataset, for a supplementary panel).
    # Centered on the fitted eclipse minimum (phase 0, range [-0.5, 0.5))
    # rather than VSX phase-zero, so the eclipse sits mid-panel, not split
    # across both edges.
    folded = [{'phase': round(wrap_dt(jd, t0_ref, args.period) / args.period, 5), 'mag': mag}
              for jd, mag, _obs in points]

    generated_at = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    data_through = jd_to_date_str(points[-1][0])
    all_observers = sorted({p[2] for p in points if p[2]})

    oc_base = os.path.join(out_dir, 'oc_data_mgab_v249')
    write_json_and_js(oc_base, 'MGAB_V249_OC_DATA', {
        'generated_at': generated_at,
        'target': TARGET,
        'data_through': data_through,
        'observers': all_observers,
        'params': {
            'P': args.period, 'T0_vsx': args.t0_vsx, 'T0_ref': t0_ref,
            'dP_sec': dP_sec, 'dT0_sec': dT0_sec, 'units': 'seconds',
        },
        'results': fitted,
        'missing': missing,
        'folded': folded,
    }, indent=2)
    print(f"Wrote {oc_base}.json and .js")

    print("\n--- Summary ---")
    print(f"New observations added:   {added}")
    print(f"Observing runs:            {len(runs)}")
    print(f"Fitted eclipse timings:    {len(fitted)}")
    print(f"Runs with insufficient data: {len(missing)}")
    if dP_sec is not None:
        print(f"Weighted O-C trend: {dT0_sec:+.2f} + ({dP_sec:+.4f})*E seconds")
    for r in fitted:
        print(f"  E={r['E']:>6d}  {r['run_start']}  O-C = {r['OC_sec']:+.2f}"
              f" ± {r['OC_err_sec']:.2f} s  (n={r['n_points']})")


if __name__ == '__main__':
    main()
