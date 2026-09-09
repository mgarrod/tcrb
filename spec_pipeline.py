#!/usr/bin/env python3
"""
spec_pipeline.py — extract, normalize, and analyze the T CrB spectral time series.

Reads every FITS spectrum in spec/AVSpec_Download/, resamples onto a common
wavelength grid (the overlap region shared by all 139 files), continuum-
normalizes each one, measures the standard symbiotic-star diagnostic lines
(Hβ, He I/Na D, Hα), and computes each spectrum's orbital phase against the
same ELL ephemeris used for the photometric O-C analysis.

Outputs (into spec/analysis/):
  spectra_data.json    -> wavelength grid + normalized flux matrix + per-epoch
                           metadata, for the timelapse HTML page
  spectra_summary.csv  -> one row per spectrum: JD, date, phase, line metrics
  variability_spectrum.json -> pixel-wise std-dev across all epochs (data-driven
                           "what's changing" spectrum, no line list assumed)

Usage:
    python3 spec_pipeline.py
"""
import glob
import json
import os
import re
from datetime import datetime, timedelta, timezone

import numpy as np
from astropy.io import fits

FITS_DIR = os.path.join(os.path.dirname(__file__), '..', 'AVSpec_Download')
OUT_DIR = os.path.dirname(__file__)
PHOTOMETRY_CSV = os.path.join(os.path.dirname(__file__), '..', '..', 'TCrB_V_1946-2026.csv')

# Same ephemeris as the O-C / light-curve pipeline (VSX oid=10602)
P = 227.5528
T0 = 2455828.9

# Common grid: the wavelength range shared by every spectrum (4859.4-6771.6 A
# per the survey), trimmed slightly for a clean round-number grid.
GRID_START, GRID_END, GRID_STEP = 4870.0, 6760.0, 1.0

# Continuum: wide running-median window (odd, in pixels at GRID_STEP=1A/px)
CONTINUUM_WINDOW = 151

# Diagnostic lines: (name, center, half-width for peak search, [blue_cont, red_cont] for EW)
LINES = {
    'Hbeta': {'center': 4861.3, 'peak_hw': 12, 'ew_window': 25},
    'HeI_NaD': {'center': 5884.0, 'peak_hw': 14, 'ew_window': 20},  # He I 5876 / Na D 5890/5896 blend
    'Halpha': {'center': 6562.8, 'peak_hw': 15, 'ew_window': 30},
}

FNAME_JD_RE = re.compile(r'_(\d+\.\d+)\.FITS$', re.IGNORECASE)


def write_json_and_js(path_no_ext, varname, obj, **json_kwargs):
    """Writes both <path>.json (plain, for any tool/interop) and <path>.js
    (`const VARNAME = {...};`, loadable via <script src> so the HTML pages
    work when opened directly via file:// — fetch() of local files is
    blocked by browsers, but <script src> of local files is not)."""
    payload = json.dumps(obj, **json_kwargs)
    with open(path_no_ext + '.json', 'w') as f:
        f.write(payload)
    with open(path_no_ext + '.js', 'w') as f:
        f.write(f'const {varname} = {payload};\n')


def jd_to_date_str(jd):
    return (datetime(1970, 1, 1, tzinfo=timezone.utc) + timedelta(days=jd - 2440587.5)).strftime('%Y-%m-%d')


def jd_to_decyear(jd):
    d = datetime(1970, 1, 1, tzinfo=timezone.utc) + timedelta(days=jd - 2440587.5)
    ys = datetime(d.year, 1, 1, tzinfo=timezone.utc)
    ye = datetime(d.year + 1, 1, 1, tzinfo=timezone.utc)
    return d.year + (d - ys).total_seconds() / (ye - ys).total_seconds()


def load_spectrum(path):
    with fits.open(path) as hdul:
        h = hdul[0].header
        flux = hdul[0].data.astype(float)
        crpix1 = float(h.get('CRPIX1', 1))
        crval1 = float(h.get('CRVAL1'))
        cdelt1 = float(h.get('CDELT1'))
        wl = crval1 + (np.arange(len(flux)) + 1 - crpix1) * cdelt1
        instrument = h.get('AAV_INST', h.get('INSTRUME', h.get('USPEC', 'unknown')))
        exptime = h.get('EXPTIME')
    return wl, flux, {'instrument': instrument, 'exptime': exptime}


def running_median_continuum(flux, window):
    """Continuum estimate: running median with edge padding (reflect)."""
    half = window // 2
    padded = np.pad(flux, half, mode='reflect')
    n = len(flux)
    cont = np.empty(n)
    # simple O(n*window) approach; fine for ~1900-point spectra x 139 files
    for i in range(n):
        cont[i] = np.median(padded[i:i + window])
    return cont


def measure_line(wl_grid, norm_flux, center, peak_hw, ew_window):
    peak_mask = (wl_grid >= center - peak_hw) & (wl_grid <= center + peak_hw)
    if not peak_mask.any():
        return None, None
    peak = float(np.nanmax(norm_flux[peak_mask]))
    ew_mask = (wl_grid >= center - ew_window) & (wl_grid <= center + ew_window)
    if ew_mask.sum() < 3:
        return peak, None
    excess = norm_flux[ew_mask] - 1.0
    ew = float(np.trapezoid(excess, wl_grid[ew_mask]))  # Angstrom, positive = net emission
    return peak, ew


def load_photometry(path):
    """(jd, mag) pairs, Johnson V, non-detections excluded — for nearest-epoch cross-reference."""
    import csv
    pts = []
    try:
        with open(path, newline='') as f:
            for row in csv.DictReader(f):
                if row.get('fainterthan') == 'True':
                    continue
                try:
                    pts.append((float(row['jd']), float(row['mag'])))
                except (ValueError, KeyError):
                    continue
    except FileNotFoundError:
        print(f"WARNING: photometry CSV not found at {path}, skipping V-mag cross-reference")
        return []
    pts.sort()
    return pts


def nearest_photometry(phot, jd, max_dt=1.0):
    if not phot:
        return None, None
    import bisect
    jds = [p[0] for p in phot]
    i = bisect.bisect_left(jds, jd)
    candidates = []
    if i < len(phot):
        candidates.append(phot[i])
    if i > 0:
        candidates.append(phot[i - 1])
    if not candidates:
        return None, None
    best = min(candidates, key=lambda p: abs(p[0] - jd))
    if abs(best[0] - jd) > max_dt:
        return None, None
    return best[1], best[0] - jd


def main():
    files = sorted(glob.glob(os.path.join(FITS_DIR, '*.FITS')) + glob.glob(os.path.join(FITS_DIR, '*.fits')))
    print(f"Found {len(files)} FITS files in {FITS_DIR}")

    grid = np.arange(GRID_START, GRID_END + GRID_STEP, GRID_STEP)
    phot = load_photometry(PHOTOMETRY_CSV)
    print(f"Loaded {len(phot)} photometric points for cross-reference")

    records = []
    flux_matrix = []
    skipped = []

    for path in files:
        base = os.path.basename(path)
        m = FNAME_JD_RE.search(base)
        if not m:
            skipped.append((base, 'no JD in filename'))
            continue
        jd = float(m.group(1))

        wl, flux, meta = load_spectrum(path)
        if wl.min() > GRID_START or wl.max() < GRID_END:
            skipped.append((base, f'coverage {wl.min():.0f}-{wl.max():.0f} does not span grid'))
            continue

        # resample onto common grid (linear interpolation)
        flux_r = np.interp(grid, wl, flux)

        # continuum-normalize
        cont = running_median_continuum(flux_r, CONTINUUM_WINDOW)
        cont[cont == 0] = np.nan
        norm = flux_r / cont

        # continuum color proxy: red/blue ratio of the (self-normalized-scale) smoothed continuum
        blue_c = np.nanmedian(cont[(grid >= 4900) & (grid <= 5000)])
        red_c = np.nanmedian(cont[(grid >= 6650) & (grid <= 6750)])
        color = float(red_c / blue_c) if blue_c else None

        line_metrics = {}
        for name, spec in LINES.items():
            peak, ew = measure_line(grid, norm, spec['center'], spec['peak_hw'], spec['ew_window'])
            line_metrics[f'{name}_peak'] = peak
            line_metrics[f'{name}_ew'] = ew

        phase = ((jd - T0) / P) % 1.0
        epoch = (jd - T0) / P
        vmag, dt_days = nearest_photometry(phot, jd)

        rec = {
            'file': base, 'jd': jd, 'date': jd_to_date_str(jd), 'year': jd_to_decyear(jd),
            'phase': phase, 'epoch': epoch,
            'instrument': meta['instrument'], 'exptime': meta['exptime'],
            'continuum_color': color,
            'vmag_nearest': vmag, 'vmag_dt_days': dt_days,
            **line_metrics,
        }
        records.append(rec)
        flux_matrix.append(norm)

    print(f"Processed {len(records)} spectra, skipped {len(skipped)}")
    for b, reason in skipped:
        print(f"  SKIP {b}: {reason}")

    # sort by JD
    order = np.argsort([r['jd'] for r in records])
    records = [records[i] for i in order]
    flux_matrix = np.array(flux_matrix)[order]

    # ---- variability spectrum: pixel-wise stats across all epochs ----
    mean_spec = np.nanmean(flux_matrix, axis=0)
    std_spec = np.nanstd(flux_matrix, axis=0)
    with open(os.path.join(OUT_DIR, 'variability_spectrum.json'), 'w') as f:
        json.dump({
            'wavelength': [round(float(w), 1) for w in grid],
            'mean': [round(float(v), 4) for v in mean_spec],
            'std': [round(float(v), 4) for v in std_spec],
            'n_spectra': len(records),
        }, f)
    print("Wrote variability_spectrum.json")

    # ---- summary CSV ----
    import csv as csvmod
    fields = ['file', 'jd', 'date', 'year', 'phase', 'epoch', 'instrument', 'exptime',
              'continuum_color', 'vmag_nearest', 'vmag_dt_days']
    for name in LINES:
        fields += [f'{name}_peak', f'{name}_ew']
    with open(os.path.join(OUT_DIR, 'spectra_summary.csv'), 'w', newline='') as f:
        w = csvmod.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in records:
            w.writerow({k: r.get(k) for k in fields})
    print("Wrote spectra_summary.csv")

    write_json_and_js(os.path.join(OUT_DIR, 'spectra_summary'), 'SPECTRA_SUMMARY', {
        'generated_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'ephemeris': {'P': P, 'T0': T0},
        'lines': {name: spec['center'] for name, spec in LINES.items()},
        'records': [{k: r.get(k) for k in fields} for r in records],
    }, separators=(',', ':'))
    print("Wrote spectra_summary.json and .js")

    # ---- main data file for the timelapse HTML (fetched at load time) ----
    # flux stored per-epoch at reduced precision to keep file size reasonable
    epochs_out = []
    for r, flux_row in zip(records, flux_matrix):
        epochs_out.append({
            **{k: r[k] for k in fields},
            'flux': [None if np.isnan(v) else round(float(v), 3) for v in flux_row],
        })
    write_json_and_js(os.path.join(OUT_DIR, 'spectra_data'), 'SPECTRA_DATA', {
        'generated_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'wavelength_grid': {'start': GRID_START, 'end': GRID_END, 'step': GRID_STEP, 'n': len(grid)},
        'ephemeris': {'P': P, 'T0': T0},
        'lines': {name: spec['center'] for name, spec in LINES.items()},
        'n_spectra': len(records),
        'epochs': epochs_out,
    }, separators=(',', ':'))
    size_mb = (os.path.getsize(os.path.join(OUT_DIR, 'spectra_data.json'))
               + os.path.getsize(os.path.join(OUT_DIR, 'spectra_data.js'))) / 1e6
    print(f"Wrote spectra_data.json and .js ({size_mb:.1f} MB combined)")


if __name__ == '__main__':
    main()
