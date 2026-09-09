# T CrB ELL O–C Pipeline

Photometry pipeline and chart pages for tracking T Coronae Borealis's
ellipsoidal (ELL) variation — AAVSO Johnson V data, weighted quadratic
minimum-fitting, and two self-updating HTML charts.

## File inventory

**6 files required** across two directories. Everything else in this repo is
either a generated/deployable output of these, or a superseded earlier draft.

| File | Location | Role |
|---|---|---|
| `TCrB_V_1946-2026.csv` | `tcrb/` | Master photometry data — `update_data.py` reads/updates this |
| `update_data.py` | `tcrb/analysis/` | The one pipeline script |
| `oc_data.json` | `tcrb/analysis/` | Generated output — fetched by the O–C page |
| `lightcurve_data.json` | `tcrb/analysis/` | Generated output — fetched by the light curve page |
| `tcrb_ell_oc_coverage.html` | `tcrb/analysis/` | O–C diagram page |
| `tcrb_2005_present_lightcurve.html` | `tcrb/analysis/` | Observed/calculated light curve page |

## Deploying

Copy the two `.html` files plus `oc_data.json` and `lightcurve_data.json` to
the same directory on your web server (served over http(s) — opening the HTML
via `file://` won't work, since the pages `fetch()` their JSON at load time).

## Keeping it current

Rerun `update_data.py` on a schedule (cron, etc.) with your AAVSO API key:

```bash
AAVSO_API_KEY=your_key python3 update_data.py
```

Each run downloads only observations newer than what's already in the CSV,
merges them in, re-fits every ELL epoch, and regenerates `oc_data.json` /
`lightcurve_data.json` in place. `--no-fetch` re-fits from the existing CSV
without hitting the network; `--skip-tls-verify` is only for networks with a
TLS-intercepting proxy (not needed on most servers).

## Known landmine

There's a stray file at `tcrb/analysis/TCrB_V_1946-2026.csv` (~1.6 KB) — same
filename as the real master CSV (`tcrb/TCrB_V_1946-2026.csv`, ~21.8 MB) but
one directory down, and it contains old O–C results table data, not
photometry. Easy to grab the wrong one by accident. Safe to delete.

## Superseded files (safe to delete)

Not read by anything in the "required" list above:

- `oc_analysis.py`, `oc_results.json` — replaced by `update_data.py` / `oc_data.json`
- `tcrb_2005_present_nightly.json`, `tcrb_2005_present_minima.json`,
  `tcrb_2005_present_template.html` — intermediate build artifacts from before
  the fetch-based rewrite
- `tcrb_2025_2026_lightcurve.html` + its `_nightly.json`/`_template.html`,
  `tcrb_2025_lightcurve.html` + its `_nightly.json` — earlier, narrower-range
  light curve versions, all superseded by `tcrb_2005_present_lightcurve.html`
- `tcrb_ell_oc.html` — the older year-indexed O–C chart (embedded-data style,
  never converted to fetch-based); keep only if you specifically want that
  variant alongside the coverage one
