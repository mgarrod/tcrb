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
| `oc_data.json` + `oc_data.js` | `tcrb/analysis/` | Generated output — loaded by the O–C page |
| `lightcurve_data.json` + `lightcurve_data.js` | `tcrb/analysis/` | Generated output — loaded by the light curve page |
| `tcrb_ell_oc_coverage.html` | `tcrb/analysis/` | O–C diagram page |
| `tcrb_2005_present_lightcurve.html` | `tcrb/analysis/` | Observed/calculated light curve page |

Each data file comes in two forms: the plain `.json` (for any other tool that
wants to read it) and a `.js` twin (`const OC_DATA = {...};`) that the HTML
pages actually load via `<script src>`. Both are regenerated together by
`update_data.py` — never edit either by hand.

## Deploying

Copy the two `.html` files plus `oc_data.json`/`.js` and
`lightcurve_data.json`/`.js` to the same directory — that's it. The pages load
their data via `<script src>`, not `fetch()`, so they work both **served over
http(s)** and **opened directly via `file://`** (double-click the HTML file).

## Keeping it current

Rerun `update_data.py` on a schedule (cron, etc.) with your AAVSO API key:

```bash
AAVSO_API_KEY=your_key python3 update_data.py
```

Each run downloads only observations newer than what's already in the CSV,
merges them in, re-fits every ELL epoch, and regenerates `oc_data.json`/`.js`
and `lightcurve_data.json`/`.js` in place. `--no-fetch` re-fits from the
existing CSV without hitting the network; `--skip-tls-verify` is only for
networks with a TLS-intercepting proxy (not needed on most servers).

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

## Second pipeline: MGAB-V249 (eclipse timing)

A separate, independent mini-pipeline for **MGAB-V249**, an ultra-compact
eclipsing hot-subdwarf+white-dwarf binary (VSX type E/WD+IBWD, P = 0.0273195 d
≈ 39.3 min — nothing like T CrB's 227.5 d ellipsoidal period). It shares no
code or data files with the T CrB pipeline above.

| File | Location | Role |
|---|---|---|
| `MGAB-V249_V.csv` | `tcrb/` | Master photometry for this target |
| `update_data_mgab_v249.py` | `tcrb/analysis/` | Pipeline script (phase-wrapped minimum timing, not nightly binning) |
| `oc_data_mgab_v249.json` + `.js` | `tcrb/analysis/` | Generated output — loaded by `mgab_v249_oc.html` |
| `mgab_v249_oc.html` | `tcrb/analysis/` | O–C diagram + phase-folded light curve page |

At a 39-minute period, a single AAVSO night spans hundreds of cycles, so
"one point per night" (T CrB's method) doesn't apply. Instead, each
observing run's points are folded onto the nearest predicted eclipse
(wrapped by period) and quadratic-fit as one unit, giving one O–C timing
per *run*, not per cycle. See the script's docstring for details.

Current data (AAVSO, single observer) spans only ~3 weeks across 4 runs —
the O–C diagram is real but thin; it sharpens as more V photometry comes
in. Rerun the same way as the main pipeline:

```bash
AAVSO_API_KEY=your_key python3 update_data_mgab_v249.py
```
