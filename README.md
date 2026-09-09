  ┌───────────────────────────────────┬────────────────┬────────────────────────────────────────────────────────────┐
  │               File                │    Location    │                            Role                            │
  ├───────────────────────────────────┼────────────────┼────────────────────────────────────────────────────────────┤
  │ TCrB_V_1946-2026.csv              │ tcrb/          │ Master photometry data — update_data.py reads/updates this │
  ├───────────────────────────────────┼────────────────┼────────────────────────────────────────────────────────────┤
  │ update_data.py                    │ tcrb/analysis/ │ The one pipeline script                                    │
  ├───────────────────────────────────┼────────────────┼────────────────────────────────────────────────────────────┤
  │ oc_data.json                      │ tcrb/analysis/ │ Generated output — fetched by the O–C page                 │
  ├───────────────────────────────────┼────────────────┼────────────────────────────────────────────────────────────┤
  │ lightcurve_data.json              │ tcrb/analysis/ │ Generated output — fetched by the light curve page         │
  ├───────────────────────────────────┼────────────────┼────────────────────────────────────────────────────────────┤
  │ tcrb_ell_oc_coverage.html         │ tcrb/analysis/ │ O–C diagram page                                           │
  ├───────────────────────────────────┼────────────────┼────────────────────────────────────────────────────────────┤
  │ tcrb_2005_present_lightcurve.html │ tcrb/analysis/ │ Observed/calculated light curve page                       │
  └───────────────────────────────────┴────────────────┴────────────────────────────────────────────────────────────┘

  To deploy: copy the two .html files + oc_data.json + lightcurve_data.json (same directory, served over http). To keep them current: rerun
  update_data.py on a schedule, which regenerates the two JSON files in place.

  Found one landmine while checking

  There's a stray file at tcrb/analysis/TCrB_V_1946-2026.csv (1.6 KB) — same filename as your real master CSV (tcrb/TCrB_V_1946-2026.csv, 21.8 MB) but
  one directory down and containing old O–C results table data, not photometry. Easy to grab the wrong one by accident. I'd delete it.

  Safe to delete (superseded, not read by anything above)

  - oc_analysis.py, oc_results.json — replaced by update_data.py / oc_data.json
  - tcrb_2005_present_nightly.json, tcrb_2005_present_minima.json, tcrb_2005_present_template.html — intermediate build artifacts from before the
  fetch-based rewrite
  - tcrb_2025_2026_lightcurve.html + its _nightly.json/_template.html, tcrb_2025_lightcurve.html + its _nightly.json — earlier, narrower-range
  versions, all superseded by tcrb_2005_present_lightcurve.html
  - tcrb_ell_oc.html — the older year-indexed O–C chart (embedded-data style, never converted to fetch-based); keep only if you specifically still
  want that variant alongside the coverage one
