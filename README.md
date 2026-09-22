# CONUS Forest Edge Dynamics and Drivers — manuscript code

Analysis code for:

> **Human-directed disturbances accelerate forest edge expansion in the U.S.**
> Hangkai You, Min Chen\*, Zhe Zhu, Shi Qiu, Ian G. Brosnan, Ramakrishna Nemani, Volker Radeloff, Ming Liu, Taejin Park\*

This repository contains **only the code used to produce the results, figures and supplementary material of the manuscript.** Exploratory and unrelated notebooks are not included (see [What was excluded](#what-was-excluded)). Notebook outputs are stripped; all code is retained.

This repository was organized by Claude Code and checked manually by the authors.

---

## Data sources

| Dataset | Description | Access |
|---|---|---|
| LCMAP Collection 1.3 | Annual 30 m land cover, CONUS, 1985–2021 | USGS (public) |
| LCMAP C1.3 reference | 25,000-plot simple random sample, annually interpreted | USGS (public) |
| GERS disturbance agent maps | Annual 30 m dominant disturbance agent, 1985–2022 | Global Environmental Remote Sensing Lab, UConn |
| Census regions | Four-region division (Northeast, Midwest, South, West) | U.S. Census Bureau |

Input rasters and intermediate products are **not** included here (they are hundreds of GB). Scripts reference local paths that must be edited to point at your own copies.

---

## Pipeline

Run in numerical order. Steps 01–03 produce the rasters and tabular products consumed by 04–09.

### `00_core_scripts/` — production scripts
Batch versions of the core per-pixel operations, applied across all years.

| File | Purpose |
|---|---|
| `LCMAP_Edge_Dynamic_Mapping.py` | Driver for annual edge mapping across the LCMAP archive |
| `Edge_Encoding.py` | Four-bit directional edge encoding (E/W/S/N) |
| `Forest_Pixel_Classification.py` | Forest depth classification into 0–30, 30–60, 60–90, 90–120, >120 m bins |
| `Adjunct_LC_Encoding.py` | Adjacent (non-forest) land-cover class for each edge pixel |
| `zonal_hist.py` | Zonal histogram helper for spatial-unit summaries |
| `Forest Edge Time Trend.py` | Edge-length time series assembly |

### `01_preprocessing/`
CONUS tile extraction, annual binary forest masks, raster organisation and compression.

### `02_edge_depth_age/`
Edge identification (four-neighbour adjacency), forest depth, and **year since edge creation**.
Edge age is stored per direction as `EdgeAge_{north,south,east,west}_YYYY.tif` (uint8; 0 = no edge; ages 1–34, where 34 denotes edges predating 1988). Creation year = `2022 − age`.
**Edge-age percentages are tallied per directional edge segment** (pooling all four direction rasters), consistent with the directional edge-length encoding.

### `03_aggregation_1km/`
30 m → 1 km aggregation of area, edge length and interior-forest proportion, and the reverse resampling used to pair the 1 km dominant-disturbance product back to 30 m.

### `04_disturbance_attribution/`
Assembly of disturbance history, the 1 km dominant-disturbance product (33 × 33 pixel = 990 m majority), and the merge of edge/area dynamics with disturbance agents. Produces the master attribution table used by steps 05–09.

### `05_efcr/`
**Exterior Forest Conversion Ratio (EFCR)** = interior forest converted to exterior forest per unit forest loss.

Computed from the transition classes as
`EFCR = (pixels 5 → {1,2,3,4}) / (pixels {1..5} → 0)`
where class 5 is interior forest (depth > 120 m), classes 1–4 are exterior forest (≤ 120 m) and class 0 is non-forest.

**Aggregation rule, applied identically at every spatial scale:**
1. within the spatial unit, sum pixels → `frag_y`, `loss_y` for each year;
2. annual ratio `EFCR_y = frag_y / loss_y`;
3. reported mean = temporal mean of `EFCR_y` over 1988–2020 (n = 33);
4. Trend estimated by linear regression on year with AR(p) errors (p = 0–3 selected by AIC; Ljung–Box residual check, all p > 0.05), reported as the coefficient with its 95% CI. Mann–Kendall and Theil–Sen were replaced because they overstate significance in temporally autocorrelated series (Ives et al. 2021); see `08_uncertainty_and_revision/mannkendall_vs_AR1_simulation.py` for the calibrated false-positive simulation.

Only step 1 changes between CONUS and a region. This makes EFCR **area-weighted by forest loss** (the denominator of the ratio), so regions contribute in proportion to their forest loss. Do not average regional EFCR values unweighted: that over-weights regions with negligible disturbance area.

### `06_area_edge_analysis/`
Forest area and edge-length time series, regional division, and area–edge divergence analysis.

### `07_adjunct_landcover/`
Adjacent land-cover composition of pre-1988 ("legacy") edges. Supports Fig. S5.

### `08_uncertainty_and_revision/`
Design-based (bias-adjusted) area estimation and the revision analyses.

| File | Purpose |
|---|---|
| `bias_adjusted_loss_by_agent.py` | Design-based forest-loss-by-agent estimates from the reference `change_process` field |
| `bias_adjusted_net_by_agent.py` | Design-based gross gain / net flux by agent |
| `map_loss_gain_net_by_agent.py` | Map (pixel-count) gross loss, gross gain and net flux by agent |
| `regional_loss_gain_net.py` | The same, split by Census region |
| `area_edge_trends.py` | Theil–Sen slopes and Mann–Kendall p for area and edge, by region. **Superseded** by the AR analyses below; retained only for the MK-vs-AR comparison |
| `efcr_regional_trends.py` | Regional and CONUS EFCR annual series (Mann–Kendall version; superseded, retained for comparison) |
| `efcr_annual_by_disturbance_extract.py` | Extracts the national annual EFCR series for each disturbance agent from the attribution table |
| `efcr_AR_trends_all_series.py` | AR(p)-error regression trends for all EFCR series (regional, national, by disturbance; 1988–2021 and 2001–2021), with residual diagnostics |
| `edge_acceleration_AR_increments.py` | Acceleration γ from the annual-increment model ΔE_t = α + γt + ε_t with AR(p) errors; Newey–West and quadratic-fit cross-checks |
| `stock_series_differencing_check.py` | ADF unit-root tests and first-difference diagnostics showing why area and edge *levels* are not trend-tested |
| `mannkendall_vs_AR1_simulation.py` | MK vs OLS vs AR(1) on every series, plus the zero-trend false-positive simulation calibrated to the observed autocorrelation |

**Scope of the bias adjustment.** Main-text areas are map pixel counts, for consistency with edge length and EFCR (geometric properties of the classified map with no reference-based analogue). Design-based estimates are reported alongside as an uncertainty check. Two limits are load-bearing:
- Forest **loss** by agent is estimable, but only **Logging** (n = 1,107, ±6%) and **Fire** (n = 149, ±17%) are well constrained; minor agents have wide intervals (Stress n = 18, Water Dynamic n = 7, Natural Hazard n = 0).
- Forest **gain** cannot be attributed by agent at all: the reference labels recovery as `Growth/Recovery` (978 of 1,248 gain plots) and only ~12 carry a disturbance label. Per-agent gain and net flux are therefore map-based only.

### `09_figures/`
Main figures 1–3 and supplementary figures S1–S5, plus the revision figures.

| File | Output |
|---|---|
| `Figure 1.ipynb` … `Figure S5.ipynb` | Manuscript figures |
| `figure3_pooled_efcr_compute.py` | Panel (a) data: area-weighted (pooled) mean EFCR by disturbance type |
| `figure3b_extract_region_x_disturbance_efcr.py` | Panel (b) input: annual EFCR series for each region × disturbance combination (32 series) |
| `figure3b_AR_trend_fits.py` | Panel (b) statistics: AR(p)-error trend, 95% CI and p for each of the 32 series, plus national and regional aggregates |
| `figure3_AR_plot.py` | **Figure 3 as published.** Panel (a) unchanged; panel (b) plots mean EFCR against the AR-regression trend, filled markers p ≤ 0.05 |
| `figure_efcr_region_vs_conus.py` | Regional vs CONUS EFCR comparison |
| `figure_loss_gain_net_bars.py` / `_annual.py` | Loss / gain / net flux by agent |
| `figure_biasadj_total_flux.py` | Total flux, map vs bias-adjusted |

---

## Reproducing

```bash
conda env create -f environment.yml
conda activate conus-edge
```

Then edit the input paths at the top of each script or notebook and run `00`–`09` in order.

`numpy` is pinned to 1.25.2 deliberately; 2.x breaks matplotlib/pandas/scipy/cartopy in this stack.

---

## What was excluded

Removed from this branch because it does not contribute to the manuscript:

- `NLM simulation.ipynb` — neutral landscape model; no such analysis appears in the manuscript
- `Playground.ipynb`, `Untitled.ipynb` — scratch
- `Edge analysis stats single year.ipynb` — single-year prototype, superseded by the 1988–2021 series
- `Analyze the driver of the forest edge dynamics.ipynb` — machine-learning driver models, not used
- `Log edge length vs log edge area.ipynb` — exploratory
- FINESST 2025 and NSF CAREER figure notebooks — different proposals
- County-scale visualisations of area, interior rate and edge density — the manuscript reports Census regions
- 225 PNG/GIF renders, slide decks, and toy test rasters (`mock_forest.tif`, `classified_forest.tif`)

A note on repository history: the manuscript pipeline previously lived in a separate private repository and was referenced from `main` as a `Summer Intern` submodule with no `.gitmodules` entry, so it cloned as an empty directory. That broken reference is removed here and the code is committed directly.

---

## Licence

Released under the [MIT License](LICENSE). Copyright (c) 2026 Hangkai You.

## Citation

If you use this code, please cite the manuscript:

> You, H., Chen, M., Zhu, Z., Qiu, S., Brosnan, I. G., Nemani, R., Radeloff, V., Liu, M., & Park, T. Human-directed disturbances accelerate forest edge expansion in the U.S.
