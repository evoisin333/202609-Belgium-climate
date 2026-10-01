# Belgium's Changing Seasons — seasonal climate change in Belgium (1961–2024)

How has Belgium's seasonal climate shifted since the 1960s, and how is vegetation responding today?
This repository holds the full processing chain behind the project: open datasets, Python and ArcGIS Pro scripts, quality control and statistical analysis.

**Story map:** [Belgium's Changing Seasons](https://storymaps.arcgis.com/stories/296c84e52d4f4972813faf4ed28a5b3c)
**Authors:** Emma Voisin and Mattia Grasso, 2026

---

## Key findings

Change between the WMO climate normals **1961–1990 and 1991–2020**, national means weighted by area (values from figure 2 below):

| Indicator | Change |
|---|---|
| Mean summer temperature | **+1.2 °C** |
| Hottest summer day (TXx) | **+2.8 °C** |
| Summer days reaching 30 °C | **+3.5 days** |
| Spring precipitation | **−28 mm** |
| Winter precipitation | **+30 mm** |
| Winter frost days | **−7 days** |
| Spring soil moisture (top 7 cm) | **−0.016 m³/m³** |

- **Summers are warming everywhere.** Mean summer temperature rises by 0.35 to 0.44 °C per decade, significant in every grid cell (Mann-Kendall, p < 0.05).
- **Extremes are rising faster than the mean.** TXx rises by 0.5 to 1.4 °C per decade, up to three times faster than mean summer temperature.
- **Rain is shifting from spring to winter**, while annual totals barely change, and soils are drying in spring and summer.
- **Summer vegetation is declining on farmland.** Over 2001–2024, the median summer NDVI trend is −0.012 per decade, with 29 % of pixels showing a significant trend (26 % declining, 3 % increasing). The decline concentrates on cropland and grassland (−0.019 and −0.021 per decade), while forests lose about five times less. Grassland and cropland NDVI is strongly linked to summer temperature (r ≈ −0.72 to −0.77) and soil moisture (r ≈ +0.66 to +0.72); forests show no significant link with summer temperature.

![Mean seasonal temperature in Belgium, 1961–2024](06_figures_en/fig1_seasonal_temperature.png)

![Change between the two climate normals](06_figures_en/fig2_normals_change.png)

---

## Data

| Dataset | Use | Resolution | Period | Access |
|---|---|---|---|---|
| ERA5-Land monthly means (t2m, tp, swvl1, swvl2) | Mean temperature, precipitation, soil moisture | 0.1° | 1961–2024 | [Copernicus CDS](https://cds.climate.copernicus.eu/) |
| E-OBS v31.0e — seasonal ETCCDI indices (FD, SU, TR, TXx, PRCPTOT) and daily Tmax | Thermal extremes, precipitation, days ≥ 30 °C | 0.1° | 1961–2024 | [ECA&D / Copernicus](https://surfobs.climate.copernicus.eu/) |
| MODIS MOD13Q1 v061 (NDVI, pixel reliability) | Vegetation response | 250 m nominal (≈ 230 m), 16-day | 2001–2024 | [NASA AppEEARS](https://appeears.earthdatacloud.nasa.gov/) |
| CORINE Land Cover 2018 | Land-cover classes (44 → 8) | 100 m | 2018 | [Copernicus Land](https://land.copernicus.eu/) |
| NUTS 2 regions (Brussels-Capital and the 10 provinces) | Zonal statistics | — | — | [Eurostat GISCO](https://ec.europa.eu/eurostat/web/gisco) |

Indicator definitions (ETCCDI): **FD** frost days (Tmin < 0 °C), **SU** summer days (Tmax > 25 °C), **TR** tropical nights (Tmin > 20 °C), **TXx** hottest day of the season, **PRCPTOT** precipitation on wet days (≥ 1 mm). **su30** (days with Tmax ≥ 30 °C) is computed in this project from E-OBS daily Tmax.

Raw and intermediate data are not stored here (several GB). They can be downloaded from the sources above.

---

## Method

**Two time axes, kept deliberately separate:**
- **Climate axis, 1961–2024** (ERA5-Land, E-OBS): WMO normals 1961–1990 vs 1991–2020, plus the recent decade 2015–2024 for illustration.
- **Impact axis, 2001–2024** (MODIS NDVI): the satellite record is too short for the normals, so it documents contemporary vegetation response only.

**Seasons:** DJF, MAM, JJA, SON. December is assigned to the following year's winter.

**Trends:** Sen's slope with a two-sided Mann-Kendall test, per pixel, expressed per decade (ArcGIS Pro *Generate Trend Raster* for climate indicators; an equivalent Python implementation for NDVI, validated against ArcGIS).

**Regional figures:** zonal statistics by NUTS 2 region, and by region × land-cover class for NDVI. National values are area-weighted means.

---

## Processing chain

Scripts are numbered in execution order. Python scripts run in the conda environment from `environment.yml`, from the project root folder. ArcGIS scripts are pasted into the ArcGIS Pro Python window; set the folder paths at the top of each one first.

| # | Script | Environment | What it does |
|---|---|---|---|
| 01 | `01_decoupe.py` | Python | Clips all European NetCDF files to the Belgian bounding box |
| 02 | `02_agregation.py` | Python | ERA5 unit conversions, seasonal aggregation, days ≥ 30 °C from daily Tmax |
| 03 | `03_modis_ndvi.py` | Python | MODIS quality masking, scaling, seasonal NDVI |
| 04 | `04_controles.py` | Python | Consistency checks on all preprocessed outputs |
| 05–06 | `05_diagnostic_pixel_max.py`, `06_diagnostic_txx.py` | Python | Diagnostics that revealed the E-OBS station artefact |
| 07 | `07_controle_spatial_qc.py` | Python | Spatial-coherence quality control of TXx, SU and su30 |
| 08 | `08_masque_artefact_lux.py` | Python | Masks the residual artefact halo (summers 2003, 2004, 2006, 2007) |
| 09–11 | `09_export_era5.py`, `10_export_ndvi_su30.py`, `11_export_eobs_indices.py` | Python | One CF-compliant NetCDF per season, readable by ArcGIS Pro |
| 12 | `12_normales.py` | Python | Climate normals and change maps (GeoTIFF) |
| — | *manual step* | ArcGIS Pro | NetCDF → CRF (*Copy Raster*, multidimensional), then *Generate Trend Raster* (Mann-Kendall) |
| 13 | `13_tendances_arcgis.py` | ArcGIS Pro | Slope per decade, p-values, significance masking |
| 14 | `14_stats_zonales.py` | ArcGIS Pro | Zonal statistics by region, ERA5 and E-OBS |
| 15 | `15_ndvi_corine.py` | ArcGIS Pro | NDVI by region × land-cover class |
| 16 | `16_analyse.py` | Python | Normals by region, national series, NDVI–climate correlations |
| 17 | `17_figures.py`, `17_figures_en.py` | Python | Figures, in French and English |
| 18 | `18_tendance_ndvi.py` | Python | Per-pixel NDVI trend (Sen + Mann-Kendall) |
| 19 | `19_couche_signif_ndvi.py` | ArcGIS Pro | Significant NDVI trend areas larger than 20 km² |

Code comments are in French.

### Folder layout expected by the scripts

```
Project Belgium/          ← project root (this repository)
├── *.nc                  raw ERA5-Land and E-OBS files, as downloaded
├── 01_raw/               MODIS tiles, CORINE
├── 02_working/           clipped, aggregated and quality-controlled NetCDF
├── 03_arcgis/            one NetCDF per season (and the CRF made from them)
├── 04_normales/          normals and change maps (GeoTIFF)
├── 04_stats/             zonal statistics (CSV)
├── 05_resultats/         national series, normals, correlations (CSV)
└── 06_figures/, 06_figures_en/
Belgium/                  ← ArcGIS Pro project, next to the root
├── tendances/            Generate Trend Raster outputs
├── 05_tendances/         trend GeoTIFFs
└── Belgium.gdb
```

### Reproducing

```bash
conda env create -f environment.yml
conda activate geo
python 01_decoupe.py
```

The ArcGIS steps (13, 14, 15, 19 and the manual step) require ArcGIS Pro with the Spatial Analyst and Image Analyst extensions. Steps 12 and 16–18 rely only on open-source Python.

---

## Known issues and limitations

- **E-OBS station artefact.** A faulty station in the south of the Grand Duchy of Luxembourg produces impossible summer maxima (TXx up to 47 °C) that E-OBS interpolation spreads to neighbouring cells, including in the Belgian province of Luxembourg. A spatial-coherence test (TXx more than 4 °C above the median of its neighbours within ±0.5°) removed 40 cells, mostly around that station plus a faulty winter 2022 signal in Flanders, while keeping the genuine 2019 and 2022 heatwaves. A second test, cross-checked against ERA5-Land, masked a residual halo of 67 pixels for the summers of 2003, 2004, 2006 and 2007. Heat indices in the south of the province of Luxembourg should still be read with caution.
- **Different grids.** ERA5-Land and E-OBS cell centres are offset by about 0.05° (≈ 5.5 km). Data were not resampled; cross-source comparisons rely on regional zonal statistics.
- **Precipitation sources differ.** ERA5-Land `tp` is total precipitation; E-OBS PRCPTOT counts wet days (≥ 1 mm) only.
- **Rare events and Sen's slope.** Where days ≥ 30 °C are rare (coast, Ardennes), most summers score zero, so Sen's slope collapses towards zero even when Mann-Kendall detects an increase. The trend map shows *where* the increase is strongest, but underestimates its magnitude there.
- **No pre-whitening and no multiple-testing correction.** The Mann-Kendall test does not correct for serial autocorrelation, and it is repeated over thousands of pixels. This mainly concerns precipitation, soil moisture and NDVI, where significance is closer to the threshold.
- **Slopes in the figures.** Figures 1 and 3 show least-squares slopes on national means; the maps use per-pixel Sen's slopes. The two are not expected to be identical.
- **Short NDVI record.** 24 years (2001–2024) allow trend and correlation analysis, but not a comparison of climate normals.

---

## Credits

- Muñoz-Sabater, J. et al. (2021). ERA5-Land: a state-of-the-art global reanalysis dataset for land applications. *Earth System Science Data*, 13, 4349–4383.
- Cornes, R. et al. (2018). An Ensemble Version of the E-OBS Temperature and Precipitation Datasets. *Journal of Geophysical Research: Atmospheres*, 123. We acknowledge the E-OBS dataset and the data providers in the ECA&D project.
- Didan, K. (2021). MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V061. NASA EOSDIS Land Processes DAAC.
- CORINE Land Cover 2018, Copernicus Land Monitoring Service, European Environment Agency.
- NUTS boundaries © EuroGeographics for the administrative boundaries, Eurostat GISCO.
