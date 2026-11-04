# Input data

This repository bundles **one dataset: Supplementary Data 1.** Everything else the
notebooks read is published alongside the paper and is downloaded rather than
redistributed here. **No individual-level clinical record is required to reproduce any
figure or statistic in this repository.**

| File | Bundled | Description |
|---|---|---|
| `Supplementary_Data_1.xlsx` | yes | Metadata for all 170 serotype 2 genomes — sample ID, country, year, source, GPSC, ST, resistance profile, CbpA/PspC paralog status |
| `Source_Data.xlsx` | no — download | Aggregate values underlying every display item, published with the paper |
| `naturalearth_lowres/` | no — download | Natural Earth 1:110m Admin 0 – Countries shapefile |

## Getting the two downloaded files

**Source Data.** Published with the article on nature.com as the Source Data file. Save it
here as `data/Source_Data.xlsx`. It contains one sheet per display item:

| Sheet | Drives |
|---|---|
| `Figure 1` | annual serotype 2 counts by detection method, and as a percentage of meningeal and non-meningeal IPD |
| `Figure 2a`, `Figure 2b` | per-serotype counts, IPD versus otitis media and versus carriage |
| `Figure 5` | genome counts per country by sequence type |
| `Supplementary Figure 2`, `Supplementary Figure 3` | enlarged regions of Figure 2 |
| `Supplementary Figure 4` | age in months and serotype group, per case |
| `Table 1` | full contingency counts with exact two-sided *P* values |

**Natural Earth basemap.** Download the 1:110m "Admin 0 – Countries" shapefile from
https://www.naturalearthdata.com (public domain) and unzip it into
`data/naturalearth_lowres/`.

## Why no patient-level data

Figures 1, 2a, 2b and Supplementary Figures 2–4 are built from counts and from a single
age-and-group column, and every statistic in Table 1 can be recomputed from the published
contingency counts — `Table_01_statistics.ipynb` does exactly that and reproduces each
exact *P* value. The surveillance line list behind those counts contains individual-level
clinical records for children and is therefore not distributed here or in the repository
history.

## Known issues in the source counts

Both corrected upstream, so the values in Source Data are already clean. Recorded here
because they affected earlier versions of this analysis code:

1. The carriage counts used `NT` for non-typeable isolates while the code filtered `NTD`,
   so non-typeables were not excluded from the serotype comparison. `NT`, `NTD` and `QNS`
   are now all excluded.
2. The carriage counts contained `14.0` and `3.0` as rows separate from `14` and `3`;
   these are merged.
