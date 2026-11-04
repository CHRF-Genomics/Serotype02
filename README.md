# Reemergence and global distribution of an invasive lineage of *Streptococcus pneumoniae* serotype 2

Analysis code accompanying:

> Hooda Y, Tanmoy AM, Pushpita KB, Kanon N, Rahman H, Naziat H, Hung HCH, Malaker R,
> Hasanuzzaman M, Malaker AR, Keya DP, Nath SD, Hossain B, Saha S, Uddin MJ, Klugman KP,
> Santosham M, McGee L, Bentley SD, Lo SW, Saha S, Saha SK.
> *Reemergence and global distribution of an invasive lineage of Streptococcus pneumoniae
> serotype 2.* Nature Communications (2026). Manuscript NCOMMS-26-018270A.

Maintained by the [Child Health Research Foundation](https://chrfbd.org), Dhaka, Bangladesh.

## Contents

```
data/        Supplementary Data 1, plus two downloaded files (see data/README.md)
notebooks/   Jupyter notebooks, one per display item
scripts/     shell scripts for RAxML and BEAST, and iTOL annotation files
results/     generated figures
```

| Notebook | Display item |
|---|---|
| `Figure_01_temporal_trends.ipynb` | Figure 1 — annual serotype 2 IPD counts, 1993–2022 |
| `Figure_02a_IPD_vs_otitis_media.ipynb` | Figure 2a — serotype distribution, IPD vs otitis media |
| `Figure_02b_IPD_vs_carriage.ipynb` | Figure 2b — serotype distribution, IPD vs nasopharyngeal carriage |
| `Figure_05_global_distribution.ipynb` | Figure 5 — global distribution of genomes by sequence type |
| `Supplementary_Figure_04_age_distribution.ipynb` | Supplementary Figure 4 — age distribution by serotype |
| `Table_01_statistics.ipynb` | Table 1 — recomputes every exact *P* value from the published counts |

Figures 3 and 4 and Supplementary Figure 5 are phylogenies rendered in
[iTOL v6](https://itol.embl.de); the trees and their annotation files are in `scripts/`.

## Running the analyses

```bash
git clone https://github.com/CHRF-Genomics/<repository>.git
cd <repository>
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter lab
```

Run notebooks from the repository root — all paths are relative to it.

This repository bundles one dataset, Supplementary Data 1. Two further files are
downloaded rather than redistributed: the **Source Data file** published with the paper
(save as `data/Source_Data.xlsx`) and the public-domain **Natural Earth** basemap used by
Figure 5. See `data/README.md`.

**No individual-level clinical record is needed.** Every figure runs from aggregate
counts, and `Table_01_statistics.ipynb` reproduces each exact *P* value in Table 1 from
the published contingency counts alone. The surveillance line list behind those counts
contains clinical records for children and is not distributed here.

## Sequencing data

| Dataset | Accession |
|---|---|
| This study (68 Bangladeshi serotype 2 isolates, Illumina + ONT) | [PRJEB81323](https://www.ebi.ac.uk/ena/browser/view/PRJEB81323) |
| Global Pneumococcal Sequencing project | [PRJEB3084](https://www.ebi.ac.uk/ena/browser/view/PRJEB3084) |
| *S. pneumoniae* D39 reference chromosome | [CP000410.2](https://www.ncbi.nlm.nih.gov/nuccore/CP000410.2) |
| D39 RefSeq record | [NC_008533.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_008533.2) |
| GPSC96/ST74 isolate used for paralog reconstruction | [ERR750825](https://www.ebi.ac.uk/ena/browser/view/ERR750825) |

Genome assemblies and metadata for all 170 serotype 2 isolates are published as
Supplementary Data 1. Underlying values for every display item are published as
the Source Data file.

## Genomic pipeline

Read QC, assembly, serotyping, GPSC/ST assignment and AMR prediction used the
[GPS Pipeline v rc11](https://github.com/GlobalPneumoSeq/gps-pipeline), which pins the
versions of Bowtie2, SAMtools, BCFtools, Gubbins and RAxML used here.

## Notes on reproducibility

Three issues were identified while preparing this repository and are corrected in the
code published here:

1. **Non-typeable isolates.** The serotype comparison previously filtered the label
   `NTD`, but the count files use `NT`, so non-typeable isolates were not excluded.
   They are now removed along with `QNS` and missing values.
2. **Float-coerced serotype labels.** Serotype labels such as `14.0` and `3.0` appeared
   as rows separate from `14` and `3`. These are now merged before plotting.
3. **Deprecated APIs.** `gpd.datasets.get_path()` was removed in GeoPandas 1.0 and
   `sns.kdeplot(shade=...)` in seaborn 0.13. The notebooks now read a local Natural Earth
   shapefile and use `fill=`.

## Basemap attribution

Figure 5 uses the Natural Earth 1:110m "Admin 0 – Countries" dataset, which is in the
public domain (https://www.naturalearthdata.com). It is not redistributed here; see
`data/README.md`.

## Licence

Code is released under the MIT Licence (`LICENSE`).
