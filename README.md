# Sarcoma myogenic axis analysis

R code and preserved inputs for **Leiomyosarcoma-Enriched Myogenic Axis in Public Sarcoma Cohorts** by Sanghyun Roh and Jeong-An Gim. This is the repository content prepared for **v1.1.0**.

## Start here

Download and extract the complete repository archive from the v1.1.0 GitHub release or its corresponding Zenodo record. Keep `S1_Data.xlsx`, `S1_Code/` and `Data/` together. In RStudio, set the working directory to `S1_Code`, open the required script and run its numbered modules in order. Use a fresh R session for each script. No Python or launcher script is required for the scientific analyses or figures.

The three execution routes are independent; running all five scripts in numerical order is unnecessary.

| Purpose | Scripts | Inputs and outputs |
| --- | --- | --- |
| Reproduce submitted figures and reported calculations | `01_Figures.R` and `02_Reported_statistics.R` | Reads the supplied `S1_Data.xlsx`; writes PDF/TIFF figures, clinical models, paired comparisons and GSEA results under `S1_Code/Results/`. |
| Reconstruct the historical primary expression analysis | `03_Primary_expression.R` | Restores nine preserved GEO matrices, downloads fixed TCGA files, and performs representative-record, clustering, signature and projection analyses. |
| Repeat the separate online reconstruction | `04_Online_inputs.R`, then `05_Online_analysis.R` | Downloads GEO CEL files and metadata, prepares RMA matrices, and runs the online analysis. Both scripts must use the same `SARCOMA_ONLINE_DIR`. |

## Input distribution

`Data/Primary/` contains ten binary parts of the exact historical input ZIP. Each part is at most 20,000,000 bytes, so it can be committed as an ordinary Git file and is below GitHub's 25 MiB browser-upload limit. All ten parts must remain in the repository at the release tag. They are included in the source archive that Zenodo preserves, without depending on Git LFS settings or separately attached release assets.

Module 01a of `03_Primary_expression.R` joins these local parts in R, checks the complete archive and each matrix, and extracts the nine required CSV files. Do not open or decompress individual parts. The reconstructed ZIP is 198,666,176 bytes; the expression values and IDs are unchanged.

If only the small journal `S1_Code.zip` is available, Module 01a downloads missing parts from the fixed **v1.1.0 tag**, using URLs in `Inputs/Primary_archive_parts.csv`. Those URLs become available after that tag is published with the supplied repository files. The complete repository download can restore the historical matrices locally before publication. TCGA downloads still require an internet connection on first use.

An already available `Sarcoma_Primary_Inputs_20261116.zip` may be placed beside `S1_Code`, or its path may be set with `SARCOMA_PRIMARY_INPUT_ARCHIVE`. `SARCOMA_PRIMARY_INPUT_URL` accepts an HTTPS mirror of the exact full ZIP. `SARCOMA_PRIMARY_PARTS_BASE_URL` accepts an HTTPS directory containing the same ten parts. Checksums remain mandatory for all sources.

## Environment

The scripts install missing CRAN and Bioconductor packages in the active R library. The package lists appear at the start of each script. Use a writable library and an R/Bioconductor combination compatible with those packages; Windows installations may require matching Rtools when building packages from source.

The reported online enrichment was recalculated with **R 4.5.0, fgsea 1.38.0 and BiocParallel 1.40.0**, using the included **MSigDB 2026.1.Hs** collections and collection-specific seeds 124 and 125. An installation of the newest packages is not a locked reproduction environment; stochastic results may differ across versions. See the manuscript for the environments used in the earlier analyses.

For figures and reported calculations, use the supplied workbook. Primary expression reconstruction requires roughly 426 MB for the extracted GEO matrices, the 199 MB restored archive, additional TCGA inputs and space for results. The raw CEL route downloads substantially more data and needs additional memory, disk space and processing time.

## Interpretation

The historical primary analysis and the online RMA reconstruction are separate analyses. Freshly normalized CEL files must not replace the preserved primary matrices. Historical normalization and annotation versions are incompletely documented. The 629 connected-record groups are not confirmed unique patients. The myogenic axis describes histology-associated expression variation; it is not established as an independent prognostic classifier or stable binary sarcoma taxonomy.

## Citation and release

Software authorship follows the original repository citation, with Sanghyun Roh as the code author. `CITATION.cff` supplies GitHub citation metadata and `.zenodo.json` supplies matching Zenodo metadata. The MIT license applies to the original code and documentation; see `DATA_SOURCES.md` for third-party materials.

The existing DOI [10.5281/zenodo.19876316](https://doi.org/10.5281/zenodo.19876316) identifies **v1.0.0**, which does not contain this revision or the historical matrices. The concept DOI is [10.5281/zenodo.19876315](https://doi.org/10.5281/zenodo.19876315). A new version DOI is assigned after v1.1.0 is released through the enabled GitHub–Zenodo integration. Cite that new version DOI for this package; no new DOI or release date is asserted in the prepared metadata.

## Funding

This research was supported by the development of a real-world big data-based clinic-genomic information management system for sarcoma through the National Research Foundation of Korea (NRF), funded by the Ministry of Education (RS-2024-00349335), awarded to J.-A.G. This research was supported by the Global Learning & Academic Research Institution for Master’s–PhD Students and Postdocs (G-LAMP) Program of the National Research Foundation of Korea (NRF), funded by the Ministry of Education (No. RS-2025-25441283). This research was supported by the Soonchunhyang University Research Fund. The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.
