S1 Code. Sarcoma analysis and figure-generation scripts

Preparation
Extract this archive as S1_Code and keep the separately supplied S1_Data.xlsx beside that folder. Set the RStudio working directory to S1_Code. Run numbered modules in the required script, using a fresh R session for each script. All scientific analyses and figures use R; no launcher or source() chain is required. Missing packages are installed in the current R library.

Five scripts
01_Figures.R recreates Figs 1–7 and S1 Fig from S1 Data as PDF and 600-dpi TIFF. It plots the reported numerical results without refitting expression models.
02_Reported_statistics.R refits the reported Cox models and sensitivity analyses, performs paired score/AUC comparisons, and recalculates the online Hallmark/Reactome enrichment. Its small fixed ranks and MSigDB 2026.1.Hs gene sets are included.
03_Primary_expression.R performs the historical primary expression analysis, record linkage, representative selection, ComBat, signature construction, TCGA projection and reviewer-requested sensitivity analyses. It restores the exact historical GEO inputs from the companion archive and downloads the exact TCGA expression and clinical files from a fixed cBioPortal commit. SHA-256 checks prevent substitution with changed files.
04_Online_inputs.R downloads GEO metadata and CEL archives for nine series, prepares RMA matrices, and obtains the same fixed TCGA files. This route downloads its large raw inputs at runtime. Included gene sets retain the reported enrichment release.
05_Online_analysis.R runs after 04 with the same PROJECT_DIR and performs the separate online reconstruction, enrichment and TCGA projection. It does not replace the historical primary analysis.

Historical primary input archive
The full GitHub/Zenodo repository contains Data/Primary beside S1_Code. Module 01a of 03 joins its ten binary parts in R, verifies the exact 198,666,176-byte historical archive, and extracts the nine required matrices. Keep all ten parts; do not decompress them individually.

When only this small S1 Code archive is available, Module 01a downloads missing parts from the v1.1.0 GitHub tag using Inputs/Primary_archive_parts.csv. The default URLs become available when that release/tag and its data parts are published. Until publication, use the supplied full repository set or the exact local companion archive. A failed download stops with a message and does not substitute newly processed GEO data.

An existing Sarcoma_Primary_Inputs_20261116.zip can be placed beside S1_Code. An alternative path can be set using SARCOMA_PRIMARY_INPUT_ARCHIVE. SARCOMA_PRIMARY_INPUT_URL accepts an HTTPS mirror of the exact full ZIP; SARCOMA_PRIMARY_PARTS_BASE_URL accepts an HTTPS directory containing the ten unchanged parts. All restored and downloaded inputs are checked before analysis.

TCGA inputs
No TCGA expression or clinical data files are embedded in S1 Code. Scripts 03 and 04 download only the three required files, pinned to cBioPortal/datahub commit 065ec89bc152961dd3b15eacafd77edcc05950b5. Required_input_files.csv supplies the exact URLs, sizes and SHA-256 checksums. A completed download is reused only if it matches the recorded checksum. Internet access is required on the first run unless those exact files are already present.

Execution routes
Run 01 and 02 to reproduce the reported figures/statistical calculations from S1 Data. Run 03 with the full repository, a local companion archive or access to the published v1.1.0 input URLs for primary expression reconstruction. Run 04 then 05 for the separate online reconstruction. It is not necessary to execute all five scripts in sequence. SARCOMA_CODE_DIR, SARCOMA_PROJECT_DIR and SARCOMA_ONLINE_DIR can override the corresponding paths without editing the scripts.

Interpretation and versions
The primary and online analyses remain separate. Original historical normalization and annotation versions remain incompletely documented. The preserved matrices retain the original values and IDs; no rounding or single-precision conversion has been applied. Linked-record groups are not confirmed unique patients.
Reported online enrichment uses MSigDB 2026.1.Hs, fgsea 1.38.0, BiocParallel 1.40.0 and R 4.5.0, with collection-specific seeds 124 and 125. Other package versions can affect stochastic estimates. The earlier discovery/online analyses used R 4.6.1/4.6.2 as described in the manuscript.
Original public code archive: doi:10.5281/zenodo.19876316 (v1.0.0). It does not identify this revised release or the historical matrices. Cite the v1.1.0 version DOI once issued through GitHub–Zenodo. This package prepares the files but does not itself publish a release.
