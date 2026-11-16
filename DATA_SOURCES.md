# Data sources and reuse

    ## Historical GEO matrices

    The archived matrices are the prepared CSV inputs used for the historical primary analysis. They preserve the original numerical precision, identifiers and ordering. The exact historical normalization and annotation versions could not be fully recovered, so these files are not described as raw GEO data or as newly reproduced RMA matrices.

    | Series | Original public record |
    | --- | --- |
    | GSE71118 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71118 |
| GSE21050 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21050 |
| GSE16382 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16382 |
| GSE34620 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34620 |
| GSE142162 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142162 |
| GSE66533 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66533 |
| GSE20196 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20196 |
| GSE34800 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34800 |
| GSE94321 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94321 |

    `S1_Code/Inputs/Required_input_files.csv` identifies the nine expected extracted matrices and three TCGA inputs. `S1_Code/Inputs/Primary_archive_parts.csv` controls reconstruction of the historical ZIP. These files are required by the download code.

    ## TCGA SARC

    Three files from `sarc_tcga_pub` are downloaded from [cBioPortal datahub commit 065ec89bc152961dd3b15eacafd77edcc05950b5](https://github.com/cBioPortal/datahub/tree/065ec89bc152961dd3b15eacafd77edcc05950b5/public/sarc_tcga_pub): `data_mrna_seq_v2_rsem.txt`, `data_clinical_sample.txt` and `data_clinical_patient.txt`. They are not embedded here. A newer datahub snapshot may change gene identifiers and must not silently substitute for this version.

    ## Gene sets and reported data

    The included Hallmark and Reactome collections preserve the MSigDB 2026.1.Hs gene-set membership used for enrichment. Compact source tables retain `gs_name`, `gene_symbol` and `db_version`. The fixed online differential-expression ranks retain `GeneSymbol` and `t`. `S1_Data.xlsx` contains the reported values used by scripts 01 and 02; it does not replace expression inputs for script 03.

    Cite the original GEO studies, TCGA/cBioPortal, MSigDB and Reactome when reusing the relevant materials. GEO-derived matrices, cBioPortal data and gene-set definitions retain the terms and attribution requirements of their original sources. The repository MIT license does not relicense those third-party materials. Sources: https://www.ncbi.nlm.nih.gov/geo/ ; https://www.cbioportal.org/ ; https://www.gsea-msigdb.org/gsea/msigdb/ ; https://reactome.org/ .
