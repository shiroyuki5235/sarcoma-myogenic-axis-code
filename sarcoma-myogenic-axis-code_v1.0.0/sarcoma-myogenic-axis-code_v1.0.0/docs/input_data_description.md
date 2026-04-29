# Input data description

## General folder structure

The scripts assume that the working directory is the project root and use the following directory structure:

```text
Input/
  GEO_prepared/
  GEO_raw/
  TCGA/
Output/
```

## GEO discovery cohorts

The discovery workflow uses public GEO sarcoma transcriptome cohorts. The scripts refer to the following GEO Series:

```text
GSE71118
GSE21050
GSE16382
GSE34620
GSE142162
GSE66533
GSE20196
GSE34800
GSE94321
```

For each GEO Series, the preferred prepared input files are:

```text
Input/GEO_prepared/<GSE_ID>/<GSE_ID>_ExpressionMatrix_Annotated.csv
Input/GEO_prepared/<GSE_ID>/<GSE_ID>_SampleMetadata.csv
```

If these files are absent, the discovery script attempts to download supplementary GEO files and process available CEL files where possible.

## TCGA-SARC external validation input

The TCGA validation script uses the cBioPortal TCGA-SARC public study archive:

```text
Input/TCGA/sarc_tcga_pub.tar.gz
Input/TCGA/sarc_tcga_pub/
```

The script expects the following files after extraction:

```text
data_mrna_seq_v2_rsem.txt
data_clinical_sample.txt
data_clinical_patient.txt
```

## Data redistribution note

Raw GEO/TCGA-derived files and large generated outputs should not be committed to GitHub unless redistribution is explicitly permitted. Use `.gitignore` to exclude `Input/` and `Output/` by default.
