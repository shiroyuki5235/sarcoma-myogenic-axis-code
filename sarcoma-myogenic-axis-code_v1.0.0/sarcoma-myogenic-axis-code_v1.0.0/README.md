# Analysis code for Leiomyosarcoma-Enriched Myogenic Axis in Public Sarcoma Cohorts

## Overview

This repository contains the R scripts used for the manuscript **“Leiomyosarcoma-Enriched Myogenic Axis in Public Sarcoma Cohorts.”**

The workflow supports reproducible computational analysis of public sarcoma transcriptome cohorts, including input preparation, discovery cohort curation, batch correction, unsupervised clustering, differential expression analysis, enrichment analysis, additional robustness analyses, and TCGA-SARC external validation.

## Repository structure

```text
R/
  01_prepare_inputs.R
  02_discovery_main.R
  03_discovery_additional_analyses.R
  04_tcga_external_validation.R
docs/
  script_manifest.csv
  input_data_description.md
  output_file_description.md
session_info/
  session_info_placeholder.txt
example_output/
  README.md
```

## Script order

Run the scripts from the repository root in the following order:

```r
source("R/01_prepare_inputs.R")
source("R/02_discovery_main.R")
source("R/03_discovery_additional_analyses.R")
source("R/04_tcga_external_validation.R")
```

Each script assumes that the current working directory is the repository/project root. The scripts create and use `Input/` and `Output/` directories under the project root.

## Input data

Raw expression and phenotype data are not included in this repository unless redistribution is explicitly permitted. The scripts are designed for public GEO sarcoma datasets and the TCGA-SARC cBioPortal study archive. See `docs/input_data_description.md` for details.

## Output

The workflow writes tables, figures, logs, and serialized R objects to the `Output/` directory. See `docs/output_file_description.md` for the expected output structure.

## Software environment

The scripts were prepared for R-based bioinformatics analysis. Before running the workflow, install the packages listed in the Methods section of the manuscript and in the scripts. A package check can be performed using `requireNamespace()` or `sessionInfo()` after installation.

## Release and citation

For manuscript submission, cite the archived release rather than only the live GitHub URL. After creating a GitHub release, archive the release in Zenodo and cite the resulting DOI.

Suggested citation before DOI assignment:

> Roh S. Analysis code for Leiomyosarcoma-Enriched Myogenic Axis in Public Sarcoma Cohorts. Version 1.0.0. GitHub. 2026.

After Zenodo archiving, replace the GitHub-only citation with the Zenodo DOI citation.

## License

This code is released under the MIT License unless otherwise stated.

## Contact

Sanghyun Roh  
shiroyuki5235@sch.ac.kr  
Department of Medical Science, Soonchunhyang University, Republic of Korea
