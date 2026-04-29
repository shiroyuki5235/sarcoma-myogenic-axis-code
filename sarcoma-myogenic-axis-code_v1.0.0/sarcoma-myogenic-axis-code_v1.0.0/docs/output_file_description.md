# Output file description

The scripts write outputs under the project-level `Output/` directory.

## Main output directories

```text
Output/
  tables/
  plots_png/
  plots_pdf/
  cache/
  logs/
  tables_addon/
  plots_png_addon/
  plots_pdf_addon/
  tables_validation/
  plots_png_validation/
  plots_pdf_validation/
```

## Major expected outputs

| Output path or pattern | Description |
|---|---|
| `Output/00_prepare_inputs_summary.csv` | Summary of prepared GEO input files. |
| `Output/tables/cohort_accounting_by_GSE.csv` | Per-GEO accounting after loading and preprocessing. |
| `Output/tables/sample_master_metadata_with_curation.csv` | Sample-level curation metadata and keep/drop decisions. |
| `Output/tables/sample_master_metadata_curated_with_cluster.csv` | Curated samples with assigned discovery clusters. |
| `Output/exprs_all_curated_prefiltered.rds` | Merged curated expression matrix before batch correction. |
| `Output/exprs_combat_curated.rds` | Batch-corrected curated discovery expression matrix. |
| `Output/tables/cluster_k_support_metrics.csv` | Silhouette and bootstrap ARI-based support for candidate k values. |
| `Output/tables/cluster_final_k_selected.csv` | Final selected k for discovery clustering. |
| `Output/tables/DEG_*.csv` | Differential expression analysis results for cluster contrasts. |
| `Output/tables/GSEA_*.csv` | Gene set enrichment outputs. |
| `Output/tables_addon/*.csv` | Additional robustness, composition, signature score, and sensitivity analysis outputs. |
| `Output/tables_validation/*.csv` | TCGA-SARC external validation outputs. |
| `Output/plots_*/*.png` and `Output/plots_*/*.pdf` | Main, additional, and validation figures. |

Large output files should generally not be deposited in the GitHub repository unless required and redistribution is allowed. For manuscript transparency, selected final tables/figures can be archived separately if needed.
