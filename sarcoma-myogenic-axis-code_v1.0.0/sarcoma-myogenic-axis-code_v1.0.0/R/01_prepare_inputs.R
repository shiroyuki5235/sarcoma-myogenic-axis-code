rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(123)

PROJECT_DIR <- getwd()
INPUT_DIR <- file.path(PROJECT_DIR, "Input")
INPUT_GEO_PREP_DIR <- file.path(INPUT_DIR, "GEO_prepared")
INPUT_GEO_RAW_DIR <- file.path(INPUT_DIR, "GEO_raw")
INPUT_TCGA_DIR <- file.path(INPUT_DIR, "TCGA")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Output")
LOG_DIR <- file.path(OUTPUT_DIR, "logs")

dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_GEO_PREP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_GEO_RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_TCGA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "prepare_inputs.log")

gse_list <- c(
  "GSE71118","GSE21050","GSE16382","GSE34620","GSE142162",
  "GSE66533","GSE20196","GSE34800","GSE94321"
)

copy_status <- vector("list", length(gse_list))
names(copy_status) <- gse_list

for (i in seq_along(gse_list)) {
  gse_id <- gse_list[i]
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste("Preparing", gse_id))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)

  legacy_dir <- file.path(PROJECT_DIR, gse_id)
  prep_dir <- file.path(INPUT_GEO_PREP_DIR, gse_id)
  dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)

  legacy_expr <- file.path(legacy_dir, paste0(gse_id, "_ExpressionMatrix_Annotated.csv"))
  legacy_sample <- file.path(legacy_dir, paste0(gse_id, "_SampleMetadata.csv"))

  prep_expr <- file.path(prep_dir, paste0(gse_id, "_ExpressionMatrix_Annotated.csv"))
  prep_sample <- file.path(prep_dir, paste0(gse_id, "_SampleMetadata.csv"))

  copied_expr <- FALSE
  copied_sample <- FALSE
  if (file.exists(legacy_expr) && !file.exists(prep_expr)) copied_expr <- file.copy(legacy_expr, prep_expr, overwrite = FALSE)
  if (file.exists(legacy_sample) && !file.exists(prep_sample)) copied_sample <- file.copy(legacy_sample, prep_sample, overwrite = FALSE)

  copy_status[[i]] <- data.frame(
    GSE_ID = gse_id,
    prep_expr_exists = file.exists(prep_expr),
    prep_sample_exists = file.exists(prep_sample),
    copied_expr = copied_expr,
    copied_sample = copied_sample,
    stringsAsFactors = FALSE
  )
}

copy_status_df <- dplyr::bind_rows(copy_status)
write.csv(copy_status_df, file.path(OUTPUT_DIR, "00_prepare_inputs_summary.csv"), row.names = FALSE)

TCGA_TAR <- file.path(INPUT_TCGA_DIR, "sarc_tcga_pub.tar.gz")
if (!file.exists(TCGA_TAR)) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Downloading TCGA-SARC cBioPortal study tar.gz")
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
  download.file(
    url = "https://datahub.assets.cbioportal.org/sarc_tcga_pub.tar.gz",
    destfile = TCGA_TAR,
    mode = "wb",
    quiet = FALSE
  )
}

msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "00_prepare_inputs finished")
cat(msg, "\n")
cat(msg, "\n", file = LOG_FILE, append = TRUE)
