rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(123)

PROJECT_DIR <- getwd()
INPUT_DIR <- file.path(PROJECT_DIR, "Input")
INPUT_GEO_PREP_DIR <- file.path(INPUT_DIR, "GEO_prepared")
INPUT_GEO_RAW_DIR <- file.path(INPUT_DIR, "GEO_raw")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Output")

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
PLOT_PNG_DIR <- file.path(OUTPUT_DIR, "plots_png")
PLOT_PDF_DIR <- file.path(OUTPUT_DIR, "plots_pdf")
CACHE_DIR <- file.path(OUTPUT_DIR, "cache")
GEO_CACHE <- file.path(CACHE_DIR, "geo_pheno")
LOG_DIR <- file.path(OUTPUT_DIR, "logs")
ADDON_TAB_DIR <- file.path(OUTPUT_DIR, "tables_addon")
ADDON_PLOT_PNG_DIR <- file.path(OUTPUT_DIR, "plots_png_addon")
ADDON_PLOT_PDF_DIR <- file.path(OUTPUT_DIR, "plots_pdf_addon")

dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_GEO_PREP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_GEO_RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_PNG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GEO_CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ADDON_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ADDON_PLOT_PNG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ADDON_PLOT_PDF_DIR, recursive = TRUE, showWarnings = FALSE)

PLOT_DPI <- 600L
LOG_FILE <- file.path(LOG_DIR, "discovery_main.log")

library(data.table)
library(Biobase)
library(GEOquery)
library(affy)
library(limma)
library(sva)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(R.utils)
library(matrixStats)
library(cluster)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(umap)
library(Rtsne)
library(msigdbr)
library(fgsea)
library(GSVA)
library(clusterProfiler)
library(enrichplot)
library(igraph)
library(ggraph)

data.table::setDTthreads(0L)

has_rcolorbrewer <- TRUE
has_umap <- TRUE
has_rtsne <- TRUE
has_ggrepel <- TRUE
has_msigdbr <- TRUE
has_fgsea <- TRUE
has_clusterprofiler <- TRUE
has_enrichplot <- TRUE
has_gsva <- TRUE
has_igraph <- TRUE
has_ggraph <- TRUE

gse_list <- c(
  "GSE71118","GSE21050","GSE16382","GSE34620","GSE142162",
  "GSE66533","GSE20196","GSE34800","GSE94321"
)

analysis_gse_tbl <- data.frame(
  gse_id = gse_list,
  cohort_role = c(
    "broad sarcoma tissue cohort",
    "broad sarcoma tissue cohort",
    "broad sarcoma tissue cohort",
    "Ewing sarcoma tissue cohort",
    "Ewing sarcoma tissue cohort",
    "rhabdomyosarcoma tissue cohort",
    "synovial sarcoma tissue cohort",
    "included in discovery; also excluded in sensitivity",
    "included in discovery after sample-level filtering; also excluded in sensitivity"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  analysis_gse_tbl,
  file.path(TABLE_DIR, "analysis_gse_manifest_9cohorts.csv"),
  row.names = FALSE
)

K_GRID               <- 2:6
K_FINAL_OVERRIDE     <- NA_integer_
BOOT_B               <- 30L
BOOT_PROP            <- 0.80
KSUPPORT_TOP_VAR     <- 3000L
KSUPPORT_N_PC        <- 30L
KSUPPORT_LOG_EVERY   <- 5L

LOW_VAR_QUANTILE     <- 0.05
DEG_FDR_CUTOFF       <- 0.05
DEG_LOGFC_CUTOFF     <- 1
TOPN_DEG_FOR_SIG     <- 50L
MIN_SAMPLES_PER_GSE  <- 5L

RUN_NETWORK          <- TRUE
NET_OVERALL_THR      <- 0.70
NET_RESIDUAL_THR     <- 0.30

RUN_UMAP             <- TRUE
RUN_TSNE             <- TRUE

COMP_MIN_N          <- 10L
COMP_TOP_N          <- 12L
COMP_SIM_B          <- 20000L
TOPN_MARKER_EXPORT  <- 100L
RUN_DOWN_ORA_SUPP   <- TRUE

RUN_SENSITIVITY_EXCLUDE_GSE <- TRUE
SENSITIVITY_EXCLUDE_GSE     <- c("GSE34800")

RUN_ARTICLE_ADDON              <- TRUE
K_GRID_ADDON                   <- 2:6
BOOT_B_ADDON                   <- 20L
BOOT_PROP_ADDON                <- 0.80
KSUPPORT_TOP_VAR_ADDON         <- 3000L
KSUPPORT_N_PC_ADDON            <- 30L
LEAVE_ONE_OUT_MIN_SAMPLES      <- 100L
SENSITIVITY_GSE_SET            <- c("GSE34800", "GSE94321")
TOP_MARKER_N_FOR_SCORE         <- 30L

write.csv(
  analysis_gse_tbl,
  file.path(TABLE_DIR, "analysis_gse_manifest_9cohorts.csv"),
  row.names = FALSE
)

sarcoma_keep_pattern <- paste(
  c(
    "sarcoma",
    "ewing",
    "synovial",
    "leiomyosarcoma",
    "leiomyosarcome",
    "liposarcoma",
    "rhabdomyosarcoma",
    "undifferentiated sarcoma",
    "pleomorphic sarcoma",
    "ups",
    "mpnst",
    "malignant peripheral nerve sheath",
    "osteosarcoma",
    "chondrosarcoma",
    "myxofibrosarcoma",
    "fibrosarcoma",
    "angiosarcoma",
    "epithelioid sarcoma",
    "alveolar soft part sarcoma",
    "clear cell sarcoma",
    "desmoplastic small round cell",
    "ddlps",
    "lms"
  ),
  collapse = "|"
)

exclude_non_sarcoma_pattern <- paste(
  c(
    "ovarian",
    "fallopian",
    "peritoneum",
    "kidney",
    "renal",
    "ccrcc",
    "brain",
    "meningi",
    "medulloblast",
    "astrocyt",
    "ependym",
    "craniopharyngi",
    "breast cancer",
    "colorectal cancer",
    "carcinoma",
    "normal brain",
    "adj.normal",
    "adj normal",
    "normal kidney",
    "endothelial",
    "lymphatic endothelial",
    "mesenchymal stem",
    "hmsc",
    "hpmsc",
    "wnt3a",
    "ews-fli"
  ),
  collapse = "|"
)

exclude_model_pattern <- paste(
  c(
    "cell line",
    "sorted",
    "side population",
    "\\bsp\\b",
    "\\bnsp\\b",
    "ssea-4",
    "xenograft",
    "tumorgraft",
    "orthotopic",
    "mouse",
    "murine",
    "pdx",
    "bone marrow"
  ),
  collapse = "|"
)

exclude_treatment_pattern <- paste(
  c(
    "\\btreated\\b",
    "\\bvehicle\\b",
    "\\bazacitidine\\b",
    "\\bbay1436032\\b",
    "\\bzstk\\b",
    "\\bpre-treatment\\b",
    "\\bpost-treatment\\b",
    "\\bpretreatment\\b",
    "\\bposttreatment\\b",
    "\\bpre treatment\\b",
    "\\bpost treatment\\b"
  ),
  collapse = "|"
)

exclude_normal_pattern <- paste(
  c(
    "normal",
    "adj normal",
    "adj.normal",
    "non-tumor",
    "non tumor"
  ),
  collapse = "|"
)

exprs_list <- list()
meta_list  <- list()
cohort_accounting <- list()

for (gse_id in gse_list) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Processing %s", gse_id))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)

  prep_dir <- file.path(INPUT_GEO_PREP_DIR, gse_id)
  raw_dir <- file.path(INPUT_GEO_RAW_DIR, gse_id)
  legacy_dir <- file.path(PROJECT_DIR, gse_id)
  exprs_file  <- file.path(prep_dir, paste0(gse_id, "_ExpressionMatrix_Annotated.csv"))
  sample_file <- file.path(prep_dir, paste0(gse_id, "_SampleMetadata.csv"))
  dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(exprs_file) && file.exists(file.path(legacy_dir, paste0(gse_id, "_ExpressionMatrix_Annotated.csv")))) {
    file.copy(file.path(legacy_dir, paste0(gse_id, "_ExpressionMatrix_Annotated.csv")), exprs_file, overwrite = FALSE)
  }
  if (!file.exists(sample_file) && file.exists(file.path(legacy_dir, paste0(gse_id, "_SampleMetadata.csv")))) {
    file.copy(file.path(legacy_dir, paste0(gse_id, "_SampleMetadata.csv")), sample_file, overwrite = FALSE)
  }

  loaded_from <- "reuse"
  mat <- NULL
  sample_df <- NULL

  if (file.exists(exprs_file)) {
    mat0 <- try(
      read.csv(exprs_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
      silent = TRUE
    )
    if (inherits(mat0, "try-error")) mat0 <- NULL
    if (!is.null(mat0) && nrow(mat0) > 0) {
      if ("Gene_Symbol" %in% colnames(mat0)) {
        rn <- mat0$Gene_Symbol
        mat0 <- mat0[, setdiff(colnames(mat0), "Gene_Symbol"), drop = FALSE]
        rownames(mat0) <- rn
      }
      mat <- as.matrix(mat0)
    }
  }

  if (file.exists(sample_file)) {
    sample_df <- try(
      read.csv(sample_file, check.names = FALSE, stringsAsFactors = FALSE),
      silent = TRUE
    )
    if (inherits(sample_df, "try-error")) sample_df <- NULL
  }

  if (is.null(mat) || ncol(mat) == 0) {
    loaded_from <- "download"
    msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("  Downloading raw CEL for %s", gse_id))
    cat(msg, "\n")
    cat(msg, "\n", file = LOG_FILE, append = TRUE)

    ok <- FALSE
    for (attempt in 1:5) {
      ok_try <- try({
        GEOquery::getGEOSuppFiles(GEO = gse_id, baseDir = INPUT_GEO_RAW_DIR)
        TRUE
      }, silent = TRUE)
      ok <- !inherits(ok_try, "try-error")
      if (!ok) Sys.sleep(5)
      if (ok) break
    }

    if (!ok) {
      msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("  FAILED download: %s", gse_id))
      cat(msg, "\n")
      cat(msg, "\n", file = LOG_FILE, append = TRUE)
      cohort_accounting[[gse_id]] <- data.frame(
        GSE_ID = gse_id, loaded_from = loaded_from, status = "download_fail",
        n_expr_samples = 0, n_genes = 0, n_samplefile_rows = ifelse(is.null(sample_df), NA, nrow(sample_df)),
        stringsAsFactors = FALSE
      )
      next
    }

    tarfiles <- list.files(raw_dir, pattern = "\\.tar$", full.names = TRUE, recursive = TRUE)
    if (length(tarfiles) > 0) {
      for (tf in tarfiles) {
        try(untar(tf, exdir = dirname(tf)), silent = TRUE)
      }
    }

    cels <- list.files(raw_dir, pattern = "\\.[cC][eE][lL]\\.gz$", full.names = TRUE, recursive = TRUE)
    if (length(cels) == 0) cels <- list.files(raw_dir, pattern = "\\.[cC][eE][lL]$", full.names = TRUE, recursive = TRUE)

    if (length(cels) == 0) {
      tar_probe <- list.files(raw_dir, full.names = TRUE, recursive = TRUE)
      tar_probe <- tar_probe[grepl("\\.tar$|\\.[cC][eE][lL](\\.gz)?$", tar_probe)]
      if (length(tar_probe) > 0) {
        msg2 <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("  Files found but CEL still unavailable for %s: %s", gse_id, paste(basename(tar_probe), collapse = "; ")))
        cat(msg2, "\n")
        cat(msg2, "\n", file = LOG_FILE, append = TRUE)
      }
      msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("  No CEL files for %s", gse_id))
      cat(msg, "\n")
      cat(msg, "\n", file = LOG_FILE, append = TRUE)
      cohort_accounting[[gse_id]] <- data.frame(
        GSE_ID = gse_id, loaded_from = loaded_from, status = "no_cel",
        n_expr_samples = 0, n_genes = 0, n_samplefile_rows = ifelse(is.null(sample_df), NA, nrow(sample_df)),
        stringsAsFactors = FALSE
      )
      next
    }
    cels2 <- character(length(cels))
    for (i_cel in seq_along(cels)) {
      f <- cels[i_cel]
      if (grepl("\\.gz$", f, ignore.case = TRUE)) {
        out <- sub("\\.gz$", "", f, ignore.case = TRUE)
        if (!file.exists(out)) try(R.utils::gunzip(f, remove = FALSE, overwrite = FALSE), silent = TRUE)
        cels2[i_cel] <- out
      } else {
        cels2[i_cel] <- f
      }
    }

    cels2 <- cels2[file.exists(cels2)]
    if (length(cels2) == 0) {
      msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("  Gunzip/readable CEL absent for %s", gse_id))
      cat(msg, "\n")
      cat(msg, "\n", file = LOG_FILE, append = TRUE)
      next
    }

    eset <- try({
      affy_raw <- affy::ReadAffy(filenames = cels2, cdfname = "hgu133plus2")
      affy::rma(affy_raw)
    }, silent = TRUE)
    if (inherits(eset, "try-error")) eset <- NULL

    if (is.null(eset)) {
      msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("  RMA failed for %s", gse_id))
      cat(msg, "\n")
      cat(msg, "\n", file = LOG_FILE, append = TRUE)
      next
    }

    probe_ids <- try(Biobase::featureNames(eset), silent = TRUE)
    if (inherits(probe_ids, "try-error")) probe_ids <- NULL
    if (is.null(probe_ids)) next

    anno <- try(
      AnnotationDbi::select(
        hgu133plus2.db,
        keys = probe_ids, columns = "SYMBOL", keytype = "PROBEID"
      ),
      silent = TRUE
    )
    if (inherits(anno, "try-error")) anno <- NULL

    ex <- Biobase::exprs(eset)
    sym <- if (!is.null(anno)) anno$SYMBOL[match(probe_ids, anno$PROBEID)] else NA_character_
    keep <- !is.na(sym) & sym != ""
    ex <- ex[keep, , drop = FALSE]
    sym <- sym[keep]

    exprs_agg <- try(
      aggregate(ex, by = list(Gene_Symbol = sym), FUN = mean, na.rm = TRUE),
      silent = TRUE
    )
    if (inherits(exprs_agg, "try-error")) exprs_agg <- NULL

    if (is.null(exprs_agg) || nrow(exprs_agg) == 0) next

    write.csv(exprs_agg, exprs_file, row.names = FALSE)
    sample_df <- try(as.data.frame(Biobase::pData(eset), stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(sample_df, "try-error")) sample_df <- NULL
    if (is.null(sample_df)) sample_df <- data.frame(sample_row_order = seq_len(ncol(exprs_agg) - 1))
    write.csv(sample_df, sample_file, row.names = FALSE)

    mat <- as.matrix(exprs_agg[, -1, drop = FALSE])
    rownames(mat) <- exprs_agg$Gene_Symbol
  }

  if (is.null(mat) || ncol(mat) == 0) {
    cohort_accounting[[gse_id]] <- data.frame(
      GSE_ID = gse_id, loaded_from = loaded_from, status = "expr_missing",
      n_expr_samples = 0, n_genes = 0, n_samplefile_rows = ifelse(is.null(sample_df), NA, nrow(sample_df)),
      stringsAsFactors = FALSE
    )
    next
  }

  orig_colnames <- colnames(mat)
  colnames(mat) <- paste(gse_id, orig_colnames, sep = "_")
  exprs_list[[gse_id]] <- mat

  geo_pheno <- NULL
  rds_file <- file.path(GEO_CACHE, paste0(gse_id, "_geo_pheno.rds"))
  if (file.exists(rds_file)) {
    geo_pheno <- try(readRDS(rds_file), silent = TRUE)
    if (inherits(geo_pheno, "try-error")) geo_pheno <- NULL
  } else {
    msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("GEO pheno fetch: %s", gse_id))
    cat(msg, "\n")
    cat(msg, "\n", file = LOG_FILE, append = TRUE)

    gse_obj <- try(GEOquery::getGEO(gse_id, GSEMatrix = TRUE), silent = TRUE)
    if (inherits(gse_obj, "try-error")) gse_obj <- NULL
    if (!is.null(gse_obj)) {
      if (length(gse_obj) > 1) gse_obj <- gse_obj[[1]] else gse_obj <- gse_obj[[1]]
      pheno <- try(Biobase::pData(gse_obj), silent = TRUE)
      if (inherits(pheno, "try-error")) pheno <- NULL

      if (!is.null(pheno) && nrow(pheno) > 0) {
        char_cols <- grep("^characteristics", names(pheno), value = TRUE, ignore.case = TRUE)
        char_df <- NULL

        if (length(char_cols) > 0) {
          out_list <- vector("list", nrow(pheno))
          for (i in seq_len(nrow(pheno))) {
            kv <- list()
            vals <- unlist(pheno[i, char_cols, drop = TRUE], use.names = FALSE)
            vals <- vals[!is.na(vals) & vals != ""]
            if (length(vals) > 0) {
              for (v in vals) {
                pieces <- trimws(unlist(strsplit(as.character(v), ";", fixed = TRUE)))
                pieces <- pieces[pieces != ""]
                if (length(pieces) > 0) {
                  for (p in pieces) {
                    if (grepl(":", p, fixed = TRUE)) {
                      sp <- strsplit(p, ":", fixed = TRUE)[[1]]
                      k <- trimws(sp[1])
                      vv <- trimws(paste(sp[-1], collapse = ":"))
                      if (k != "" && vv != "") {
                        k2 <- tolower(gsub("[^A-Za-z0-9]+", "_", k))
                        kv[[k2]] <- vv
                      }
                    }
                  }
                }
              }
            }
            out_list[[i]] <- as.data.frame(kv, stringsAsFactors = FALSE)
          }

          char_df <- dplyr::bind_rows(out_list)
          if (!is.null(char_df) && nrow(char_df) == 0) char_df <- NULL
        }

        geo_pheno <- cbind(
          data.frame(geo_row_order = seq_len(nrow(pheno)), stringsAsFactors = FALSE),
          pheno,
          char_df
        )
        saveRDS(geo_pheno, rds_file)
      }
    }
  }

  expr_names <- as.character(orig_colnames)
  n_expr <- length(expr_names)

  base_df <- data.frame(
    SampleID = paste(gse_id, expr_names, sep = "_"),
    expr_col_original = expr_names,
    expr_norm = as.character(expr_names),
    GSE_ID = gse_id,
    expr_order = seq_len(n_expr),
    stringsAsFactors = FALSE
  )
  base_df$expr_norm <- trimws(base_df$expr_norm)
  base_df$expr_norm <- tolower(base_df$expr_norm)
  base_df$expr_norm <- gsub("\\.cel(\\.gz)?$", "", base_df$expr_norm, ignore.case = TRUE)
  base_df$expr_norm <- gsub("[^a-z0-9]+", "", base_df$expr_norm)

  meta_gse <- base_df
  match_summary <- c(sample_file_match = "none", geo_match = "none")

  if (!is.null(sample_df) && nrow(sample_df) > 0) {
    sample_df <- as.data.frame(sample_df, stringsAsFactors = FALSE)
    sample_df$sample_row_order <- seq_len(nrow(sample_df))

    cand_sample <- intersect(
      names(sample_df),
      c("SampleID", "sampleid", "sample_id", "geo_accession", "title", "sample_name",
        "Sample_Name", "filename", "file_name", "cel_file", "array_data_file",
        "source_name_ch1", "gsm", "GSM", "supplementary_file")
    )

    matched <- NULL
    if (length(cand_sample) > 0) {
      best_hits <- -1L
      best_col <- NULL

      for (cc in cand_sample) {
        x <- as.character(sample_df[[cc]])
        x <- trimws(x)
        x <- tolower(x)
        x <- gsub("\\.cel(\\.gz)?$", "", x, ignore.case = TRUE)
        x <- gsub("[^a-z0-9]+", "", x)

        hits <- sum(base_df$expr_norm %in% x, na.rm = TRUE)
        if (hits > best_hits) {
          best_hits <- hits
          best_col <- cc
        }
      }

      if (!is.null(best_col) && best_hits >= max(3L, floor(0.2 * n_expr))) {
        x <- as.character(sample_df[[best_col]])
        x <- trimws(x)
        x <- tolower(x)
        x <- gsub("\\.cel(\\.gz)?$", "", x, ignore.case = TRUE)
        x <- gsub("[^a-z0-9]+", "", x)

        idx <- match(base_df$expr_norm, x)
        if (sum(!is.na(idx)) > 0) {
          matched <- sample_df[idx, , drop = FALSE]
          match_summary["sample_file_match"] <- paste0("key:", best_col)
        }
      }
    }

    if (is.null(matched) && nrow(sample_df) == n_expr) {
      matched <- sample_df
      match_summary["sample_file_match"] <- "row_order"
    }

    if (!is.null(matched)) {
      names(matched) <- make.names(names(matched), unique = TRUE)
      keep_cols <- setdiff(names(matched), c("SampleID", "sampleid", "sample_id"))
      meta_gse <- cbind(meta_gse, matched[, keep_cols, drop = FALSE])
    }
  }

  if (!is.null(geo_pheno) && nrow(geo_pheno) > 0) {
    geo_pheno <- as.data.frame(geo_pheno, stringsAsFactors = FALSE)
    geo_pheno$geo_row_order2 <- seq_len(nrow(geo_pheno))

    cand_geo <- intersect(
      names(geo_pheno),
      c("geo_accession", "title", "source_name_ch1", "sample_name",
        "SampleID", "sampleid", "sample_id", "supplementary_file")
    )

    matched_geo <- NULL
    if (length(cand_geo) > 0) {
      key_pool <- list(
        expr_col_original = meta_gse$expr_col_original,
        sample_title = if ("title" %in% names(meta_gse)) meta_gse$title else NULL,
        sample_geo_accession = if ("geo_accession" %in% names(meta_gse)) meta_gse$geo_accession else NULL
      )

      best_hits <- -1L
      best_col <- NULL
      best_idx <- NULL

      for (kk in names(key_pool)) {
        if (is.null(key_pool[[kk]])) next

        lhs <- as.character(key_pool[[kk]])
        lhs <- trimws(lhs)
        lhs <- tolower(lhs)
        lhs <- gsub("\\.cel(\\.gz)?$", "", lhs, ignore.case = TRUE)
        lhs <- gsub("[^a-z0-9]+", "", lhs)

        for (cc in cand_geo) {
          rhs <- as.character(geo_pheno[[cc]])
          rhs <- trimws(rhs)
          rhs <- tolower(rhs)
          rhs <- gsub("\\.cel(\\.gz)?$", "", rhs, ignore.case = TRUE)
          rhs <- gsub("[^a-z0-9]+", "", rhs)

          idx <- match(lhs, rhs)
          hits <- sum(!is.na(idx))
          if (hits > best_hits) {
            best_hits <- hits
            best_col <- paste0(kk, "->", cc)
            best_idx <- idx
          }
        }
      }

      if (!is.null(best_idx) && best_hits >= max(3L, floor(0.2 * n_expr))) {
        matched_geo <- geo_pheno[best_idx, , drop = FALSE]
        match_summary["geo_match"] <- paste0("key:", best_col)
      }
    }

    if (is.null(matched_geo) && nrow(geo_pheno) == n_expr) {
      matched_geo <- geo_pheno
      match_summary["geo_match"] <- "row_order"
    }

    if (!is.null(matched_geo)) {
      names(matched_geo) <- make.names(names(matched_geo), unique = TRUE)
      add_cols <- setdiff(names(matched_geo), names(meta_gse))
      if (length(add_cols) > 0) {
        meta_gse <- cbind(meta_gse, matched_geo[, add_cols, drop = FALSE])
      }
    }
  }

  meta_gse$sample_file_match <- match_summary["sample_file_match"]
  meta_gse$geo_match <- match_summary["geo_match"]

  meta_list[[gse_id]] <- meta_gse

  cohort_accounting[[gse_id]] <- data.frame(
    GSE_ID = gse_id,
    loaded_from = loaded_from,
    status = "ok",
    n_expr_samples = ncol(mat),
    n_genes = nrow(mat),
    n_samplefile_rows = ifelse(is.null(sample_df), NA, nrow(sample_df)),
    n_geo_pheno_rows = ifelse(is.null(geo_pheno), NA, nrow(geo_pheno)),
    sample_file_match = unique(meta_gse$sample_file_match)[1],
    geo_match = unique(meta_gse$geo_match)[1],
    stringsAsFactors = FALSE
  )
}

cohort_accounting_df <- dplyr::bind_rows(cohort_accounting)
write.csv(cohort_accounting_df, file.path(TABLE_DIR, "cohort_accounting_by_GSE.csv"), row.names = FALSE)

stopifnot(length(exprs_list) > 0)

meta_list_chr <- vector("list", length(meta_list))
names(meta_list_chr) <- names(meta_list)
for (ii_meta in seq_along(meta_list)) {
  tmp_meta <- meta_list[[ii_meta]]
  if (!is.null(tmp_meta) && nrow(tmp_meta) > 0) {
    for (jj_meta in seq_len(ncol(tmp_meta))) tmp_meta[[jj_meta]] <- as.character(tmp_meta[[jj_meta]])
  }
  meta_list_chr[[ii_meta]] <- tmp_meta
}
meta_df <- dplyr::bind_rows(meta_list_chr)
meta_df <- dplyr::mutate(meta_df, dplyr::across(dplyr::everything(), as.character))
nm1 <- trimws(tolower(names(meta_df)))
nm1 <- gsub("[^a-z0-9_]+", "_", nm1)
nm1 <- gsub("_+", "_", nm1)
nm1 <- gsub("^_|_$", "", nm1)
names(meta_df) <- make.unique(nm1, sep = "__dup")

meta_df$sampleid <- ifelse(is.na(meta_df$sampleid), "", meta_df$sampleid)
meta_df <- meta_df[!duplicated(meta_df$sampleid), , drop = FALSE]

all_sample_ids <- unlist(lapply(exprs_list, colnames), use.names = FALSE)

sample_master <- data.frame(
  sampleid = all_sample_ids,
  gse_id   = sub("_.*$", "", all_sample_ids),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(meta_df, by = c("sampleid", "gse_id"))

nm1 <- trimws(tolower(names(sample_master)))
nm1 <- gsub("[^a-z0-9_]+", "_", nm1)
nm1 <- gsub("_+", "_", nm1)
nm1 <- gsub("^_|_$", "", nm1)
names(sample_master) <- make.unique(nm1, sep = "__dup")
sample_master <- dplyr::left_join(sample_master, analysis_gse_tbl, by = "gse_id")

if (!"diagnosis" %in% names(sample_master)) sample_master$diagnosis <- NA_character_
if (!"diagnosis_ch1" %in% names(sample_master)) sample_master$diagnosis_ch1 <- NA_character_
if (!"tumor_type" %in% names(sample_master)) sample_master$tumor_type <- NA_character_
if (!"tumor_type_ch1" %in% names(sample_master)) sample_master$tumor_type_ch1 <- NA_character_
if (!"cancer" %in% names(sample_master)) sample_master$cancer <- NA_character_
if (!"cancer_ch1" %in% names(sample_master)) sample_master$cancer_ch1 <- NA_character_
if (!"source_name_ch1" %in% names(sample_master)) sample_master$source_name_ch1 <- NA_character_
if (!"title" %in% names(sample_master)) sample_master$title <- NA_character_
if (!"subtype" %in% names(sample_master)) sample_master$subtype <- NA_character_
if (!"subtype_ch1" %in% names(sample_master)) sample_master$subtype_ch1 <- NA_character_
if (!"cinsarc_group" %in% names(sample_master)) sample_master$cinsarc_group <- NA_character_
if (!"cinsarc_group_ch1" %in% names(sample_master)) sample_master$cinsarc_group_ch1 <- NA_character_
if (!"molecular_subtype" %in% names(sample_master)) sample_master$molecular_subtype <- NA_character_
if (!"subclassification" %in% names(sample_master)) sample_master$subclassification <- NA_character_
if (!"histological_subtype" %in% names(sample_master)) sample_master$histological_subtype <- NA_character_
if (!"reclassified_histological_subtype" %in% names(sample_master)) sample_master$reclassified_histological_subtype <- NA_character_
if (!"histology" %in% names(sample_master)) sample_master$histology <- NA_character_
if (!"histology_ch1" %in% names(sample_master)) sample_master$histology_ch1 <- NA_character_
if (!"histologic_subtype" %in% names(sample_master)) sample_master$histologic_subtype <- NA_character_
if (!"histologic_type" %in% names(sample_master)) sample_master$histologic_type <- NA_character_
if (!"tissue" %in% names(sample_master)) sample_master$tissue <- NA_character_
if (!"tissue_ch1" %in% names(sample_master)) sample_master$tissue_ch1 <- NA_character_
if (!"treatment" %in% names(sample_master)) sample_master$treatment <- NA_character_
if (!"treatment_ch1" %in% names(sample_master)) sample_master$treatment_ch1 <- NA_character_
if (!"treatment_protocol_ch1" %in% names(sample_master)) sample_master$treatment_protocol_ch1 <- NA_character_
if (!"primary_treatment" %in% names(sample_master)) sample_master$primary_treatment <- NA_character_
if (!"primary_treatment_ch1" %in% names(sample_master)) sample_master$primary_treatment_ch1 <- NA_character_
if (!"grade" %in% names(sample_master)) sample_master$grade <- NA_character_
if (!"grade_ch1" %in% names(sample_master)) sample_master$grade_ch1 <- NA_character_
if (!"tumor_grade" %in% names(sample_master)) sample_master$tumor_grade <- NA_character_
if (!"stage" %in% names(sample_master)) sample_master$stage <- NA_character_
if (!"stage_ch1" %in% names(sample_master)) sample_master$stage_ch1 <- NA_character_
if (!"tumor_stage" %in% names(sample_master)) sample_master$tumor_stage <- NA_character_
if (!"organism2" %in% names(sample_master)) sample_master$organism2 <- NA_character_
if (!"organism_ch1" %in% names(sample_master)) sample_master$organism_ch1 <- NA_character_
if (!"organism" %in% names(sample_master)) sample_master$organism <- NA_character_

sample_master$diagnosis_h <- as.character(sample_master$diagnosis)
sample_master$diagnosis_h[trimws(sample_master$diagnosis_h) == ""] <- NA_character_
sample_master$diagnosis_h[is.na(sample_master$diagnosis_h)] <- as.character(sample_master$diagnosis_ch1[is.na(sample_master$diagnosis_h)])
sample_master$diagnosis_h[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""] <- as.character(sample_master$tumor_type[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""])
sample_master$diagnosis_h[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""] <- as.character(sample_master$tumor_type_ch1[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""])
sample_master$diagnosis_h[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""] <- as.character(sample_master$cancer[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""])
sample_master$diagnosis_h[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""] <- as.character(sample_master$cancer_ch1[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""])
sample_master$diagnosis_h[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""] <- as.character(sample_master$source_name_ch1[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""])
sample_master$diagnosis_h[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""] <- as.character(sample_master$title[is.na(sample_master$diagnosis_h) | trimws(sample_master$diagnosis_h) == ""])

sample_master$subtype_h <- as.character(sample_master$subtype)
sample_master$subtype_h[trimws(sample_master$subtype_h) == ""] <- NA_character_
sample_master$subtype_h[is.na(sample_master$subtype_h)] <- as.character(sample_master$subtype_ch1[is.na(sample_master$subtype_h)])
sample_master$subtype_h[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""] <- as.character(sample_master$cinsarc_group[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""])
sample_master$subtype_h[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""] <- as.character(sample_master$cinsarc_group_ch1[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""])
sample_master$subtype_h[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""] <- as.character(sample_master$molecular_subtype[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""])
sample_master$subtype_h[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""] <- as.character(sample_master$subclassification[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""])
sample_master$subtype_h[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""] <- as.character(sample_master$histological_subtype[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""])
sample_master$subtype_h[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""] <- as.character(sample_master$reclassified_histological_subtype[is.na(sample_master$subtype_h) | trimws(sample_master$subtype_h) == ""])

sample_master$histology_h <- as.character(sample_master$histology)
sample_master$histology_h[trimws(sample_master$histology_h) == ""] <- NA_character_
sample_master$histology_h[is.na(sample_master$histology_h)] <- as.character(sample_master$histology_ch1[is.na(sample_master$histology_h)])
sample_master$histology_h[is.na(sample_master$histology_h) | trimws(sample_master$histology_h) == ""] <- as.character(sample_master$histologic_subtype[is.na(sample_master$histology_h) | trimws(sample_master$histology_h) == ""])
sample_master$histology_h[is.na(sample_master$histology_h) | trimws(sample_master$histology_h) == ""] <- as.character(sample_master$histologic_type[is.na(sample_master$histology_h) | trimws(sample_master$histology_h) == ""])

sample_master$tissue_h <- as.character(sample_master$tissue)
sample_master$tissue_h[trimws(sample_master$tissue_h) == ""] <- NA_character_
sample_master$tissue_h[is.na(sample_master$tissue_h)] <- as.character(sample_master$tissue_ch1[is.na(sample_master$tissue_h)])
sample_master$tissue_h[is.na(sample_master$tissue_h) | trimws(sample_master$tissue_h) == ""] <- as.character(sample_master$source_name_ch1[is.na(sample_master$tissue_h) | trimws(sample_master$tissue_h) == ""])
sample_master$tissue_h[is.na(sample_master$tissue_h) | trimws(sample_master$tissue_h) == ""] <- as.character(sample_master$title[is.na(sample_master$tissue_h) | trimws(sample_master$tissue_h) == ""])

sample_master$treatment_h <- as.character(sample_master$treatment)
sample_master$treatment_h[trimws(sample_master$treatment_h) == ""] <- NA_character_
sample_master$treatment_h[is.na(sample_master$treatment_h)] <- as.character(sample_master$treatment_ch1[is.na(sample_master$treatment_h)])
sample_master$treatment_h[is.na(sample_master$treatment_h) | trimws(sample_master$treatment_h) == ""] <- as.character(sample_master$treatment_protocol_ch1[is.na(sample_master$treatment_h) | trimws(sample_master$treatment_h) == ""])
sample_master$treatment_h[is.na(sample_master$treatment_h) | trimws(sample_master$treatment_h) == ""] <- as.character(sample_master$primary_treatment[is.na(sample_master$treatment_h) | trimws(sample_master$treatment_h) == ""])
sample_master$treatment_h[is.na(sample_master$treatment_h) | trimws(sample_master$treatment_h) == ""] <- as.character(sample_master$primary_treatment_ch1[is.na(sample_master$treatment_h) | trimws(sample_master$treatment_h) == ""])

sample_master$grade_h <- as.character(sample_master$grade)
sample_master$grade_h[trimws(sample_master$grade_h) == ""] <- NA_character_
sample_master$grade_h[is.na(sample_master$grade_h)] <- as.character(sample_master$grade_ch1[is.na(sample_master$grade_h)])
sample_master$grade_h[is.na(sample_master$grade_h) | trimws(sample_master$grade_h) == ""] <- as.character(sample_master$tumor_grade[is.na(sample_master$grade_h) | trimws(sample_master$grade_h) == ""])

sample_master$stage_h <- as.character(sample_master$stage)
sample_master$stage_h[trimws(sample_master$stage_h) == ""] <- NA_character_
sample_master$stage_h[is.na(sample_master$stage_h)] <- as.character(sample_master$stage_ch1[is.na(sample_master$stage_h)])
sample_master$stage_h[is.na(sample_master$stage_h) | trimws(sample_master$stage_h) == ""] <- as.character(sample_master$tumor_stage[is.na(sample_master$stage_h) | trimws(sample_master$stage_h) == ""])

sample_master$organism2 <- as.character(sample_master$organism2)
sample_master$organism2[trimws(sample_master$organism2) == ""] <- NA_character_
sample_master$organism2[is.na(sample_master$organism2)] <- as.character(sample_master$organism_ch1[is.na(sample_master$organism2)])
sample_master$organism2[is.na(sample_master$organism2) | trimws(sample_master$organism2) == ""] <- as.character(sample_master$organism[is.na(sample_master$organism2) | trimws(sample_master$organism2) == ""])

sample_master$diagnosis_h <- gsub("\\s+", " ", trimws(sample_master$diagnosis_h))
sample_master$subtype_h <- gsub("\\s+", " ", trimws(sample_master$subtype_h))
sample_master$histology_h <- gsub("\\s+", " ", trimws(sample_master$histology_h))
sample_master$tissue_h <- gsub("\\s+", " ", trimws(sample_master$tissue_h))
sample_master$treatment_h <- gsub("\\s+", " ", trimws(sample_master$treatment_h))
sample_master$grade_h <- gsub("\\s+", " ", trimws(sample_master$grade_h))
sample_master$stage_h <- gsub("\\s+", " ", trimws(sample_master$stage_h))
sample_master$organism2 <- gsub("\\s+", " ", trimws(sample_master$organism2))
sample_master$diagnosis_h[sample_master$diagnosis_h == ""] <- NA_character_
sample_master$subtype_h[sample_master$subtype_h == ""] <- NA_character_
sample_master$histology_h[sample_master$histology_h == ""] <- NA_character_
sample_master$tissue_h[sample_master$tissue_h == ""] <- NA_character_
sample_master$treatment_h[sample_master$treatment_h == ""] <- NA_character_
sample_master$grade_h[sample_master$grade_h == ""] <- NA_character_
sample_master$stage_h[sample_master$stage_h == ""] <- NA_character_
sample_master$organism2[sample_master$organism2 == ""] <- NA_character_

text_cols <- intersect(
  c(
    "title", "source_name_ch1", "diagnosis_h", "histology_h", "subtype_h",
    "tissue_h", "treatment_h", "grade_h", "stage_h"
  ),
  names(sample_master)
)

sample_text <- rep("", nrow(sample_master))
for (i_text in seq_len(nrow(sample_master))) {
  z <- unlist(sample_master[i_text, text_cols, drop = FALSE], use.names = FALSE)
  z <- z[!is.na(z) & z != ""]
  sample_text[i_text] <- paste(z, collapse = " | ")
}
sample_text <- tolower(sample_text)

sample_master$curation_text <- sample_text

dx_txt <- sample_master$curation_text
sample_master$diagnosis_group <- NA_character_
sample_master$diagnosis_group[grepl("leiomyosarcoma|leiomyosarcome|\\blms\\b", dx_txt)] <- "Leiomyosarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("ewing", dx_txt)] <- "Ewing sarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("synovial", dx_txt)] <- "Synovial sarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("rhabdomyosarcoma", dx_txt)] <- "Rhabdomyosarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("liposarcoma.*dediff|dediff.*liposarcoma|dedifferentiated liposarcoma", dx_txt)] <- "Dedifferentiated liposarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("liposarcoma", dx_txt)] <- "Liposarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("undifferentiated sarcoma|pleomorphic sarcoma|\\bups\\b", dx_txt)] <- "Undifferentiated/pleomorphic sarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("mpnst|malignant peripheral nerve sheath", dx_txt)] <- "MPNST"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("osteosarcoma", dx_txt)] <- "Osteosarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("chondrosarcoma", dx_txt)] <- "Chondrosarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("sarcoma|ewing|rhabdo|synovial|lipo|mpnst|osteosarcoma|chondrosarcoma", dx_txt)] <- "Other sarcoma"
x <- as.character(sample_master$diagnosis_group)
x <- trimws(x)
x[x == ""] <- NA_character_
x <- gsub("\\s+", " ", x)
sample_master$diagnosis_group <- x

tx_txt <- tolower(paste(sample_master$tissue_h, sample_master$source_name_ch1, sample_master$title, sep = " | "))
sample_master$tissue_group <- NA_character_
sample_master$tissue_group[grepl("extremit", tx_txt)] <- "Extremities"
sample_master$tissue_group[is.na(sample_master$tissue_group) & grepl("internal trunk", tx_txt)] <- "Internal trunk"
sample_master$tissue_group[is.na(sample_master$tissue_group) & grepl("trunk wall", tx_txt)] <- "Trunk wall"
sample_master$tissue_group[is.na(sample_master$tissue_group) & grepl("pelvis|pelvic", tx_txt)] <- "Pelvis"
sample_master$tissue_group[is.na(sample_master$tissue_group) & grepl("abdomen|abdominal|retroperitone", tx_txt)] <- "Abdomen/retroperitoneum"
sample_master$tissue_group[is.na(sample_master$tissue_group) & grepl("head|neck", tx_txt)] <- "Head/neck"
sample_master$tissue_group[is.na(sample_master$tissue_group) & grepl("bone|osseous", tx_txt)] <- "Bone"
sample_master$tissue_group[is.na(sample_master$tissue_group) & grepl("soft tissue", tx_txt)] <- "Soft tissue NOS"
x <- as.character(sample_master$tissue_group)
x <- trimws(x)
x[x == ""] <- NA_character_
x <- gsub("\\s+", " ", x)
sample_master$tissue_group <- x

sample_master$is_human <- ifelse(
  is.na(sample_master$organism2) | sample_master$organism2 == "",
  TRUE,
  grepl("homo sapiens|human", tolower(sample_master$organism2))
)

sample_master$hit_keep_sarcoma <- grepl(sarcoma_keep_pattern, sample_master$curation_text, ignore.case = TRUE)
sample_master$hit_exclude_nonsarcoma <- grepl(exclude_non_sarcoma_pattern, sample_master$curation_text, ignore.case = TRUE)
sample_master$hit_exclude_model <- grepl(exclude_model_pattern, sample_master$curation_text, ignore.case = TRUE)
sample_master$hit_exclude_treatment <- grepl(exclude_treatment_pattern, sample_master$curation_text, ignore.case = TRUE)
sample_master$hit_exclude_normal <- grepl(exclude_normal_pattern, sample_master$curation_text, ignore.case = TRUE)

sample_master$sample_decision <- "DROP"
sample_master$sample_reason <- ""

always_keep_gse <- c(
  "GSE71118","GSE21050","GSE16382","GSE34620",
  "GSE142162","GSE66533","GSE20196","GSE34800"
)
needs_positive_keyword_gse <- c("GSE94321")

idx_eval <- sample_master$gse_id %in% analysis_gse_tbl$gse_id
positive_keyword_required <- sample_master$gse_id[idx_eval] %in% needs_positive_keyword_gse

sample_master$sample_decision[idx_eval] <- ifelse(
  sample_master$is_human[idx_eval] &
    !sample_master$hit_exclude_nonsarcoma[idx_eval] &
    !sample_master$hit_exclude_model[idx_eval] &
    !sample_master$hit_exclude_treatment[idx_eval] &
    !sample_master$hit_exclude_normal[idx_eval] &
    (
      sample_master$gse_id[idx_eval] %in% always_keep_gse |
        (positive_keyword_required & sample_master$hit_keep_sarcoma[idx_eval])
    ),
  "KEEP", "DROP"
)

reason_eval <- data.frame(
  nonhuman = ifelse(!sample_master$is_human[idx_eval], "nonhuman_or_unknown", ""),
  nonsarcoma = ifelse(sample_master$hit_exclude_nonsarcoma[idx_eval], "nonsarcoma_keyword", ""),
  model = ifelse(sample_master$hit_exclude_model[idx_eval], "model_or_sorted_keyword", ""),
  treatment = ifelse(sample_master$hit_exclude_treatment[idx_eval], "treatment_keyword", ""),
  normal = ifelse(sample_master$hit_exclude_normal[idx_eval], "normal_keyword", ""),
  no_positive = ifelse(
    positive_keyword_required & !sample_master$hit_keep_sarcoma[idx_eval],
    "no_positive_sarcoma_keyword_in_keyword_required_cohort",
    ""
  ),
  stringsAsFactors = FALSE
)

sample_reason_eval <- rep("", nrow(reason_eval))
for (i_reason in seq_len(nrow(reason_eval))) {
  z <- unlist(reason_eval[i_reason, , drop = FALSE], use.names = FALSE)
  sample_reason_eval[i_reason] <- paste(z[z != ""], collapse = ";")
}

sample_master$sample_reason[idx_eval] <- ifelse(
  sample_master$sample_decision[idx_eval] == "KEEP",
  "passed_conservative_sample_filter",
  ifelse(sample_reason_eval == "", "filtered_by_conservative_sample_rules", sample_reason_eval)
)

sample_master$sample_reason <- gsub("(^;|;;+|;$)", "", sample_master$sample_reason)
sample_master$sample_reason <- gsub(";+", ";", sample_master$sample_reason)

gse_post <- sample_master %>%
  dplyr::filter(gse_id %in% analysis_gse_tbl$gse_id) %>%
  dplyr::group_by(gse_id) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_keep_initial = sum(sample_decision == "KEEP", na.rm = TRUE),
    .groups = "drop"
  )

small_gse <- gse_post$gse_id[gse_post$n_keep_initial > 0 & gse_post$n_keep_initial < MIN_SAMPLES_PER_GSE]
if (length(small_gse) > 0) {
  sample_master$sample_decision[sample_master$gse_id %in% small_gse & sample_master$sample_decision == "KEEP"] <- "DROP"
  sample_master$sample_reason[sample_master$gse_id %in% small_gse & grepl("^passed_conservative_sample_filter$", sample_master$sample_reason)] <-
    paste0("dropped_posthoc_small_gse_lt_", MIN_SAMPLES_PER_GSE)
}

gse_post2 <- sample_master %>%
  dplyr::filter(gse_id %in% analysis_gse_tbl$gse_id) %>%
  dplyr::group_by(gse_id, cohort_role) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_keep = sum(sample_decision == "KEEP", na.rm = TRUE),
    n_drop = sum(sample_decision == "DROP", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_keep), gse_id)

write.csv(sample_master, file.path(TABLE_DIR, "sample_master_metadata_with_curation.csv"), row.names = FALSE)
write.csv(gse_post2, file.path(TABLE_DIR, "gse_post_curation_summary.csv"), row.names = FALSE)

write.csv(
  subset(gse_post2, n_keep > 0),
  file.path(TABLE_DIR, "gse_final_discovery_cohort_summary.csv"),
  row.names = FALSE
)

keep_samples <- sample_master$sampleid[sample_master$sample_decision == "KEEP"]
if (length(keep_samples) < 30) {
  stop("Too few curated samples remained. Review the sample-filtering rules before proceeding.")
}
msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Final discovery cohort samples kept: %d", length(keep_samples)))
cat(msg, "\n")
cat(msg, "\n", file = LOG_FILE, append = TRUE)

exprs_keep_list <- list()
gse_keep_map <- split(keep_samples, sub("_.*$", "", keep_samples))

for (gse_id in names(exprs_list)) {
  want <- intersect(colnames(exprs_list[[gse_id]]), gse_keep_map[[gse_id]])
  if (length(want) > 0) {
    exprs_keep_list[[gse_id]] <- exprs_list[[gse_id]][, want, drop = FALSE]
  }
}

if (length(exprs_keep_list) == 0) stop("No curated expression matrices remained.")
exprs_all <- do.call(cbind, exprs_keep_list)
batch_info <- rep(names(exprs_keep_list), times = vapply(exprs_keep_list, ncol, integer(1)))
stopifnot(length(batch_info) == ncol(exprs_all))

pre_filter <- data.frame(
  step = "raw_curated_merged",
  n_genes = nrow(exprs_all),
  n_samples = ncol(exprs_all),
  stringsAsFactors = FALSE
)

nz <- rowSums(exprs_all != 0, na.rm = TRUE) > 0
exprs_all <- exprs_all[nz, , drop = FALSE]

rv <- matrixStats::rowVars(exprs_all)
keep_var <- rv > stats::quantile(rv, LOW_VAR_QUANTILE, na.rm = TRUE)
exprs_all <- exprs_all[keep_var, , drop = FALSE]

post_filter <- data.frame(
  step = c("post_zero_filter", "post_low_variance_filter"),
  n_genes = c(sum(nz), nrow(exprs_all)),
  n_samples = c(ncol(exprs_all), ncol(exprs_all)),
  stringsAsFactors = FALSE
)

write.csv(rbind(pre_filter, post_filter), file.path(TABLE_DIR, "expression_filtering_summary.csv"), row.names = FALSE)
write.csv(data.frame(SampleID = colnames(exprs_all), Batch = batch_info), file.path(TABLE_DIR, "sample_batch_info.csv"), row.names = FALSE)

if (length(unique(batch_info)) >= 2) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Running ComBat on curated matrix"))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
  exprs_combat <- try(
    sva::ComBat(dat = exprs_all, batch = batch_info, par.prior = TRUE, prior.plots = FALSE),
    silent = TRUE
  )
  if (inherits(exprs_combat, "try-error")) {
    msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "ComBat failed")
    cat(msg, "\n")
    cat(msg, "\n", file = LOG_FILE, append = TRUE)
    exprs_combat <- NULL
  }

  if (is.null(exprs_combat)) {
    msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Using limma::removeBatchEffect fallback"))
    cat(msg, "\n")
    cat(msg, "\n", file = LOG_FILE, append = TRUE)
    exprs_combat <- try(
      limma::removeBatchEffect(exprs_all, batch = batch_info),
      silent = TRUE
    )
    if (inherits(exprs_combat, "try-error")) exprs_combat <- NULL
  }

  if (is.null(exprs_combat)) stop("Batch correction failed.")
} else {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Single batch after curation -> batch correction skipped"))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
  exprs_combat <- exprs_all
}

storage.mode(exprs_combat) <- "double"
saveRDS(exprs_all, file.path(OUTPUT_DIR, "exprs_all_curated_prefiltered.rds"))
saveRDS(exprs_combat, file.path(OUTPUT_DIR, "exprs_combat_curated.rds"))

msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Preparing clustering objects"))
cat(msg, "\n")
cat(msg, "\n", file = LOG_FILE, append = TRUE)

dist_mat_full <- dist(t(exprs_combat))
hc_full <- hclust(dist_mat_full, method = "ward.D2")

rv_k <- matrixStats::rowVars(exprs_combat)
ord_k <- order(rv_k, decreasing = TRUE)
n_k_features <- min(KSUPPORT_TOP_VAR, nrow(exprs_combat))
exprs_k <- exprs_combat[ord_k[seq_len(n_k_features)], , drop = FALSE]

pc_fit <- prcomp(t(exprs_k), center = TRUE, scale. = FALSE)
n_pc_use <- min(KSUPPORT_N_PC, ncol(pc_fit$x))
if (n_pc_use < 2L) stop("Too few PCs available for k-support evaluation")
exprs_k_eval <- t(pc_fit$x[, seq_len(n_pc_use), drop = FALSE])

write.csv(
  data.frame(
    metric = c("n_ksupport_genes", "n_ksupport_pcs", "bootstrap_B", "bootstrap_prop"),
    value  = c(n_k_features, n_pc_use, BOOT_B, BOOT_PROP),
    stringsAsFactors = FALSE
  ),
  file.path(TABLE_DIR, "cluster_k_support_runinfo.csv"),
  row.names = FALSE
)

dist_mat_k <- dist(t(exprs_k_eval))
hc_k <- hclust(dist_mat_k, method = "ward.D2")

k_metrics <- list()
for (k in K_GRID) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("k-support start: k=%d", k))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
  cl_k <- cutree(hc_k, k = k)
  sil <- try(cluster::silhouette(cl_k, dist_mat_k), silent = TRUE)
  if (inherits(sil, "try-error")) sil <- NULL
  sil_avg <- if (!is.null(sil)) mean(sil[, "sil_width"], na.rm = TRUE) else NA_real_

  set.seed(123 + k)
  n_sub <- max(5L, floor(ncol(exprs_k_eval) * BOOT_PROP))
  ari_vec <- rep(NA_real_, BOOT_B)

  for (b in seq_len(BOOT_B)) {
    idx_b <- sort(sample.int(ncol(exprs_k_eval), size = n_sub, replace = FALSE))
    sub_expr <- exprs_k_eval[, idx_b, drop = FALSE]
    sub_hc <- try(hclust(dist(t(sub_expr)), method = "ward.D2"), silent = TRUE)
    if (inherits(sub_hc, "try-error")) sub_hc <- NULL

    if (!is.null(sub_hc)) {
      sub_cl <- cutree(sub_hc, k = k)

      ok <- complete.cases(cl_k[idx_b], sub_cl)
      x_ari <- as.factor(cl_k[idx_b][ok])
      y_ari <- as.factor(sub_cl[ok])
      n_ari <- length(x_ari)

      if (n_ari > 1) {
        tab_ari <- table(x_ari, y_ari)
        nij <- as.vector(tab_ari)
        ai <- rowSums(tab_ari)
        bj <- colSums(tab_ari)

        nij2 <- ifelse(nij >= 2, choose(nij, 2), 0)
        ai2 <- ifelse(ai >= 2, choose(ai, 2), 0)
        bj2 <- ifelse(bj >= 2, choose(bj, 2), 0)
        total2 <- ifelse(n_ari >= 2, choose(n_ari, 2), 0)

        expected <- if (total2 == 0) 0 else (sum(ai2) * sum(bj2)) / total2
        max_idx <- 0.5 * (sum(ai2) + sum(bj2))
        denom <- max_idx - expected

        if (denom != 0) {
          ari_vec[b] <- (sum(nij2) - expected) / denom
        }
      }
    }

    if (b %% max(1L, KSUPPORT_LOG_EVERY) == 0L || b == BOOT_B) {
      msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("k-support bootstrap [%s]: %d / %d complete", paste0("k=", k), b, BOOT_B))
      cat(msg, "\n")
      cat(msg, "\n", file = LOG_FILE, append = TRUE)
    }
  }

  stab_df <- data.frame(
    k = k,
    bootstrap_n = BOOT_B,
    bootstrap_prop = BOOT_PROP,
    ari_mean = mean(ari_vec, na.rm = TRUE),
    ari_median = median(ari_vec, na.rm = TRUE),
    ari_sd = stats::sd(ari_vec, na.rm = TRUE),
    ari_n_valid = sum(!is.na(ari_vec)),
    stringsAsFactors = FALSE
  )

  k_metrics[[as.character(k)]] <- data.frame(
    k = k,
    avg_silhouette = sil_avg,
    ari_mean = stab_df$ari_mean,
    ari_median = stab_df$ari_median,
    ari_sd = stab_df$ari_sd,
    ari_n_valid = stab_df$ari_n_valid,
    stringsAsFactors = FALSE
  )
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("k-support done: k=%d | silhouette=%.4f | ARI_mean=%.4f", k, sil_avg, stab_df$ari_mean))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

k_metrics_df <- dplyr::bind_rows(k_metrics)
k_metrics_df$rank_silhouette <- rank(-k_metrics_df$avg_silhouette, ties.method = "min")
k_metrics_df$rank_ari_mean   <- rank(-k_metrics_df$ari_mean, ties.method = "min")
k_metrics_df$rank_sum        <- k_metrics_df$rank_silhouette + k_metrics_df$rank_ari_mean

if (is.na(K_FINAL_OVERRIDE)) {
  k_pick <- k_metrics_df %>%
    dplyr::arrange(rank_sum, dplyr::desc(ari_mean), dplyr::desc(avg_silhouette), k) %>%
    dplyr::slice(1)
  K_FINAL <- k_pick$k[1]
  k_pick_reason <- "auto_by_rank_sum"
} else {
  K_FINAL <- as.integer(K_FINAL_OVERRIDE)
  k_pick_reason <- "manual_override"
}

write.csv(k_metrics_df, file.path(TABLE_DIR, "cluster_k_support_metrics.csv"), row.names = FALSE)
write.csv(
  data.frame(K_FINAL = K_FINAL, selection_mode = k_pick_reason, stringsAsFactors = FALSE),
  file.path(TABLE_DIR, "cluster_final_k_selected.csv"),
  row.names = FALSE
)

p_k1 <- ggplot2::ggplot(k_metrics_df, ggplot2::aes(k, avg_silhouette)) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_continuous(breaks = K_GRID) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Support for a 2-cluster structure: silhouette", x = "Number of clusters (k)", y = "Average silhouette")
try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, "cluster_k_support_silhouette.png"), plot = p_k1, width = 6, height = 4, dpi = PLOT_DPI), silent = TRUE)
try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, "cluster_k_support_silhouette.pdf"), plot = p_k1, width = 6, height = 4, device = "pdf"), silent = TRUE)

p_k2 <- ggplot2::ggplot(k_metrics_df, ggplot2::aes(k, ari_mean)) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_continuous(breaks = K_GRID) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Support for a 2-cluster structure: bootstrap stability", x = "Number of clusters (k)", y = "Mean ARI across bootstrap resamples")
try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, "cluster_k_support_bootstrap_ARI.png"), plot = p_k2, width = 6, height = 4, dpi = PLOT_DPI), silent = TRUE)
try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, "cluster_k_support_bootstrap_ARI.pdf"), plot = p_k2, width = 6, height = 4, device = "pdf"), silent = TRUE)

cluster_assignments <- cutree(hc_full, k = K_FINAL)
cluster_df <- data.frame(
  SampleID = names(cluster_assignments),
  Cluster  = factor(cluster_assignments, levels = seq_len(K_FINAL), labels = paste0("C", seq_len(K_FINAL))),
  stringsAsFactors = FALSE
)
write.csv(cluster_df, file.path(TABLE_DIR, paste0("cluster_assignments_k", K_FINAL, ".csv")), row.names = FALSE)

sample_master2 <- sample_master %>%
  dplyr::filter(sample_decision == "KEEP") %>%
  dplyr::left_join(cluster_df, by = c("sampleid" = "SampleID"))

write.csv(sample_master2, file.path(TABLE_DIR, "sample_master_metadata_curated_with_cluster.csv"), row.names = FALSE)

cluster_sizes <- as.data.frame(table(cluster_df$Cluster), stringsAsFactors = FALSE)
names(cluster_sizes) <- c("Cluster", "N")
write.csv(cluster_sizes, file.path(TABLE_DIR, paste0("cluster_sizes_k", K_FINAL, ".csv")), row.names = FALSE)

if (RUN_SENSITIVITY_EXCLUDE_GSE && length(SENSITIVITY_EXCLUDE_GSE) > 0) {
  sens_drop <- unique(SENSITIVITY_EXCLUDE_GSE)
  sens_keep_idx <- !sub("_.*$", "", colnames(exprs_all)) %in% sens_drop

  if (sum(sens_keep_idx) >= 30) {
    msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Sensitivity clustering: excluding %s from the 9-GEO discovery set", paste(sens_drop, collapse = ", ")))
    cat(msg, "\n")
    cat(msg, "\n", file = LOG_FILE, append = TRUE)

    sens_exprs_all <- exprs_all[, sens_keep_idx, drop = FALSE]
    sens_batch_info <- batch_info[sens_keep_idx]

    if (length(unique(sens_batch_info)) >= 2) {
      sens_exprs <- try(
        sva::ComBat(
          dat = sens_exprs_all,
          batch = sens_batch_info,
          par.prior = TRUE,
          prior.plots = FALSE
        ),
        silent = TRUE
      )
      if (inherits(sens_exprs, "try-error")) {
        msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Sensitivity ComBat failed")
        cat(msg, "\n")
        cat(msg, "\n", file = LOG_FILE, append = TRUE)
        sens_exprs <- limma::removeBatchEffect(sens_exprs_all, batch = sens_batch_info)
      }
    } else {
      sens_exprs <- sens_exprs_all
    }

    storage.mode(sens_exprs) <- "double"

    sens_rv_k <- matrixStats::rowVars(sens_exprs)
    sens_ord_k <- order(sens_rv_k, decreasing = TRUE)
    sens_n_k_features <- min(KSUPPORT_TOP_VAR, nrow(sens_exprs))
    sens_exprs_k <- sens_exprs[sens_ord_k[seq_len(sens_n_k_features)], , drop = FALSE]

    sens_pc_fit <- prcomp(t(sens_exprs_k), center = TRUE, scale. = FALSE)
    sens_n_pc_use <- min(KSUPPORT_N_PC, ncol(sens_pc_fit$x))

    if (sens_n_pc_use >= 2L) {
      sens_exprs_k_eval <- t(sens_pc_fit$x[, seq_len(sens_n_pc_use), drop = FALSE])
      sens_dist_mat_k <- dist(t(sens_exprs_k_eval))
      sens_hc_k <- hclust(sens_dist_mat_k, method = "ward.D2")

      sens_k_metrics <- list()
      for (k in K_GRID) {
        msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Sensitivity k-support start: k=%d", k))
        cat(msg, "\n")
        cat(msg, "\n", file = LOG_FILE, append = TRUE)

        sens_cl_k <- cutree(sens_hc_k, k = k)
        sens_sil <- try(cluster::silhouette(sens_cl_k, sens_dist_mat_k), silent = TRUE)
        if (inherits(sens_sil, "try-error")) sens_sil <- NULL
        sens_sil_avg <- if (!is.null(sens_sil)) mean(sens_sil[, "sil_width"], na.rm = TRUE) else NA_real_

        set.seed(123 + k)
        n_sub <- max(5L, floor(ncol(sens_exprs_k_eval) * BOOT_PROP))
        ari_vec <- rep(NA_real_, BOOT_B)

        for (b in seq_len(BOOT_B)) {
          idx_b <- sort(sample.int(ncol(sens_exprs_k_eval), size = n_sub, replace = FALSE))
          sub_expr <- sens_exprs_k_eval[, idx_b, drop = FALSE]
          sub_hc <- try(hclust(dist(t(sub_expr)), method = "ward.D2"), silent = TRUE)
    if (inherits(sub_hc, "try-error")) sub_hc <- NULL

          if (!is.null(sub_hc)) {
            sub_cl <- cutree(sub_hc, k = k)

            ok <- complete.cases(sens_cl_k[idx_b], sub_cl)
            x_ari <- as.factor(sens_cl_k[idx_b][ok])
            y_ari <- as.factor(sub_cl[ok])
            n_ari <- length(x_ari)

            if (n_ari > 1) {
              tab_ari <- table(x_ari, y_ari)
              nij <- as.vector(tab_ari)
              ai <- rowSums(tab_ari)
              bj <- colSums(tab_ari)

              nij2 <- ifelse(nij >= 2, choose(nij, 2), 0)
              ai2 <- ifelse(ai >= 2, choose(ai, 2), 0)
              bj2 <- ifelse(bj >= 2, choose(bj, 2), 0)
              total2 <- ifelse(n_ari >= 2, choose(n_ari, 2), 0)

              expected <- if (total2 == 0) 0 else (sum(ai2) * sum(bj2)) / total2
              max_idx <- 0.5 * (sum(ai2) + sum(bj2))
              denom <- max_idx - expected

              if (denom != 0) {
                ari_vec[b] <- (sum(nij2) - expected) / denom
              }
            }
          }

          if (b %% max(1L, KSUPPORT_LOG_EVERY) == 0L || b == BOOT_B) {
            msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("k-support bootstrap [%s]: %d / %d complete", paste0("sensitivity_k=", k), b, BOOT_B))
            cat(msg, "\n")
            cat(msg, "\n", file = LOG_FILE, append = TRUE)
          }
        }

        sens_stab_df <- data.frame(
          k = k,
          bootstrap_n = BOOT_B,
          bootstrap_prop = BOOT_PROP,
          ari_mean = mean(ari_vec, na.rm = TRUE),
          ari_median = median(ari_vec, na.rm = TRUE),
          ari_sd = stats::sd(ari_vec, na.rm = TRUE),
          ari_n_valid = sum(!is.na(ari_vec)),
          stringsAsFactors = FALSE
        )

        sens_k_metrics[[as.character(k)]] <- data.frame(
          k = k,
          avg_silhouette = sens_sil_avg,
          ari_mean = sens_stab_df$ari_mean,
          ari_median = sens_stab_df$ari_median,
          ari_sd = sens_stab_df$ari_sd,
          ari_n_valid = sens_stab_df$ari_n_valid,
          stringsAsFactors = FALSE
        )
      }

      sens_k_metrics_df <- dplyr::bind_rows(sens_k_metrics)
      sens_k_metrics_df$rank_silhouette <- rank(-sens_k_metrics_df$avg_silhouette, ties.method = "min")
      sens_k_metrics_df$rank_ari_mean   <- rank(-sens_k_metrics_df$ari_mean, ties.method = "min")
      sens_k_metrics_df$rank_sum        <- sens_k_metrics_df$rank_silhouette + sens_k_metrics_df$rank_ari_mean

      sens_pick <- sens_k_metrics_df %>%
        dplyr::arrange(rank_sum, dplyr::desc(ari_mean), dplyr::desc(avg_silhouette), k) %>%
        dplyr::slice(1)

      sens_K_FINAL <- sens_pick$k[1]
      sens_stub <- paste0("sensitivity_excluding_", paste(sens_drop, collapse = "_"))

      write.csv(
        sens_k_metrics_df,
        file.path(TABLE_DIR, paste0(sens_stub, "_cluster_k_support_metrics.csv")),
        row.names = FALSE
      )

      write.csv(
        data.frame(
          K_FINAL = sens_K_FINAL,
          selection_mode = "auto_by_rank_sum",
          excluded_gse = paste(sens_drop, collapse = ";"),
          stringsAsFactors = FALSE
        ),
        file.path(TABLE_DIR, paste0(sens_stub, "_cluster_final_k_selected.csv")),
        row.names = FALSE
      )

      sens_hc_full <- hclust(dist(t(sens_exprs)), method = "ward.D2")
      sens_cluster_assignments <- cutree(sens_hc_full, k = sens_K_FINAL)

      sens_cluster_df <- data.frame(
        SampleID = names(sens_cluster_assignments),
        Cluster  = factor(
          sens_cluster_assignments,
          levels = seq_len(sens_K_FINAL),
          labels = paste0("C", seq_len(sens_K_FINAL))
        ),
        stringsAsFactors = FALSE
      )

      write.csv(
        sens_cluster_df,
        file.path(TABLE_DIR, paste0(sens_stub, "_cluster_assignments_k", sens_K_FINAL, ".csv")),
        row.names = FALSE
      )

      sens_cluster_sizes <- as.data.frame(table(sens_cluster_df$Cluster), stringsAsFactors = FALSE)
      names(sens_cluster_sizes) <- c("Cluster", "N")
      write.csv(
        sens_cluster_sizes,
        file.path(TABLE_DIR, paste0(sens_stub, "_cluster_sizes_k", sens_K_FINAL, ".csv")),
        row.names = FALSE
      )
    }
  }
}

pal_cluster <- if (has_rcolorbrewer) {
  base_pal <- RColorBrewer::brewer.pal(max(3, min(9, K_FINAL)), "Set1")
  base_pal[seq_len(K_FINAL)]
} else {
  grDevices::hcl.colors(K_FINAL, "Dark 3")
}

try({
  png(file.path(PLOT_PNG_DIR, paste0("dendrogram_k", K_FINAL, ".png")), width = 7, height = 6, units = "in", res = PLOT_DPI)
  plot(hc_full, labels = FALSE, main = "Hierarchical clustering of the curated 9-GEO discovery cohort")
  rect.hclust(hc_full, k = K_FINAL, border = pal_cluster)
  dev.off()
}, silent = TRUE)

try({
  pdf(file.path(PLOT_PDF_DIR, paste0("dendrogram_k", K_FINAL, ".pdf")), width = 7, height = 6)
  plot(hc_full, labels = FALSE, main = "Hierarchical clustering of the curated 9-GEO discovery cohort")
  rect.hclust(hc_full, k = K_FINAL, border = pal_cluster)
  dev.off()
}, silent = TRUE)

pca_res <- prcomp(t(exprs_combat), scale. = TRUE)
pca_df <- data.frame(
  SampleID = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  stringsAsFactors = FALSE
) %>% dplyr::left_join(cluster_df, by = "SampleID")

p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = Cluster)) +
  ggplot2::geom_point(size = 2, alpha = 0.8) +
  ggplot2::scale_color_manual(values = pal_cluster) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "PCA of the curated 9-GEO discovery cohort", color = "Cluster")
try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, "PCA_by_cluster.png"), plot = p_pca, width = 7, height = 5, dpi = PLOT_DPI), silent = TRUE)
try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, "PCA_by_cluster.pdf"), plot = p_pca, width = 7, height = 5, device = "pdf"), silent = TRUE)

if (RUN_UMAP && has_umap) {
  set.seed(42)
  um <- try(umap::umap(t(exprs_combat), config = umap::umap.defaults), silent = TRUE)
  if (inherits(um, "try-error")) um <- NULL
} else {
  um <- NULL
}
if (!is.null(um)) {
  umap_df <- data.frame(
    SampleID = colnames(exprs_combat),
    UMAP1 = um$layout[, 1],
    UMAP2 = um$layout[, 2],
    stringsAsFactors = FALSE
  ) %>% dplyr::left_join(cluster_df, by = "SampleID")

  p_umap <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = Cluster)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = pal_cluster) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "UMAP of the curated 9-GEO discovery cohort", color = "Cluster")
  try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, "UMAP_by_cluster.png"), plot = p_umap, width = 7, height = 5, dpi = PLOT_DPI), silent = TRUE)
  try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, "UMAP_by_cluster.pdf"), plot = p_umap, width = 7, height = 5, device = "pdf"), silent = TRUE)
}

if (RUN_TSNE && has_rtsne) {
  tsne_in <- t(exprs_combat)
  if (any(duplicated(as.data.frame(tsne_in)))) {
    set.seed(43)
    tsne_in <- tsne_in + matrix(rnorm(length(tsne_in), sd = 1e-6), nrow = nrow(tsne_in), ncol = ncol(tsne_in))
  }
  perplexity <- min(30, floor((nrow(tsne_in) - 1) / 3))
  set.seed(44)
  ts <- try(Rtsne::Rtsne(tsne_in, perplexity = perplexity, check_duplicates = FALSE, verbose = FALSE), silent = TRUE)
  if (inherits(ts, "try-error")) ts <- NULL
} else {
  ts <- NULL
}
if (!is.null(ts)) {
  tsne_df <- data.frame(
    SampleID = rownames(tsne_in),
    tSNE1 = ts$Y[, 1],
    tSNE2 = ts$Y[, 2],
    stringsAsFactors = FALSE
  ) %>% dplyr::left_join(cluster_df, by = "SampleID")

  p_tsne <- ggplot2::ggplot(tsne_df, ggplot2::aes(tSNE1, tSNE2, color = Cluster)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = pal_cluster) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "t-SNE of the curated 9-GEO discovery cohort", color = "Cluster")
  try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, "tSNE_by_cluster.png"), plot = p_tsne, width = 7, height = 5, dpi = PLOT_DPI), silent = TRUE)
  try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, "tSNE_by_cluster.pdf"), plot = p_tsne, width = 7, height = 5, device = "pdf"), silent = TRUE)
}

comp_input <- sample_master2

comp_input$gse_id_comp <- as.character(comp_input$gse_id)
comp_input$gse_id_comp <- trimws(comp_input$gse_id_comp)
comp_input$gse_id_comp[comp_input$gse_id_comp == ""] <- NA_character_
comp_input$gse_id_comp <- gsub("\\s+", " ", comp_input$gse_id_comp)
comp_input$diagnosis_comp <- ifelse(
  !is.na(comp_input$diagnosis_group) & comp_input$diagnosis_group != "",
  comp_input$diagnosis_group,
  comp_input$diagnosis_h
)
comp_input$tissue_comp <- ifelse(
  !is.na(comp_input$tissue_group) & comp_input$tissue_group != "",
  comp_input$tissue_group,
  comp_input$tissue_h
)

for (nm in c("diagnosis_comp", "tissue_comp")) {
  x <- as.character(comp_input[[nm]])
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x <- gsub("\\s+", " ", x)
  tab <- sort(table(x), decreasing = TRUE)
  keep <- names(tab)[tab >= COMP_MIN_N]
  if (length(keep) > COMP_TOP_N) keep <- keep[seq_len(COMP_TOP_N)]
  comp_input[[nm]] <- ifelse(is.na(x), NA_character_, ifelse(x %in% keep, x, "Other"))
}

composition_stats <- list()

v <- "gse_id_comp"
out_nm <- "gse_id"
tmp <- comp_input[, c("sampleid", "Cluster", v), drop = FALSE]
names(tmp)[3] <- "Category"
tmp <- tmp[!is.na(tmp$Cluster) & !is.na(tmp$Category) & tmp$Category != "", , drop = FALSE]
if (nrow(tmp) > 0 && length(unique(tmp$Category)) >= 2) {
  tab <- table(tmp$Cluster, tmp$Category)
  write.csv(as.data.frame.matrix(tab), file.path(TABLE_DIR, paste0("composition_count_", out_nm, "_by_cluster.csv")))
  prop_tab <- prop.table(tab, margin = 1)
  write.csv(as.data.frame.matrix(round(prop_tab, 4)), file.path(TABLE_DIR, paste0("composition_prop_", out_nm, "_by_cluster.csv")))
  if (length(dim(tab)) != 2 || nrow(tab) < 2 || ncol(tab) < 2) {
    test_res <- data.frame(test = "not_applicable", statistic = NA_real_, df = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
  } else {
    x0 <- try(suppressWarnings(chisq.test(tab, correct = FALSE)), silent = TRUE)
    if (inherits(x0, "try-error")) x0 <- NULL
    if (is.null(x0)) {
      x1 <- try(chisq.test(tab, simulate.p.value = TRUE, B = COMP_SIM_B), silent = TRUE)
      if (inherits(x1, "try-error")) x1 <- NULL
      if (is.null(x1)) {
        test_res <- data.frame(test = "failed", statistic = NA_real_, df = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
      } else {
        test_res <- data.frame(test = "chisq_sim", statistic = unname(x1$statistic), df = if (is.null(x1$parameter)) NA_real_ else unname(x1$parameter), p_value = x1$p.value, stringsAsFactors = FALSE)
      }
    } else {
      use_sim <- any(x0$expected < 5) || any(tab == 0) || prod(dim(tab)) > 20L
      if (use_sim) {
        x1 <- try(chisq.test(tab, simulate.p.value = TRUE, B = COMP_SIM_B), silent = TRUE)
        if (inherits(x1, "try-error")) x1 <- NULL
      } else {
        x1 <- NULL
      }
      if (!is.null(x1)) {
        test_res <- data.frame(test = "chisq_sim", statistic = unname(x1$statistic), df = if (is.null(x1$parameter)) NA_real_ else unname(x1$parameter), p_value = x1$p.value, stringsAsFactors = FALSE)
      } else {
        test_res <- data.frame(test = "chisq", statistic = unname(x0$statistic), df = if (is.null(x0$parameter)) NA_real_ else unname(x0$parameter), p_value = x0$p.value, stringsAsFactors = FALSE)
      }
    }
  }
  test_res$variable <- out_nm
  test_res$n_samples_used <- nrow(tmp)
  test_res$n_categories <- ncol(tab)
  composition_stats[[out_nm]] <- test_res[, c("variable", "test", "statistic", "df", "p_value", "n_samples_used", "n_categories")]
  plot_df <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(plot_df) <- c("Cluster", "Category", "N")
  plot_label <- "Dataset (GEO Series)"
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(Cluster, N, fill = Category)) + ggplot2::geom_col(position = "fill") + ggplot2::theme_bw() + ggplot2::labs(title = paste("Distribution of", tolower(plot_label), "across clusters"), x = "Cluster", y = "Proportion within cluster", fill = plot_label)
  try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, paste0("composition_", out_nm, "_stacked.png")), plot = p, width = 8, height = 5, dpi = PLOT_DPI), silent = TRUE)
  try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, paste0("composition_", out_nm, "_stacked.pdf")), plot = p, width = 8, height = 5, device = "pdf"), silent = TRUE)
}

v <- "diagnosis_comp"
out_nm <- "diagnosis"
tmp <- comp_input[, c("sampleid", "Cluster", v), drop = FALSE]
names(tmp)[3] <- "Category"
tmp <- tmp[!is.na(tmp$Cluster) & !is.na(tmp$Category) & tmp$Category != "", , drop = FALSE]
if (nrow(tmp) > 0 && length(unique(tmp$Category)) >= 2) {
  tab <- table(tmp$Cluster, tmp$Category)
  write.csv(as.data.frame.matrix(tab), file.path(TABLE_DIR, paste0("composition_count_", out_nm, "_by_cluster.csv")))
  prop_tab <- prop.table(tab, margin = 1)
  write.csv(as.data.frame.matrix(round(prop_tab, 4)), file.path(TABLE_DIR, paste0("composition_prop_", out_nm, "_by_cluster.csv")))
  if (length(dim(tab)) != 2 || nrow(tab) < 2 || ncol(tab) < 2) {
    test_res <- data.frame(test = "not_applicable", statistic = NA_real_, df = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
  } else {
    x0 <- try(suppressWarnings(chisq.test(tab, correct = FALSE)), silent = TRUE)
    if (inherits(x0, "try-error")) x0 <- NULL
    if (is.null(x0)) {
      x1 <- try(chisq.test(tab, simulate.p.value = TRUE, B = COMP_SIM_B), silent = TRUE)
      if (inherits(x1, "try-error")) x1 <- NULL
      if (is.null(x1)) {
        test_res <- data.frame(test = "failed", statistic = NA_real_, df = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
      } else {
        test_res <- data.frame(test = "chisq_sim", statistic = unname(x1$statistic), df = if (is.null(x1$parameter)) NA_real_ else unname(x1$parameter), p_value = x1$p.value, stringsAsFactors = FALSE)
      }
    } else {
      use_sim <- any(x0$expected < 5) || any(tab == 0) || prod(dim(tab)) > 20L
      if (use_sim) {
        x1 <- try(chisq.test(tab, simulate.p.value = TRUE, B = COMP_SIM_B), silent = TRUE)
        if (inherits(x1, "try-error")) x1 <- NULL
      } else {
        x1 <- NULL
      }
      if (!is.null(x1)) {
        test_res <- data.frame(test = "chisq_sim", statistic = unname(x1$statistic), df = if (is.null(x1$parameter)) NA_real_ else unname(x1$parameter), p_value = x1$p.value, stringsAsFactors = FALSE)
      } else {
        test_res <- data.frame(test = "chisq", statistic = unname(x0$statistic), df = if (is.null(x0$parameter)) NA_real_ else unname(x0$parameter), p_value = x0$p.value, stringsAsFactors = FALSE)
      }
    }
  }
  test_res$variable <- out_nm
  test_res$n_samples_used <- nrow(tmp)
  test_res$n_categories <- ncol(tab)
  composition_stats[[out_nm]] <- test_res[, c("variable", "test", "statistic", "df", "p_value", "n_samples_used", "n_categories")]
  plot_df <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(plot_df) <- c("Cluster", "Category", "N")
  plot_label <- "Diagnosis group"
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(Cluster, N, fill = Category)) + ggplot2::geom_col(position = "fill") + ggplot2::theme_bw() + ggplot2::labs(title = paste("Distribution of", tolower(plot_label), "across clusters"), x = "Cluster", y = "Proportion within cluster", fill = plot_label)
  try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, paste0("composition_", out_nm, "_stacked.png")), plot = p, width = 8, height = 5, dpi = PLOT_DPI), silent = TRUE)
  try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, paste0("composition_", out_nm, "_stacked.pdf")), plot = p, width = 8, height = 5, device = "pdf"), silent = TRUE)
}

v <- "tissue_comp"
out_nm <- "tissue"
tmp <- comp_input[, c("sampleid", "Cluster", v), drop = FALSE]
names(tmp)[3] <- "Category"
tmp <- tmp[!is.na(tmp$Cluster) & !is.na(tmp$Category) & tmp$Category != "", , drop = FALSE]
if (nrow(tmp) > 0 && length(unique(tmp$Category)) >= 2) {
  tab <- table(tmp$Cluster, tmp$Category)
  write.csv(as.data.frame.matrix(tab), file.path(TABLE_DIR, paste0("composition_count_", out_nm, "_by_cluster.csv")))
  prop_tab <- prop.table(tab, margin = 1)
  write.csv(as.data.frame.matrix(round(prop_tab, 4)), file.path(TABLE_DIR, paste0("composition_prop_", out_nm, "_by_cluster.csv")))
  if (length(dim(tab)) != 2 || nrow(tab) < 2 || ncol(tab) < 2) {
    test_res <- data.frame(test = "not_applicable", statistic = NA_real_, df = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
  } else {
    x0 <- try(suppressWarnings(chisq.test(tab, correct = FALSE)), silent = TRUE)
    if (inherits(x0, "try-error")) x0 <- NULL
    if (is.null(x0)) {
      x1 <- try(chisq.test(tab, simulate.p.value = TRUE, B = COMP_SIM_B), silent = TRUE)
      if (inherits(x1, "try-error")) x1 <- NULL
      if (is.null(x1)) {
        test_res <- data.frame(test = "failed", statistic = NA_real_, df = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
      } else {
        test_res <- data.frame(test = "chisq_sim", statistic = unname(x1$statistic), df = if (is.null(x1$parameter)) NA_real_ else unname(x1$parameter), p_value = x1$p.value, stringsAsFactors = FALSE)
      }
    } else {
      use_sim <- any(x0$expected < 5) || any(tab == 0) || prod(dim(tab)) > 20L
      if (use_sim) {
        x1 <- try(chisq.test(tab, simulate.p.value = TRUE, B = COMP_SIM_B), silent = TRUE)
        if (inherits(x1, "try-error")) x1 <- NULL
      } else {
        x1 <- NULL
      }
      if (!is.null(x1)) {
        test_res <- data.frame(test = "chisq_sim", statistic = unname(x1$statistic), df = if (is.null(x1$parameter)) NA_real_ else unname(x1$parameter), p_value = x1$p.value, stringsAsFactors = FALSE)
      } else {
        test_res <- data.frame(test = "chisq", statistic = unname(x0$statistic), df = if (is.null(x0$parameter)) NA_real_ else unname(x0$parameter), p_value = x0$p.value, stringsAsFactors = FALSE)
      }
    }
  }
  test_res$variable <- out_nm
  test_res$n_samples_used <- nrow(tmp)
  test_res$n_categories <- ncol(tab)
  composition_stats[[out_nm]] <- test_res[, c("variable", "test", "statistic", "df", "p_value", "n_samples_used", "n_categories")]
  plot_df <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(plot_df) <- c("Cluster", "Category", "N")
  plot_label <- "Tissue group"
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(Cluster, N, fill = Category)) + ggplot2::geom_col(position = "fill") + ggplot2::theme_bw() + ggplot2::labs(title = paste("Distribution of", tolower(plot_label), "across clusters"), x = "Cluster", y = "Proportion within cluster", fill = plot_label)
  try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, paste0("composition_", out_nm, "_stacked.png")), plot = p, width = 8, height = 5, dpi = PLOT_DPI), silent = TRUE)
  try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, paste0("composition_", out_nm, "_stacked.pdf")), plot = p, width = 8, height = 5, device = "pdf"), silent = TRUE)
}

if (length(composition_stats) > 0) {
  composition_stats_df <- dplyr::bind_rows(composition_stats)
  write.csv(
    composition_stats_df,
    file.path(TABLE_DIR, "cluster_composition_association_tests.csv"),
    row.names = FALSE
  )
}

dom_df <- comp_input %>%
  dplyr::filter(!is.na(Cluster), !is.na(gse_id_comp), gse_id_comp != "") %>%
  dplyr::count(Cluster, gse_id_comp, name = "n") %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::arrange(Cluster, dplyr::desc(prop)) %>%
  dplyr::ungroup()

names(dom_df)[names(dom_df) == "gse_id_comp"] <- "gse_id"

write.csv(dom_df, file.path(TABLE_DIR, "cluster_dataset_distribution.csv"), row.names = FALSE)
write.csv(
  dom_df %>% dplyr::group_by(Cluster) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup(),
  file.path(TABLE_DIR, "cluster_dataset_top_driver.csv"),
  row.names = FALSE
)

groups <- factor(
  cluster_df$Cluster[match(colnames(exprs_combat), cluster_df$SampleID)],
  levels = paste0("C", seq_len(K_FINAL))
)
stopifnot(!any(is.na(groups)))

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

fit <- limma::lmFit(exprs_combat, design)

contrast_strings <- character(length(levels(groups)))
for (i_contrast in seq_along(levels(groups))) {
  lhs <- levels(groups)[i_contrast]
  rhs <- setdiff(levels(groups), lhs)
  contrast_strings[i_contrast] <- paste0(lhs, " - (", paste(rhs, collapse = " + "), ")/", length(rhs))
}
names(contrast_strings) <- NULL
names(contrast_strings) <- paste0(levels(groups), "_vsRest")

contr <- do.call(limma::makeContrasts, c(as.list(contrast_strings), list(levels = design)))
fit2 <- limma::eBayes(limma::contrasts.fit(fit, contr))

deg_tables <- list()
sig_lists  <- list()
rank_lists <- list()

for (coef_name in colnames(contr)) {
  tt <- limma::topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  tt$GeneSymbol <- rownames(tt)
  tt$Contrast   <- coef_name
  tt$direction  <- ifelse(tt$logFC > 0, "Up", "Down")
  deg_tables[[coef_name]] <- tt
  rank_lists[[coef_name]] <- setNames(tt$t, tt$GeneSymbol)

  write.csv(tt, file.path(TABLE_DIR, paste0("DEG_", coef_name, ".csv")), row.names = FALSE)

  sig_lists[[coef_name]] <- tt$GeneSymbol[tt$adj.P.Val < DEG_FDR_CUTOFF & abs(tt$logFC) > DEG_LOGFC_CUTOFF]

  eps <- .Machine$double.xmin
  tt$neglog10P <- -log10(pmax(tt$P.Value, eps))
  tt$signif <- ifelse(tt$adj.P.Val < DEG_FDR_CUTOFF & abs(tt$logFC) > DEG_LOGFC_CUTOFF, "sig", "ns")
  lab_idx <- unique(c(
    order(tt$adj.P.Val)[1:min(20, nrow(tt))],
    order(abs(tt$logFC), decreasing = TRUE)[1:min(20, nrow(tt))]
  ))
  tt$label <- ""
  tt$label[lab_idx] <- tt$GeneSymbol[lab_idx]

  p <- ggplot2::ggplot(tt, ggplot2::aes(logFC, neglog10P, color = signif)) +
    ggplot2::geom_point(size = 1, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = c(-DEG_LOGFC_CUTOFF, DEG_LOGFC_CUTOFF), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(DEG_FDR_CUTOFF), linetype = "dashed") +
    ggplot2::scale_color_manual(values = c(ns = "grey70", sig = "#00BFC4")) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste0(coef_name, " volcano"), x = "logFC", y = expression(-log[10](p)))

  if (has_ggrepel) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(tt, label != ""),
      ggplot2::aes(label = label),
      max.overlaps = Inf, size = 3, show.legend = FALSE
    )
  }
  try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, paste0("Volcano_", coef_name, ".png")), plot = p, width = 7, height = 5, dpi = PLOT_DPI), silent = TRUE)
  try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, paste0("Volcano_", coef_name, ".pdf")), plot = p, width = 7, height = 5, device = "pdf"), silent = TRUE)
}

sig_lists_up <- vector("list", length(deg_tables))
names(sig_lists_up) <- names(deg_tables)
for (nm_sig in names(deg_tables)) {
  tt <- deg_tables[[nm_sig]]
  sig_lists_up[[nm_sig]] <- unique(tt$GeneSymbol[tt$adj.P.Val < DEG_FDR_CUTOFF & tt$logFC > DEG_LOGFC_CUTOFF])
}
sig_lists_down <- vector("list", length(deg_tables))
names(sig_lists_down) <- names(deg_tables)
for (nm_sig in names(deg_tables)) {
  tt <- deg_tables[[nm_sig]]
  sig_lists_down[[nm_sig]] <- unique(tt$GeneSymbol[tt$adj.P.Val < DEG_FDR_CUTOFF & tt$logFC < -DEG_LOGFC_CUTOFF])
}

deg_summary <- data.frame(
  Contrast = names(deg_tables),
  n_sig_any = vapply(sig_lists, length, integer(1)),
  n_up = vapply(sig_lists_up, length, integer(1)),
  n_down = vapply(sig_lists_down, length, integer(1)),
  stringsAsFactors = FALSE
)
write.csv(deg_summary, file.path(TABLE_DIR, "DEG_summary_counts.csv"), row.names = FALSE)

marker_export <- list()
for (nm in names(deg_tables)) {
  tt <- deg_tables[[nm]]

  up_tt <- tt[tt$adj.P.Val < DEG_FDR_CUTOFF & tt$logFC > DEG_LOGFC_CUTOFF,
              c("GeneSymbol", "logFC", "adj.P.Val", "P.Value", "t"),
              drop = FALSE]
  down_tt <- tt[tt$adj.P.Val < DEG_FDR_CUTOFF & tt$logFC < -DEG_LOGFC_CUTOFF,
                c("GeneSymbol", "logFC", "adj.P.Val", "P.Value", "t"),
                drop = FALSE]

  up_tt$Contrast <- nm
  up_tt$Direction <- "UP"
  down_tt$Contrast <- nm
  down_tt$Direction <- "DOWN"

  marker_export[[paste0(nm, "_UP")]] <- head(up_tt, TOPN_MARKER_EXPORT)
  marker_export[[paste0(nm, "_DOWN")]] <- head(down_tt, TOPN_MARKER_EXPORT)

  write.csv(
    head(up_tt, TOPN_MARKER_EXPORT),
    file.path(TABLE_DIR, paste0("DEG_top_up_", nm, ".csv")),
    row.names = FALSE
  )
  write.csv(
    head(down_tt, TOPN_MARKER_EXPORT),
    file.path(TABLE_DIR, paste0("DEG_top_down_", nm, ".csv")),
    row.names = FALSE
  )
}

if (length(marker_export) > 0) {
  write.csv(
    dplyr::bind_rows(marker_export),
    file.path(TABLE_DIR, "DEG_top_markers_directional_combined.csv"),
    row.names = FALSE
  )
}

all_gene_list <- vector("list", length(deg_tables))
names(all_gene_list) <- names(deg_tables)
for (nm_all in names(deg_tables)) all_gene_list[[nm_all]] <- deg_tables[[nm_all]]$GeneSymbol
all_genes <- unique(unlist(all_gene_list))
sig_mat <- data.frame(GeneSymbol = all_genes, stringsAsFactors = FALSE)

for (nm in names(deg_tables)) {
  tt <- deg_tables[[nm]][, c("GeneSymbol", "logFC", "adj.P.Val", "direction")]
  names(tt)[-1] <- paste0(names(tt)[-1], "_", nm)
  sig_mat <- merge(sig_mat, tt, by = "GeneSymbol", all.x = TRUE, all.y = FALSE)
  sig_mat[[paste0("sig_", nm)]] <- with(sig_mat, !is.na(get(paste0("adj.P.Val_", nm))) & get(paste0("adj.P.Val_", nm)) < DEG_FDR_CUTOFF & abs(get(paste0("logFC_", nm))) > DEG_LOGFC_CUTOFF)
}

sig_cols <- grep("^sig_", names(sig_mat), value = TRUE)
sig_mat$n_sig_contrasts <- rowSums(sig_mat[, sig_cols, drop = FALSE], na.rm = TRUE)
write.csv(sig_mat, file.path(TABLE_DIR, "DEG_signature_matrix.csv"), row.names = FALSE)

pairwise_overlap <- list()
if (length(sig_lists) >= 2) {
  cmb <- utils::combn(names(sig_lists), 2)
  for (i in seq_len(ncol(cmb))) {
    a <- cmb[1, i]; b <- cmb[2, i]
    sa <- unique(sig_lists[[a]])
    sb <- unique(sig_lists[[b]])
    inter <- length(intersect(sa, sb))
    uni   <- length(unique(c(sa, sb)))
    pairwise_overlap[[i]] <- data.frame(
      ContrastA = a,
      ContrastB = b,
      nA = length(sa),
      nB = length(sb),
      n_intersection = inter,
      jaccard = ifelse(uni == 0, NA_real_, inter / uni),
      stringsAsFactors = FALSE
    )
  }
}
if (length(pairwise_overlap) > 0) {
  write.csv(dplyr::bind_rows(pairwise_overlap), file.path(TABLE_DIR, "DEG_pairwise_overlap_summary.csv"), row.names = FALSE)
}

if (K_FINAL == 3) {
  all3 <- names(sig_lists)
  shared3 <- Reduce(intersect, sig_lists[all3])
  write.csv(data.frame(GeneSymbol = shared3), file.path(TABLE_DIR, "DEG_threeway_shared_sig_genes.csv"), row.names = FALSE)

  overlap_counts <- data.frame(
    set = c("C1_only","C2_only","C3_only","C1_C2","C1_C3","C2_C3","C1_C2_C3"),
    n = c(
      sum(sig_mat[[paste0("sig_", all3[1])]] & !sig_mat[[paste0("sig_", all3[2])]] & !sig_mat[[paste0("sig_", all3[3])]], na.rm = TRUE),
      sum(!sig_mat[[paste0("sig_", all3[1])]] & sig_mat[[paste0("sig_", all3[2])]] & !sig_mat[[paste0("sig_", all3[3])]], na.rm = TRUE),
      sum(!sig_mat[[paste0("sig_", all3[1])]] & !sig_mat[[paste0("sig_", all3[2])]] & sig_mat[[paste0("sig_", all3[3])]], na.rm = TRUE),
      sum(sig_mat[[paste0("sig_", all3[1])]] & sig_mat[[paste0("sig_", all3[2])]] & !sig_mat[[paste0("sig_", all3[3])]], na.rm = TRUE),
      sum(sig_mat[[paste0("sig_", all3[1])]] & !sig_mat[[paste0("sig_", all3[2])]] & sig_mat[[paste0("sig_", all3[3])]], na.rm = TRUE),
      sum(!sig_mat[[paste0("sig_", all3[1])]] & sig_mat[[paste0("sig_", all3[2])]] & sig_mat[[paste0("sig_", all3[3])]], na.rm = TRUE),
      sum(sig_mat[[paste0("sig_", all3[1])]] & sig_mat[[paste0("sig_", all3[2])]] & sig_mat[[paste0("sig_", all3[3])]], na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(overlap_counts, file.path(TABLE_DIR, "DEG_overlap_counts.csv"), row.names = FALSE)
}

signature_gene_list <- vector("list", length(deg_tables))
names(signature_gene_list) <- names(deg_tables)
for (nm_sig2 in names(deg_tables)) {
  tt <- deg_tables[[nm_sig2]]
  tt <- tt[tt$adj.P.Val < DEG_FDR_CUTOFF & tt$logFC > DEG_LOGFC_CUTOFF, , drop = FALSE]
  signature_gene_list[[nm_sig2]] <- head(tt$GeneSymbol, TOPN_DEG_FOR_SIG)
}
signature_genes <- unique(unlist(signature_gene_list))

signature_genes <- intersect(signature_genes, rownames(exprs_combat))

if (length(signature_genes) < 10) {
  signature_gene_list2 <- vector("list", length(deg_tables))
  names(signature_gene_list2) <- names(deg_tables)
  for (nm_sig3 in names(deg_tables)) {
    tt <- deg_tables[[nm_sig3]]
    tt <- tt[tt$logFC > 0, , drop = FALSE]
    signature_gene_list2[[nm_sig3]] <- head(tt$GeneSymbol, TOPN_DEG_FOR_SIG)
  }
  signature_genes <- unique(unlist(signature_gene_list2))
  signature_genes <- intersect(signature_genes, rownames(exprs_combat))
}

write.csv(
  data.frame(GeneSymbol = signature_genes),
  file.path(TABLE_DIR, "downstream_signature_genes_used.csv"),
  row.names = FALSE
)

msig_h <- NULL
msig_r <- NULL

if (has_msigdbr) {
  f_msig <- try(get("msigdbr", envir = asNamespace("msigdbr")), silent = TRUE)
  if (inherits(f_msig, "try-error")) f_msig <- NULL

  if (!is.null(f_msig)) {
    fn_msig <- names(formals(f_msig))

    args_h <- list(species = "Homo sapiens")
    if ("collection" %in% fn_msig) args_h$collection <- "H"
    if ("category" %in% fn_msig && !("collection" %in% fn_msig)) args_h$category <- "H"
    out_h <- try(do.call(f_msig, args_h), silent = TRUE)
    if (inherits(out_h, "try-error")) out_h <- NULL
    if (!is.null(out_h) && nrow(out_h) > 0) msig_h <- split(out_h$gene_symbol, out_h$gs_name)

    args_r <- list(species = "Homo sapiens")
    if ("collection" %in% fn_msig) args_r$collection <- "C2"
    if ("subcollection" %in% fn_msig) args_r$subcollection <- "CP:REACTOME"
    if ("category" %in% fn_msig && !("collection" %in% fn_msig)) args_r$category <- "C2"
    if ("subcategory" %in% fn_msig && !("subcollection" %in% fn_msig)) args_r$subcategory <- "CP:REACTOME"
    out_r <- try(do.call(f_msig, args_r), silent = TRUE)
    if (inherits(out_r, "try-error")) out_r <- NULL
    if (!is.null(out_r) && nrow(out_r) > 0) msig_r <- split(out_r$gene_symbol, out_r$gs_name)
  }
}

gsea_h <- vector("list", length(rank_lists))
names(gsea_h) <- names(rank_lists)
if (!is.null(msig_h) && has_fgsea) {
  for (nm in names(rank_lists)) {
    stats_vec <- rank_lists[[nm]]
    pathways <- msig_h
    out <- try(
      fgsea::fgseaMultilevel(pathways = pathways, stats = stats_vec, minSize = 10, maxSize = 500),
      silent = TRUE
    )
    if (inherits(out, "try-error")) out <- NULL
    if (is.null(out)) {
      out <- try(
        fgsea::fgsea(pathways = pathways, stats = stats_vec, minSize = 10, maxSize = 500, nperm = 20000),
        silent = TRUE
      )
      if (inherits(out, "try-error")) out <- NULL
    }
    if (!is.null(out) && nrow(out) > 0) {
      out <- dplyr::arrange(out, padj, dplyr::desc(abs(NES)))
      if ("leadingEdge" %in% names(out)) {
        leading_edge_vec <- character(length(out$leadingEdge))
        for (i_le in seq_along(out$leadingEdge)) leading_edge_vec[i_le] <- paste(head(out$leadingEdge[[i_le]], 10), collapse = ",")
        out$leadingEdge <- leading_edge_vec
      }
      gsea_h[[nm]] <- out
    }
  }
}

gsea_r <- vector("list", length(rank_lists))
names(gsea_r) <- names(rank_lists)
if (!is.null(msig_r) && has_fgsea) {
  for (nm in names(rank_lists)) {
    stats_vec <- rank_lists[[nm]]
    pathways <- msig_r
    out <- try(
      fgsea::fgseaMultilevel(pathways = pathways, stats = stats_vec, minSize = 10, maxSize = 500),
      silent = TRUE
    )
    if (inherits(out, "try-error")) out <- NULL
    if (is.null(out)) {
      out <- try(
        fgsea::fgsea(pathways = pathways, stats = stats_vec, minSize = 10, maxSize = 500, nperm = 20000),
        silent = TRUE
      )
      if (inherits(out, "try-error")) out <- NULL
    }
    if (!is.null(out) && nrow(out) > 0) {
      out <- dplyr::arrange(out, padj, dplyr::desc(abs(NES)))
      if ("leadingEdge" %in% names(out)) {
        leading_edge_vec <- character(length(out$leadingEdge))
        for (i_le in seq_along(out$leadingEdge)) leading_edge_vec[i_le] <- paste(head(out$leadingEdge[[i_le]], 10), collapse = ",")
        out$leadingEdge <- leading_edge_vec
      }
      gsea_r[[nm]] <- out
    }
  }
}

for (nm in names(gsea_h)) {
  if (!is.null(gsea_h[[nm]]) && nrow(gsea_h[[nm]]) > 0) {
    write.csv(gsea_h[[nm]], file.path(TABLE_DIR, paste0("GSEA_Hallmark_", nm, ".csv")), row.names = FALSE)
    top10 <- dplyr::slice_min(gsea_h[[nm]], padj, n = 10)
    p <- ggplot2::ggplot(top10, ggplot2::aes(reorder(pathway, NES), NES)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::labs(title = paste0(nm, ": Hallmark GSEA (top terms by FDR)"), x = NULL, y = "Normalized enrichment score")
    try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, paste0("GSEA_Hallmark_", nm, "_top10.png")), plot = p, width = 8, height = 5, dpi = PLOT_DPI), silent = TRUE)
    try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, paste0("GSEA_Hallmark_", nm, "_top10.pdf")), plot = p, width = 8, height = 5, device = "pdf"), silent = TRUE)
  }
}

for (nm in names(gsea_r)) {
  if (!is.null(gsea_r[[nm]]) && nrow(gsea_r[[nm]]) > 0) {
    write.csv(gsea_r[[nm]], file.path(TABLE_DIR, paste0("GSEA_Reactome_", nm, ".csv")), row.names = FALSE)
  }
}

hallmark_focus <- c(
  "HALLMARK_COMPLEMENT",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_ALLOGRAFT_REJECTION"
)

gsea_focus_summary <- list()
for (nm in names(gsea_h)) {
  if (!is.null(gsea_h[[nm]]) && nrow(gsea_h[[nm]]) > 0) {
    tmp <- gsea_h[[nm]]
    tmp <- tmp[tmp$pathway %in% hallmark_focus, c("pathway", "NES", "padj"), drop = FALSE]
    if (nrow(tmp) > 0) {
      tmp$Contrast <- nm
      tmp$direction <- ifelse(tmp$NES > 0, "higher_in_contrast", "lower_in_contrast")
      gsea_focus_summary[[nm]] <- tmp
    }
  }
}
if (length(gsea_focus_summary) > 0) {
  write.csv(
    dplyr::bind_rows(gsea_focus_summary),
    file.path(TABLE_DIR, "GSEA_Hallmark_focus_summary.csv"),
    row.names = FALSE
  )
}

if (has_clusterprofiler) {
  ora_sets <- list(UP = sig_lists_up)
  if (isTRUE(RUN_DOWN_ORA_SUPP)) ora_sets$DOWN <- sig_lists_down

  universe_symbols <- intersect(rownames(exprs_combat), keys(org.Hs.eg.db, keytype = "SYMBOL"))
  univ_id <- try(
    clusterProfiler::bitr(
      unique(universe_symbols),
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ),
    silent = TRUE
  )
  if (inherits(univ_id, "try-error")) univ_id <- NULL

  for (dir_nm in names(ora_sets)) {
    sig_src <- ora_sets[[dir_nm]]

    for (nm in names(sig_src)) {
      sigg <- unique(sig_src[[nm]])
      sigg <- sigg[!is.na(sigg)]
      sigg <- intersect(sigg, universe_symbols)

      if (length(sigg) < 10) next

      ego <- try(
        clusterProfiler::enrichGO(
          gene = sigg,
          universe = universe_symbols,
          OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "BP",
          pAdjustMethod = "BH",
          qvalueCutoff = 0.05,
          readable = TRUE
        ),
        silent = TRUE
      )
      if (inherits(ego, "try-error")) ego <- NULL

      if (!is.null(ego)) {
        ego_df <- as.data.frame(ego)
        if (nrow(ego_df) > 0) {
          write.csv(
            ego_df,
            file.path(TABLE_DIR, paste0("GO_BP_", dir_nm, "_", nm, ".csv")),
            row.names = FALSE
          )
          if (has_enrichplot) {
            p_go <- enrichplot::dotplot(ego, showCategory = 10) +
              ggplot2::ggtitle(paste0("GO BP ", dir_nm, " (", nm, ")"))
            try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, paste0("GO_BP_", dir_nm, "_", nm, ".png")), plot = p_go, width = 8, height = 5, dpi = PLOT_DPI), silent = TRUE)
            try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, paste0("GO_BP_", dir_nm, "_", nm, ".pdf")), plot = p_go, width = 8, height = 5, device = "pdf"), silent = TRUE)
          }
        }
      }

      idmap <- try(
        clusterProfiler::bitr(
          sigg,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = org.Hs.eg.db
        ),
        silent = TRUE
      )
      if (inherits(idmap, "try-error")) idmap <- NULL

      if (!is.null(idmap) && nrow(idmap) >= 10 && !is.null(univ_id) && nrow(univ_id) >= 10) {
        ek <- try(
          clusterProfiler::enrichKEGG(
            gene = unique(idmap$ENTREZID),
            universe = unique(univ_id$ENTREZID),
            organism = "hsa",
            pAdjustMethod = "BH",
            qvalueCutoff = 0.05
          ),
          silent = TRUE
        )
        if (inherits(ek, "try-error")) ek <- NULL

        if (!is.null(ek)) {
          ek_df <- as.data.frame(ek)
          if (nrow(ek_df) > 0) {
            write.csv(
              ek_df,
              file.path(TABLE_DIR, paste0("KEGG_", dir_nm, "_", nm, ".csv")),
              row.names = FALSE
            )
            if (has_enrichplot) {
              p_kg <- enrichplot::dotplot(ek, showCategory = 10) +
                ggplot2::ggtitle(paste0("KEGG ", dir_nm, " (", nm, ")"))
              try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, paste0("KEGG_", dir_nm, "_", nm, ".png")), plot = p_kg, width = 8, height = 5, dpi = PLOT_DPI), silent = TRUE)
              try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, paste0("KEGG_", dir_nm, "_", nm, ".pdf")), plot = p_kg, width = 8, height = 5, device = "pdf"), silent = TRUE)
            }
          }
        }
      }
    }
  }
}

if (has_gsva && !is.null(msig_h)) {
  gs_focus <- msig_h[intersect(names(msig_h), hallmark_focus)]
  if (length(gs_focus) > 0) {
    gsva_param <- NULL
    ss_res <- NULL

    if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
      gsva_param <- try(GSVA::ssgseaParam(exprData = as.matrix(exprs_combat), geneSets = gs_focus), silent = TRUE)
      if (inherits(gsva_param, "try-error")) gsva_param <- NULL
      if (!is.null(gsva_param)) {
        ss_res <- try(GSVA::gsva(gsva_param, verbose = FALSE), silent = TRUE)
        if (inherits(ss_res, "try-error")) ss_res <- NULL
      }
    }

    if (is.null(ss_res)) {
      ss_res <- try(
        GSVA::gsva(as.matrix(exprs_combat), gs_focus, method = "ssgsea", verbose = FALSE),
        silent = TRUE
      )
      if (inherits(ss_res, "try-error")) ss_res <- NULL
    }

    if (!is.null(ss_res)) {
      saveRDS(ss_res, file.path(OUTPUT_DIR, "ssgsea_hallmark_focus.rds"))
      write.csv(ss_res, file.path(TABLE_DIR, "ssgsea_hallmark_focus.csv"))

      ss_z <- as.matrix(ss_res)
      ss_z <- t(scale(t(ss_z)))
      ss_z[!is.finite(ss_z)] <- 0
      ann_col <- data.frame(Cluster = groups)
      rownames(ann_col) <- colnames(ss_z)

      try({
        png(file.path(PLOT_PNG_DIR, "ssGSEA_hallmark_focus_heatmap.png"), width = 8, height = 7, units = "in", res = PLOT_DPI)
        pheatmap::pheatmap(
          ss_z,
          main = "Focused stromal/inflammatory Hallmark ssGSEA (z-score)",
          annotation_col = ann_col,
          show_rownames = TRUE,
          show_colnames = FALSE
        )
        dev.off()
      }, silent = TRUE)

      try({
        pdf(file.path(PLOT_PDF_DIR, "ssGSEA_hallmark_focus_heatmap.pdf"), width = 8, height = 7)
        pheatmap::pheatmap(
          ss_z,
          main = "Focused stromal/inflammatory Hallmark ssGSEA (z-score)",
          annotation_col = ann_col,
          show_rownames = TRUE,
          show_colnames = FALSE
        )
        dev.off()
      }, silent = TRUE)

      ss_long <- as.data.frame(t(ss_res))
      ss_long$SampleID <- rownames(ss_long)
      ss_long$Cluster <- groups
      ss_long <- tidyr::pivot_longer(
        ss_long,
        cols = setdiff(names(ss_long), c("SampleID", "Cluster")),
        names_to = "GeneSet",
        values_to = "Score"
      )
      write.csv(ss_long, file.path(TABLE_DIR, "ssgsea_hallmark_focus_long.csv"), row.names = FALSE)

      p <- ggplot2::ggplot(ss_long, ggplot2::aes(Cluster, Score, fill = Cluster)) +
        ggplot2::geom_boxplot(outlier.size = 0.3) +
        ggplot2::facet_wrap(~ GeneSet, scales = "free_y") +
        ggplot2::theme_bw() +
        ggplot2::guides(fill = "none") +
        ggplot2::labs(title = "Focused stromal/inflammatory Hallmark ssGSEA by cluster", x = "Cluster", y = "ssGSEA score")
      try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, "ssGSEA_hallmark_focus_boxplots.png"), plot = p, width = 11, height = 8, dpi = PLOT_DPI), silent = TRUE)
      try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, "ssGSEA_hallmark_focus_boxplots.pdf"), plot = p, width = 11, height = 8, device = "pdf"), silent = TRUE)
    }
  }
}

expr_sig <- exprs_combat[signature_genes, , drop = FALSE]
centroids <- matrix(NA_real_, nrow = nrow(expr_sig), ncol = length(levels(groups)), dimnames = list(rownames(expr_sig), levels(groups)))
for (i_cent in seq_along(levels(groups))) {
  cl <- levels(groups)[i_cent]
  centroids[, i_cent] <- rowMeans(expr_sig[, groups == cl, drop = FALSE])
}

common <- intersect(rownames(centroids), rownames(expr_sig))
X <- expr_sig[common, , drop = FALSE]
C <- centroids[common, , drop = FALSE]
cors <- cor(X, C, method = "pearson")
pred <- character(nrow(cors))
for (i_pred in seq_len(nrow(cors))) pred[i_pred] <- colnames(C)[which.max(cors[i_pred, ])]
pred_obj <- list(cors = cors, pred = pred)
cm <- table(True = groups, Pred = factor(pred_obj$pred, levels = levels(groups)))
write.csv(as.data.frame.matrix(cm), file.path(TABLE_DIR, "in_sample_centroid_confusion_matrix.csv"))

if (RUN_NETWORK && length(signature_genes) >= 5) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Running cluster-adjusted coexpression network"))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
  E_all <- exprs_combat[signature_genes, , drop = FALSE]
  rho_all <- suppressWarnings(cor(t(E_all), method = "spearman", use = "pairwise.complete.obs"))

  groups_resid <- as.factor(groups)
  X_resid <- model.matrix(~ groups_resid)
  E_res <- E_all
  for (i_resid in seq_len(nrow(E_all))) {
    fit_resid <- try(stats::lm.fit(x = X_resid, y = as.numeric(E_all[i_resid, ])), silent = TRUE)
    if (inherits(fit_resid, "try-error")) fit_resid <- NULL
    if (is.null(fit_resid)) {
      E_res[i_resid, ] <- NA_real_
    } else {
      E_res[i_resid, ] <- fit_resid$residuals
    }
  }
  rho_res <- suppressWarnings(cor(t(E_res), method = "spearman", use = "pairwise.complete.obs"))

  rho_by_cluster <- list()
  for (cl in levels(groups)) {
    idx <- which(groups == cl)
    if (length(idx) >= 4) {
      rho_by_cluster[[cl]] <- suppressWarnings(cor(
        t(E_all[, idx, drop = FALSE]),
        method = "spearman",
        use = "pairwise.complete.obs"
      ))
    } else {
      rho_by_cluster[[cl]] <- matrix(
        NA_real_, nrow = nrow(E_all), ncol = nrow(E_all),
        dimnames = list(rownames(E_all), rownames(E_all))
      )
    }
  }

  write.csv(rho_all, file.path(TABLE_DIR, "network_overall_spearman.csv"))
  write.csv(rho_res, file.path(TABLE_DIR, "network_cluster_residual_spearman.csv"))
  for (cl in names(rho_by_cluster)) {
    write.csv(rho_by_cluster[[cl]], file.path(TABLE_DIR, paste0("network_within_", cl, "_spearman.csv")))
  }

  gnames <- rownames(E_all)
  idx <- which(upper.tri(rho_all), arr.ind = TRUE)
  edge_df <- data.frame(
    from = gnames[idx[, 1]],
    to   = gnames[idx[, 2]],
    overall_rho = rho_all[idx],
    residual_rho = rho_res[idx],
    stringsAsFactors = FALSE
  )

  for (cl in names(rho_by_cluster)) {
    edge_df[[paste0("rho_", cl)]] <- rho_by_cluster[[cl]][idx]
  }

  edge_df$sign_agree_n <- 0L
  rho_cols <- c("residual_rho", paste0("rho_", levels(groups)))
  rho_cols <- intersect(rho_cols, names(edge_df))

  for (i in seq_len(nrow(edge_df))) {
    vals <- as.numeric(edge_df[i, rho_cols, drop = TRUE])
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      pos_n <- sum(vals > 0)
      neg_n <- sum(vals < 0)
      edge_df$sign_agree_n[i] <- max(pos_n, neg_n)
    }
  }

  edge_df$keep_robust <- with(edge_df,
                              abs(overall_rho) >= NET_OVERALL_THR &
                                abs(residual_rho) >= NET_RESIDUAL_THR &
                                sign_agree_n >= 2)

  write.csv(edge_df, file.path(TABLE_DIR, "network_edge_metrics_cluster_adjusted.csv"), row.names = FALSE)

  robust_edges <- subset(edge_df, keep_robust)
  write.csv(robust_edges, file.path(TABLE_DIR, "network_edges_robust_cluster_adjusted.csv"), row.names = FALSE)

  if (nrow(robust_edges) > 0 && has_igraph) {
    g_net <- igraph::graph_from_data_frame(
      robust_edges[, c("from", "to", "overall_rho", "residual_rho"), drop = FALSE],
      directed = FALSE
    )
    net_cent <- data.frame(
      gene = igraph::V(g_net)$name,
      degree = igraph::degree(g_net),
      betweenness = igraph::betweenness(g_net),
      closeness = igraph::closeness(g_net),
      stringsAsFactors = FALSE
    ) %>% dplyr::arrange(dplyr::desc(degree), dplyr::desc(betweenness))
    write.csv(net_cent, file.path(TABLE_DIR, "network_centrality_cluster_adjusted.csv"), row.names = FALSE)

    if (has_ggraph) {
      p_net <- ggraph::ggraph(g_net, layout = "fr") +
        ggraph::geom_edge_link(ggplot2::aes(width = abs(residual_rho), color = residual_rho), alpha = 0.4) +
        ggraph::geom_node_point(size = 4) +
        ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3) +
        ggraph::scale_edge_colour_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = "Robust coexpression network after cluster adjustment")
      try(ggplot2::ggsave(file.path(PLOT_PNG_DIR, "network_cluster_adjusted_robust.png"), plot = p_net, width = 10, height = 9, dpi = PLOT_DPI), silent = TRUE)
      try(ggplot2::ggsave(file.path(PLOT_PDF_DIR, "network_cluster_adjusted_robust.pdf"), plot = p_net, width = 10, height = 9, device = "pdf"), silent = TRUE)
    }
  }
}
