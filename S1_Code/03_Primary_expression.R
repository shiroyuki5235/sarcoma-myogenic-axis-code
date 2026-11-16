# Sarcoma primary expression analysis and reviewer-requested sensitivities
# Historical prepared matrices are restored from repository parts or the companion archive.
# Public TCGA inputs are downloaded automatically. Online RMA matrices must not replace primary GEO inputs.
# Set PROJECT_DIR once; stage outputs are connected automatically within this script.
# Run the numbered modules in order. No source() or RUN_ALL wrapper is used.

#Module 00. Environment and paths ----
options(stringsAsFactors = FALSE, timeout = 1200, warn = 1)
CODE_DIR <- normalizePath(Sys.getenv("SARCOMA_CODE_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
PROJECT_DIR <- Sys.getenv("SARCOMA_PROJECT_DIR", unset = file.path(CODE_DIR, "Inputs", "Primary_project"))
LEGACY_OUTPUT_NAME <- "Output_20260402_papercode"
LEGACY_DIR <- file.path(PROJECT_DIR, LEGACY_OUTPUT_NAME)
CONFIG_DIR <- file.path(CODE_DIR, "Inputs", "Revision_config")
HALLMARK_FILE <- file.path(PROJECT_DIR, "Input", "GeneSets", "Hallmark_source_table.rds")
INSTALL_MISSING_PACKAGES <- TRUE
USE_PREPARED_MATRICES <- TRUE; ALLOW_PREFILTERED_FALLBACK <- FALSE
BATCH_METHOD <- "ComBat"; LOW_VARIANCE_QUANTILE <- 0.05
TOP_VARIABLE_GENES <- 3000L; NUMBER_OF_PCS <- 30L; K_GRID <- 2:6
SUBSAMPLE_REPLICATES <- 100L; SUBSAMPLE_FRACTION <- 0.80
SIGNATURE_SIZE <- 30L; MIN_SIGNATURE_GENES <- 5L; MIN_EVENTS_PER_PARAMETER <- 10
AUC_BOOTSTRAPS <- 2000L; RUN_CLUSTERING <- TRUE; RUN_HALLMARK <- TRUE
PLOT_DPI <- 1200L; SAVE_LARGE_MATRICES <- TRUE
SOURCE_PRIORITY <- c("GSE71118", "GSE21050", "GSE16382", "GSE34620", "GSE142162", "GSE66533", "GSE20196", "GSE34800", "GSE94321")
metadata_file <- file.path(LEGACY_DIR, "tables", "sample_master_metadata_curated_with_cluster.csv")
if (!file.exists(metadata_file)) stop("Set PROJECT_DIR to the historical primary project containing Input and Output_20260402_papercode.")
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
PRIMARY_RESULTS <- file.path(CODE_DIR, "Results", paste0("Primary_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (dir.exists(PRIMARY_RESULTS)) stop("Results directory already exists; use a fresh run.")
dir.create(PRIMARY_RESULTS, recursive = TRUE)
log_message <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste(..., collapse = " "))
save_figure <- function(plot_object, filename, width = 8, height = 5) {
  ggplot2::ggsave(file.path(PLOT_DIR, paste0(filename, ".pdf")), plot_object, width = width, height = height, bg = "white")
}

#Module 01. Packages ----
cran_packages <- c("data.table", "matrixStats", "ggplot2", "igraph", "mclust", "cluster", "survival", "digest", "cowplot")
bioc_packages <- c("sva", "limma", "fgsea")
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (INSTALL_MISSING_PACKAGES && length(missing_cran)) install.packages(missing_cran, repos = "https://cloud.r-project.org")
if (INSTALL_MISSING_PACKAGES && length(missing_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
required_packages <- c(cran_packages, bioc_packages)
stopifnot(all(vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)))
library(survival); data.table::setDTthreads(2L)

#Module 01a. Exact historical inputs and frozen public TCGA inputs ----
# A full repository/Zenodo download contains Data/Primary beside S1_Code.
# S1 Code alone downloads the same parts from the v1.1.0 GitHub tag after publication.
PRIMARY_INPUT_ARCHIVE <- Sys.getenv("SARCOMA_PRIMARY_INPUT_ARCHIVE",
  unset = file.path(dirname(CODE_DIR), "Sarcoma_Primary_Inputs_20261116.zip"))
PRIMARY_INPUT_URL <- Sys.getenv("SARCOMA_PRIMARY_INPUT_URL", unset = "")
PRIMARY_PARTS_BASE_URL <- Sys.getenv("SARCOMA_PRIMARY_PARTS_BASE_URL", unset = "")
PRIMARY_INPUT_SHA256 <- "e4e9e5eabe550945b48b698172e25e488276aec1ee26951faf447213736f9b94"
required_inputs <- read.csv(file.path(CODE_DIR, "Inputs", "Required_input_files.csv"),
  check.names = FALSE, colClasses = "character", na.strings = NULL)
stopifnot(nrow(required_inputs) == 12L, !anyDuplicated(required_inputs$relative_path),
  all(grepl("^[a-f0-9]{64}$", required_inputs$sha256)))
historical_inputs <- required_inputs[required_inputs$source == "historical_archive", , drop = FALSE]
historical_paths <- file.path(PROJECT_DIR, historical_inputs$relative_path)
if (any(!file.exists(historical_paths))) {
  if (!file.exists(PRIMARY_INPUT_ARCHIVE)) {
    dir.create(dirname(PRIMARY_INPUT_ARCHIVE), recursive = TRUE, showWarnings = FALSE)
    primary_temporary <- tempfile("primary_", tmpdir = dirname(PRIMARY_INPUT_ARCHIVE), fileext = ".download")
    primary_output <- NULL
    primary_error <- tryCatch({
      if (nzchar(PRIMARY_INPUT_URL)) {
        # Optional direct HTTPS mirror of the exact complete archive.
        if (!grepl("^https://", PRIMARY_INPUT_URL)) stop("The primary input URL must use HTTPS.")
        utils::download.file(PRIMARY_INPUT_URL, primary_temporary, method = "libcurl", mode = "wb")
      } else {
        primary_parts <- read.csv(file.path(CODE_DIR, "Inputs", "Primary_archive_parts.csv"),
          check.names = FALSE, colClasses = "character", na.strings = NULL)
        stopifnot(nrow(primary_parts) == 10L,
          identical(as.integer(primary_parts$part), seq_len(10L)),
          !anyDuplicated(primary_parts$relative_path),
          all(grepl("^Data/Primary/Sarcoma_Primary_Inputs_20261116[.]zip[.]part[0-9]{3}$", primary_parts$relative_path)),
          all(grepl("^[a-f0-9]{64}$", primary_parts$sha256)))
        primary_output <- file(primary_temporary, open = "wb")
        for (part_index in seq_len(nrow(primary_parts))) {
          primary_local_part <- file.path(dirname(CODE_DIR), primary_parts$relative_path[part_index])
          if (!file.exists(primary_local_part)) {
            primary_cache_dir <- file.path(CODE_DIR, "Inputs", "Download_cache", "Primary_parts")
            dir.create(primary_cache_dir, recursive = TRUE, showWarnings = FALSE)
            primary_local_part <- file.path(primary_cache_dir, basename(primary_parts$relative_path[part_index]))
            if (!file.exists(primary_local_part)) {
              primary_part_url <- primary_parts$url[part_index]
              if (nzchar(PRIMARY_PARTS_BASE_URL)) primary_part_url <- paste0(
                sub("/+$", "", PRIMARY_PARTS_BASE_URL), "/", basename(primary_parts$relative_path[part_index]))
              if (!grepl("^https://", primary_part_url)) stop("Primary-part URLs must use HTTPS.")
              primary_part_temporary <- tempfile("part_", tmpdir = primary_cache_dir, fileext = ".download")
              part_error <- tryCatch({
                message("Downloading historical input part ", part_index, "/", nrow(primary_parts))
                utils::download.file(primary_part_url, primary_part_temporary, method = "libcurl", mode = "wb")
                if (file.info(primary_part_temporary)$size != as.numeric(primary_parts$bytes[part_index]) ||
                    !identical(digest::digest(file = primary_part_temporary, algo = "sha256"), primary_parts$sha256[part_index]))
                  stop("Downloaded primary part has an unexpected size or checksum.")
                if (!file.rename(primary_part_temporary, primary_local_part)) stop("Could not install the downloaded part.")
                NULL
              }, error = function(e) conditionMessage(e))
              if (file.exists(primary_part_temporary)) unlink(primary_part_temporary)
              if (!is.null(part_error)) stop("Primary input download failed: ", part_error,
                "\nConfirm that GitHub release v1.1.0 is published and reachable, or use the full repository/Zenodo archive containing Data/Primary.")
            }
          }
          if (file.info(primary_local_part)$size != as.numeric(primary_parts$bytes[part_index]) ||
              !identical(digest::digest(file = primary_local_part, algo = "sha256"), primary_parts$sha256[part_index]))
            stop("Historical input part is damaged or differs from the release: ", primary_local_part)
          primary_raw <- readBin(primary_local_part, what = "raw", n = as.integer(primary_parts$bytes[part_index]))
          stopifnot(length(primary_raw) == as.integer(primary_parts$bytes[part_index]))
          writeBin(primary_raw, primary_output); rm(primary_raw)
        }
        close(primary_output); primary_output <- NULL
      }
      if (!identical(digest::digest(file = primary_temporary, algo = "sha256"), PRIMARY_INPUT_SHA256))
        stop("Historical input archive checksum mismatch; the file was not installed.")
      if (!file.rename(primary_temporary, PRIMARY_INPUT_ARCHIVE)) stop("Could not install the historical input archive.")
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(primary_output)) close(primary_output)
    if (file.exists(primary_temporary)) unlink(primary_temporary)
    if (!is.null(primary_error)) stop(primary_error)
  }
  if (!identical(digest::digest(file = PRIMARY_INPUT_ARCHIVE, algo = "sha256"), PRIMARY_INPUT_SHA256))
    stop("Historical input archive checksum mismatch; use the exact companion archive.")
  primary_members <- utils::unzip(PRIMARY_INPUT_ARCHIVE, list = TRUE)$Name
  if (!setequal(primary_members, historical_inputs$relative_path) || anyDuplicated(primary_members))
    stop("Historical input archive contains unexpected or missing members.")
  missing_members <- historical_inputs$relative_path[!file.exists(historical_paths)]
  utils::unzip(PRIMARY_INPUT_ARCHIVE, files = missing_members, exdir = PROJECT_DIR)
}
for (input_index in seq_len(nrow(historical_inputs))) {
  if (!file.exists(historical_paths[input_index]) ||
      !identical(digest::digest(file = historical_paths[input_index], algo = "sha256"), historical_inputs$sha256[input_index]))
    stop("Historical input is absent or differs from the archived version: ", historical_inputs$relative_path[input_index])
}


# Three exact files from a fixed cBioPortal commit; no current-release substitution.
tcga_sources <- required_inputs[required_inputs$source == "cbioportal_fixed_commit", , drop = FALSE]
for (input_index in seq_len(nrow(tcga_sources))) {
  tcga_destination <- file.path(PROJECT_DIR, tcga_sources$relative_path[input_index])
  dir.create(dirname(tcga_destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(tcga_destination)) {
    tcga_part <- tempfile("tcga_", tmpdir = dirname(tcga_destination), fileext = ".part")
    tcga_error <- tryCatch({
      utils::download.file(tcga_sources$url[input_index], tcga_part, method = "libcurl", mode = "wb")
      if (!identical(digest::digest(file = tcga_part, algo = "sha256"), tcga_sources$sha256[input_index]))
        stop("Downloaded TCGA checksum mismatch; the file was not installed.")
      if (!file.rename(tcga_part, tcga_destination)) stop("Could not install the TCGA file.")
      NULL
    }, error = function(e) conditionMessage(e))
    if (file.exists(tcga_part)) unlink(tcga_part)
    if (!is.null(tcga_error)) stop(tcga_error)
  }
  if (!identical(digest::digest(file = tcga_destination, algo = "sha256"), tcga_sources$sha256[input_index]))
    stop("TCGA input differs from the frozen analysis input: ", tcga_sources$relative_path[input_index])
}
message("Inputs ready: nine exact historical GEO matrices and three fixed-release TCGA files.")

#Module 02. 01_Record_preparation output paths ----
SEED <- 20261023L; set.seed(SEED)
OUTPUT_DIR <- file.path(PRIMARY_RESULTS, "01_Record_preparation")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables"); PLOT_DIR <- file.path(OUTPUT_DIR, "figures")
CACHE_DIR <- file.path(OUTPUT_DIR, "cache"); LOG_DIR <- file.path(OUTPUT_DIR, "session")
for (directory in c(TABLE_DIR, PLOT_DIR, CACHE_DIR, LOG_DIR)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
input_files <- character()
#Module 03. Metadata and ID checks ---------------------------------------------
metadata_file <- file.path(LEGACY_DIR, "tables", "sample_master_metadata_curated_with_cluster.csv")
if (!file.exists(metadata_file)) stop("Missing metadata: ", metadata_file)
input_files <- c(input_files, metadata_file)
metadata <- data.table::fread(metadata_file, data.table = FALSE, colClasses = "character",
  na.strings = c("", "NA", "N/A", "null"))
required_columns <- c("sampleid", "gse_id", "geo_accession", "expr_col_original",
  "diagnosis_group", "platform_id")
stopifnot(all(required_columns %in% names(metadata)), nrow(metadata) > 30L)
if (anyNA(metadata$sampleid) || anyDuplicated(metadata$sampleid) ||
    anyNA(metadata$geo_accession) || anyDuplicated(metadata$geo_accession)) {
  stop("Missing or duplicated sampleid/GSM. Resolve IDs before analysis.")
}
if ("sample_decision" %in% names(metadata)) {
  stopifnot(all(metadata$sample_decision == "KEEP"))
}
expression_gsm <- regmatches(metadata$expr_col_original,
  regexpr("GSM[0-9]+", metadata$expr_col_original))
if (length(expression_gsm) != nrow(metadata) ||
    !identical(expression_gsm, metadata$geo_accession)) {
  stop("Expression-column GSM and metadata GSM do not agree. Row-order matching is prohibited.")
}
metadata$legacy_cluster <- if ("Cluster" %in% names(metadata)) metadata$Cluster else NA_character_
metadata$Cluster <- NULL
for (field in c("description", "relation", "relation_1", "time", "metastasis",
                "age", "sex", "grade_h", "stage_h", "tissue_group")) {
  if (!field %in% names(metadata)) metadata[[field]] <- NA_character_
}
# Confirmed fusion annotations from GEO GSE34800 override title keywords.
bcor_gsm <- paste0("GSM", 855583:855592)
ewing_control_gsm <- paste0("GSM", 855593:855596)
stopifnot(sum(metadata$geo_accession %in% bcor_gsm) == 10L,
          sum(metadata$geo_accession %in% ewing_control_gsm) == 4L)
metadata$diagnosis_group[metadata$geo_accession %in% bcor_gsm] <- "BCOR-CCNB3 sarcoma"
metadata$diagnosis_group[metadata$geo_accession %in% ewing_control_gsm] <- "Ewing sarcoma"

log_message("Metadata:", nrow(metadata), "records; unique sample/GSM IDs; GSM alignment verified.")

#Module 04. Explicit reanalysis links and representative records ---------------
# Links are extracted only from fields explicitly containing reanalysis wording.
link_rows <- list()
link_index <- 0L
for (row_index in seq_len(nrow(metadata))) {
  for (field in c("description", "relation", "relation_1")) {
    field_text <- metadata[[field]][row_index]
    if (is.na(field_text) || !grepl("reanaly[sz]", field_text, ignore.case = TRUE)) next
    target_ids <- unique(regmatches(field_text, gregexpr("GSM[0-9]+", field_text))[[1]])
    target_ids <- setdiff(target_ids, metadata$geo_accession[row_index])
    for (target_id in target_ids) {
      link_index <- link_index + 1L
      link_rows[[link_index]] <- data.frame(from = metadata$geo_accession[row_index],
        to = target_id, source_field = field, evidence = field_text)
    }
  }
}
links <- if (length(link_rows)) do.call(rbind, link_rows) else
  data.frame(from = character(), to = character(), source_field = character(), evidence = character())
manual_links_file <- file.path(CONFIG_DIR, "Additional_verified_links.csv")
if (file.exists(manual_links_file)) {
  manual_links <- data.table::fread(manual_links_file, data.table = FALSE)
  stopifnot(all(c("from", "to", "evidence", "verified") %in% names(manual_links)))
  manual_links <- manual_links[toupper(as.character(manual_links$verified)) == "TRUE", , drop = FALSE]
  if (nrow(manual_links)) {
    stopifnot(all(nzchar(manual_links$evidence)), !anyNA(manual_links$evidence),
      all(manual_links$from %in% metadata$geo_accession), all(manual_links$to %in% metadata$geo_accession))
    links <- rbind(links, data.frame(from = manual_links$from, to = manual_links$to,
      source_field = "manual_verified", evidence = manual_links$evidence))
  }
  input_files <- c(input_files, manual_links_file)
}
links$target_in_curated_data <- links$to %in% metadata$geo_accession
write.csv(links, file.path(TABLE_DIR, "03_reanalysis_link_evidence.csv"), row.names = FALSE)
present_links <- links[links$target_in_curated_data, c("from", "to"), drop = FALSE]
edges <- unique(data.frame(from = pmin(present_links$from, present_links$to),
  to = pmax(present_links$from, present_links$to)))
link_graph <- igraph::graph_from_data_frame(edges, directed = FALSE,
  vertices = data.frame(name = metadata$geo_accession))
membership <- igraph::components(link_graph)$membership[metadata$geo_accession]
metadata$record_group <- NA_character_
for (component in unique(membership)) {
  group_rows <- which(membership == component)
  metadata$record_group[group_rows] <- sort(metadata$geo_accession[group_rows])[1]
}
metadata$priority <- match(metadata$gse_id, SOURCE_PRIORITY)
metadata$priority[is.na(metadata$priority)] <- length(SOURCE_PRIORITY) + 1L
metadata$representative <- FALSE
metadata$selection_reason <- "not_selected"
for (group_id in unique(metadata$record_group)) {
  group_rows <- which(metadata$record_group == group_id)
  chosen <- group_rows[order(metadata$priority[group_rows], metadata$sampleid[group_rows])][1]
  metadata$representative[chosen] <- TRUE
  metadata$selection_reason[chosen] <- "prespecified_source_priority_then_sampleid"
}
override_file <- file.path(CONFIG_DIR, "Representative_overrides.csv")
if (file.exists(override_file)) {
  overrides <- data.table::fread(override_file, data.table = FALSE)
  stopifnot(all(c("record_group", "sampleid", "reason") %in% names(overrides)))
  if (nrow(overrides)) {
    stopifnot(!anyDuplicated(overrides$record_group), !anyNA(overrides$reason), all(nzchar(overrides$reason)))
    for (override_index in seq_len(nrow(overrides))) {
      group_rows <- which(metadata$record_group == overrides$record_group[override_index])
      chosen <- which(metadata$sampleid == overrides$sampleid[override_index])
      if (length(chosen) != 1L || !chosen %in% group_rows) stop("Invalid representative override.")
      metadata$representative[group_rows] <- FALSE
      metadata$representative[chosen] <- TRUE
      metadata$selection_reason[group_rows] <- "manual_override_other_record_selected"
      metadata$selection_reason[chosen] <- overrides$reason[override_index]
    }
  }
  input_files <- c(input_files, override_file)
}
conflicts <- list()
for (field in c("diagnosis_group", "tissue_group", "time", "metastasis", "age", "sex", "grade_h", "stage_h")) {
  for (group_id in unique(metadata$record_group)) {
    group_rows <- which(metadata$record_group == group_id)
    values <- unique(tolower(trimws(stats::na.omit(metadata[[field]][group_rows]))))
    values <- values[nzchar(values)]
    if (length(values) > 1L) {
      conflicts[[length(conflicts) + 1L]] <- data.frame(record_group = group_id,
        field = field, values = paste(values, collapse = " | "),
        sampleids = paste(metadata$sampleid[group_rows], collapse = " | "))
    }
  }
}
conflict_table <- if (length(conflicts)) do.call(rbind, conflicts) else
  data.frame(record_group = character(), field = character(), values = character(), sampleids = character())
write.csv(conflict_table, file.path(TABLE_DIR, "03_annotation_conflicts.csv"), row.names = FALSE)
write.csv(metadata, file.path(TABLE_DIR, "03_record_lineage_and_selection.csv"), row.names = FALSE)
discovery <- metadata[metadata$representative, , drop = FALSE]
stopifnot(!anyDuplicated(discovery$record_group))
discovery$annotation_conflict <- discovery$record_group %in% conflict_table$record_group
annotation_file <- file.path(CONFIG_DIR, "Annotation_overrides.csv")
if (file.exists(annotation_file)) {
  annotations <- data.table::fread(annotation_file, data.table = FALSE, colClasses = "character")
  stopifnot(all(c("record_group", "field", "value", "reason") %in% names(annotations)))
  if (nrow(annotations)) {
    stopifnot(!anyDuplicated(paste(annotations$record_group, annotations$field)),
      all(annotations$record_group %in% discovery$record_group),
      all(annotations$field %in% c("diagnosis_group", "tissue_group", "age", "sex", "grade_h", "stage_h")),
      !anyNA(annotations$reason), all(nzchar(annotations$reason)))
    for (annotation_index in seq_len(nrow(annotations))) {
      row_index <- match(annotations$record_group[annotation_index], discovery$record_group)
      discovery[[annotations$field[annotation_index]]][row_index] <- annotations$value[annotation_index]
    }
  }
  input_files <- c(input_files, annotation_file)
  write.csv(annotations, file.path(TABLE_DIR, "03_applied_annotation_overrides.csv"), row.names = FALSE)
}
write.csv(data.frame(original_records = nrow(metadata), explicit_edges = nrow(edges),
  retained_record_groups = nrow(discovery), removed_records = nrow(metadata) - nrow(discovery),
  units = "explicitly linked record groups; not confirmed unique patients"),
  file.path(TABLE_DIR, "03_deduplication_summary.csv"), row.names = FALSE)
log_message("Deduplication:", nrow(metadata), "records ->", nrow(discovery), "record groups.")

#Module 05. Rebuild uncorrected gene expression --------------------------------
series_ids <- unique(metadata$gse_id)
prepared_files <- file.path(PROJECT_DIR, "Input", "GEO_prepared", series_ids,
  paste0(series_ids, "_ExpressionMatrix_Annotated.csv"))
can_rebuild <- USE_PREPARED_MATRICES && all(file.exists(prepared_files))
if (can_rebuild) {
  expression_list <- list()
  for (series_index in seq_along(series_ids)) {
    series_id <- series_ids[series_index]
    log_message("Reading prepared expression:", series_id)
    source_table <- data.table::fread(prepared_files[series_index], data.table = FALSE, check.names = FALSE)
    gene_column <- if ("Gene_Symbol" %in% names(source_table)) "Gene_Symbol" else names(source_table)[1]
    wanted_metadata <- metadata[metadata$gse_id == series_id, , drop = FALSE]
    wanted_columns <- wanted_metadata$expr_col_original
    if (!all(wanted_columns %in% names(source_table))) stop("Unmatched expression columns: ", series_id)
    gene_symbols <- trimws(as.character(source_table[[gene_column]]))
    expression_matrix <- as.matrix(source_table[, wanted_columns, drop = FALSE])
    storage.mode(expression_matrix) <- "double"
    if (anyNA(gene_symbols) || any(!nzchar(gene_symbols)) || anyDuplicated(gene_symbols)) {
      stop("Prepared expression must have one nonmissing row per gene symbol: ", series_id)
    }
    rownames(expression_matrix) <- gene_symbols
    colnames(expression_matrix) <- wanted_metadata$sampleid
    if (any(!is.finite(expression_matrix))) stop("Nonfinite prepared expression: ", series_id)
    expression_list[[series_id]] <- expression_matrix
  }
  shared_genes <- Reduce(intersect, lapply(expression_list, rownames))
  if (length(shared_genes) < 1000L) stop("Unexpectedly few shared gene symbols.")
  for (series_id in names(expression_list)) {
    expression_list[[series_id]] <- expression_list[[series_id]][shared_genes, , drop = FALSE]
  }
  expression_all <- do.call(cbind, expression_list)
  input_files <- c(input_files, prepared_files)
  feature_source <- "rebuilt_prepared_matrices_then_deduplicate_then_variance_filter"
  rm(expression_list, source_table, expression_matrix)
} else {
  if (!ALLOW_PREFILTERED_FALLBACK) {
    writeLines(prepared_files[!file.exists(prepared_files)], file.path(LOG_DIR, "missing_prepared_files.txt"))
    stop("Prepared inputs missing/disabled. See README before explicitly enabling ALLOW_PREFILTERED_FALLBACK.")
  }
  prebatch_file <- file.path(LEGACY_DIR, "exprs_all_curated_prefiltered.rds")
  expression_all <- readRDS(prebatch_file)
  input_files <- c(input_files, prebatch_file)
  feature_source <- "legacy_prebatch_matrix_fixed_original_filtered_gene_universe"
  log_message("LIMITATION: legacy prebatch matrix already excluded bottom 5% variance genes.")
}
stopifnot(is.matrix(expression_all), !anyDuplicated(rownames(expression_all)),
  !anyDuplicated(colnames(expression_all)), all(metadata$sampleid %in% colnames(expression_all)))
expression_all <- expression_all[, metadata$sampleid, drop = FALSE]
stopifnot(identical(colnames(expression_all), metadata$sampleid), all(is.finite(expression_all)))
expression_quantiles <- stats::quantile(expression_all, c(0, 0.01, 0.5, 0.99, 1))
write.csv(data.frame(quantile = names(expression_quantiles), value = as.numeric(expression_quantiles)),
  file.path(TABLE_DIR, "04_expression_scale_check.csv"), row.names = FALSE)
if (expression_quantiles[4] > 100 || expression_quantiles[1] < -100) {
  stop("Unexpected scale. Prepared inputs must already be log-scale normalized expression; no auto-log conversion.")
}
# Hashes flag exact duplicates; correlation or a shared name never authorizes collapsing records.
column_hash <- vapply(seq_len(ncol(expression_all)), function(column_index)
  digest::digest(expression_all[, column_index], algo = "xxhash64"), character(1))
duplicate_hash <- duplicated(column_hash) | duplicated(column_hash, fromLast = TRUE)
hash_table <- data.frame(sampleid = metadata$sampleid[duplicate_hash],
  record_group = metadata$record_group[duplicate_hash], hash = column_hash[duplicate_hash])
write.csv(hash_table, file.path(TABLE_DIR, "04_exact_expression_duplicate_candidates.csv"), row.names = FALSE)
unlinked_exact_duplicates <- 0L
for (hash_value in unique(hash_table$hash)) {
  rows <- which(column_hash == hash_value)
  same_values <- all(vapply(rows, function(row_index)
    identical(expression_all[, rows[1]], expression_all[, row_index]), logical(1)))
  if (same_values && length(unique(metadata$record_group[rows])) > 1L) {
    unlinked_exact_duplicates <- unlinked_exact_duplicates + 1L
  }
}
if (unlinked_exact_duplicates > 0L) {
  log_message("REVIEW:", unlinked_exact_duplicates, "exact-expression groups cross explicit record groups; not collapsed.")
}
expression_raw <- expression_all[, discovery$sampleid, drop = FALSE]
rm(expression_all)
gene_variance <- matrixStats::rowVars(expression_raw)
keep_genes <- is.finite(gene_variance) & gene_variance > 0 & rowSums(expression_raw != 0) > 0
if (can_rebuild) keep_genes <- keep_genes &
  gene_variance > stats::quantile(gene_variance[is.finite(gene_variance)], LOW_VARIANCE_QUANTILE)
write.csv(data.frame(GeneSymbol = rownames(expression_raw), variance = gene_variance, retained = keep_genes),
  file.path(TABLE_DIR, "04_gene_filter.csv"), row.names = FALSE)
expression_raw <- expression_raw[keep_genes, , drop = FALSE]
stopifnot(nrow(expression_raw) >= TOP_VARIABLE_GENES,
  identical(colnames(expression_raw), discovery$sampleid))
write.csv(data.frame(source = feature_source, genes = nrow(expression_raw),
  records = ncol(expression_raw)), file.path(TABLE_DIR, "04_expression_dimensions.csv"), row.names = FALSE)

#Module 06. Batch correction and batch/diagnosis confounding --------------------
discovery$batch <- factor(discovery$gse_id)
discovery$histology <- discovery$diagnosis_group
discovery$histology[is.na(discovery$histology) | !nzchar(discovery$histology)] <- "Unknown"
discovery$histology <- factor(discovery$histology)
write.csv(as.data.frame(table(discovery$batch, discovery$histology)),
  file.path(TABLE_DIR, "05_batch_by_histology.csv"), row.names = FALSE)
confounding_design <- stats::model.matrix(~ batch + histology, data = discovery)
histology_estimable <- qr(confounding_design)$rank == ncol(confounding_design)
write.csv(data.frame(columns = ncol(confounding_design), rank = qr(confounding_design)$rank,
  histology_and_batch_separable = histology_estimable),
  file.path(TABLE_DIR, "05_batch_histology_rank.csv"), row.names = FALSE)
log_message("Batch correction:", BATCH_METHOD, "; no clinical covariates in primary unsupervised correction.")
if (nlevels(discovery$batch) < 2L) {
  expression_corrected <- expression_raw
} else if (BATCH_METHOD == "ComBat") {
  expression_corrected <- sva::ComBat(expression_raw, batch = discovery$batch, mod = NULL,
    par.prior = TRUE, prior.plots = FALSE, mean.only = any(table(discovery$batch) == 1L))
} else {
  expression_corrected <- limma::removeBatchEffect(expression_raw, batch = discovery$batch)
}
stopifnot(all(is.finite(expression_corrected)), identical(dimnames(expression_raw), dimnames(expression_corrected)))
if (SAVE_LARGE_MATRICES) {
  saveRDS(expression_raw, file.path(CACHE_DIR, "discovery_deduplicated_uncorrected.rds"))
  saveRDS(expression_corrected, file.path(CACHE_DIR, "discovery_deduplicated_corrected.rds"))
}

#Module 07. Consistent clustering, k selection and subsampling ------------------
# Both selection and final labels use top-variable genes -> PCs -> Ward.D2.
cluster_genes <- head(rownames(expression_corrected)[
  order(matrixStats::rowVars(expression_corrected), decreasing = TRUE)], TOP_VARIABLE_GENES)
cluster_pca <- stats::prcomp(t(expression_corrected[cluster_genes, ]), center = TRUE,
  scale. = FALSE, rank. = min(NUMBER_OF_PCS, nrow(discovery) - 1L))
cluster_distance <- stats::dist(cluster_pca$x)
cluster_tree <- stats::hclust(cluster_distance, method = "ward.D2")
k_support <- data.frame(k = K_GRID, silhouette = NA_real_, mean_ARI = NA_real_,
  ARI_q025 = NA_real_, ARI_q975 = NA_real_, successful_subsamples = 0L)
for (k_index in seq_along(K_GRID)) {
  labels <- stats::cutree(cluster_tree, k = K_GRID[k_index])
  k_support$silhouette[k_index] <- mean(cluster::silhouette(labels, cluster_distance)[, "sil_width"])
}
selected_k <- k_support$k[order(-k_support$silhouette, k_support$k)][1]
selected_labels <- stats::cutree(cluster_tree, k = selected_k)
k2_labels <- stats::cutree(cluster_tree, k = 2L)

# External, prespecified smooth-muscle panel orients the two labels only.
# It does not choose k or optimize a survival split.
reference_panels <- list(
  Smooth_muscle_reference = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "LMOD1", "MYOCD", "DES"),
  Fibroblast_ECM_reference = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A1", "COL6A2"),
  T_cell_reference = c("CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK", "CD2"),
  Cytotoxic_reference = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "CTSW"),
  Myeloid_reference = c("LST1", "TYROBP", "FCER1G", "AIF1", "CSF1R", "CD68", "CTSS"),
  Endothelial_reference = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "ESAM", "CDH5"))
orientation_genes <- intersect(reference_panels$Smooth_muscle_reference, rownames(expression_corrected))
if (length(orientation_genes) < 5L) stop("Too few smooth-muscle reference genes to orient k=2 labels.")
orientation_z <- t(scale(t(expression_corrected[orientation_genes, , drop = FALSE])))
orientation_score <- colMeans(orientation_z)
orientation_means <- tapply(orientation_score, k2_labels, mean)
myogenic_label <- as.integer(names(orientation_means)[which.max(orientation_means)])
discovery$C1_C2 <- factor(ifelse(k2_labels == myogenic_label, "C2", "C1"), levels = c("C1", "C2"))
discovery$selected_cluster <- paste0("K", selected_labels)
write.csv(data.frame(original_k2_label = names(orientation_means),
  mean_reference_score = as.numeric(orientation_means),
  oriented_label = ifelse(as.integer(names(orientation_means)) == myogenic_label, "C2", "C1")),
  file.path(TABLE_DIR, "07_k2_label_orientation.csv"), row.names = FALSE)
write.csv(data.frame(selected_k = selected_k, k2_is_selected = selected_k == 2L,
  interpretation = if (selected_k == 2L) "k2_supported_by_maximum_silhouette" else
    "k2_is_prespecified_sensitivity_not_best_supported_partition"),
  file.path(TABLE_DIR, "07_selected_k.csv"), row.names = FALSE)
if (selected_k != 2L) log_message("REVIEW: silhouette selected k=", selected_k,
  "; C1/C2 is a prespecified two-group contrast, not the selected partition.")
write.csv(discovery, file.path(TABLE_DIR, "07_discovery_metadata_revised.csv"), row.names = FALSE)
#Module 08. Source-family and correction sensitivities -------------------------
series_edges <- data.frame(from = metadata$gse_id[match(edges$from, metadata$geo_accession)],
  to = metadata$gse_id[match(edges$to, metadata$geo_accession)])
series_graph <- igraph::graph_from_data_frame(series_edges, directed = FALSE,
  vertices = data.frame(name = unique(metadata$gse_id)))
series_membership <- igraph::components(series_graph)$membership
family_map <- data.frame(gse_id = names(series_membership), family = as.integer(series_membership))
write.csv(family_map, file.path(TABLE_DIR, "08_connected_series_families.csv"), row.names = FALSE)
stopifnot(nrow(metadata) == 1115L, nrow(discovery) == 629L,
  sum(discovery$legacy_cluster == "C2") == 37L)
V20_DIR <- OUTPUT_DIR

#Module 09. 02_Fixed_target output paths ----
SEED <- 20261024L; set.seed(SEED)
OUTPUT_DIR <- file.path(PRIMARY_RESULTS, "02_Fixed_target")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables"); PLOT_DIR <- file.path(OUTPUT_DIR, "figures")
CACHE_DIR <- file.path(OUTPUT_DIR, "cache"); LOG_DIR <- file.path(OUTPUT_DIR, "session")
for (directory in c(TABLE_DIR, PLOT_DIR, CACHE_DIR, LOG_DIR)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
input_files <- character()
PARENT_DIR <- V20_DIR
#Module 10. Load the completed run and verify IDs ------------------------------
parent_files <- file.path(PARENT_DIR, c("tables/07_discovery_metadata_revised.csv",
  "tables/03_record_lineage_and_selection.csv", "cache/discovery_deduplicated_uncorrected.rds",
  "cache/discovery_deduplicated_corrected.rds", "tables/07_selected_k.csv",
  "tables/08_connected_series_families.csv"))
if (!all(file.exists(parent_files))) stop("Missing parent input: ", paste(parent_files[!file.exists(parent_files)], collapse = "; "))
input_files <- c(input_files, parent_files)
discovery <- data.table::fread(parent_files[1], data.table = FALSE, check.names = FALSE)
metadata <- data.table::fread(parent_files[2], data.table = FALSE, check.names = FALSE)
expression_raw <- readRDS(parent_files[3]); expression_corrected <- readRDS(parent_files[4])
selected_k <- data.table::fread(parent_files[5])$selected_k[1]
family_map <- data.table::fread(parent_files[6], data.table = FALSE)
stopifnot(all(c("sampleid", "record_group", "legacy_cluster", "C1_C2", "selected_cluster", "gse_id", "histology") %in% names(discovery)),
  !anyNA(discovery$sampleid), !anyDuplicated(discovery$sampleid),
  !anyNA(discovery$record_group), !anyDuplicated(discovery$record_group),
  !anyDuplicated(colnames(expression_raw)), !anyDuplicated(rownames(expression_raw)),
  identical(dimnames(expression_raw), dimnames(expression_corrected)),
  setequal(colnames(expression_raw), discovery$sampleid),
  all(is.finite(expression_raw)), all(is.finite(expression_corrected)))
expression_raw <- expression_raw[, discovery$sampleid, drop = FALSE]
expression_corrected <- expression_corrected[, discovery$sampleid, drop = FALSE]
stopifnot(identical(colnames(expression_raw), discovery$sampleid))
metadata$priority <- as.numeric(metadata$priority)
discovery$batch <- factor(discovery$gse_id); discovery$histology <- factor(discovery$histology)
discovery$v2_k2_group <- paste0("v2_", discovery$C1_C2)
discovery$C1_C2 <- factor(discovery$legacy_cluster, levels = c("C1", "C2"))
if (anyNA(discovery$C1_C2) || any(table(discovery$C1_C2) < 3L)) stop("Original C1/C2 labels missing or too few records.")
# Ensure linked records did not have conflicting original target membership.
for (group_id in discovery$record_group) {
  linked_labels <- unique(stats::na.omit(metadata$legacy_cluster[metadata$record_group == group_id]))
  if (length(linked_labels) != 1L) stop("Conflicting original cluster labels in record group: ", group_id)
  if (linked_labels != as.character(discovery$C1_C2[discovery$record_group == group_id])) stop("Original label mismatch.")
}
write.csv(discovery, file.path(TABLE_DIR, "02_original_target_metadata.csv"), row.names = FALSE)
write.csv(as.data.frame(table(discovery$legacy_cluster, discovery$v2_k2_group)),
  file.path(TABLE_DIR, "02_original_vs_v2_k2.csv"), row.names = FALSE)
write.csv(as.data.frame(table(discovery$legacy_cluster, discovery$selected_cluster)),
  file.path(TABLE_DIR, "02_original_vs_selected_k.csv"), row.names = FALSE)
write.csv(as.data.frame(table(discovery$C1_C2, discovery$histology)),
  file.path(TABLE_DIR, "02_original_target_histology.csv"), row.names = FALSE)
discovery$source_family <- family_map$family[match(discovery$gse_id, family_map$gse_id)]
if (anyNA(discovery$source_family)) stop("Unmatched series-family mapping.")
write.csv(as.data.frame(table(discovery$C1_C2, discovery$source_family)),
  file.path(TABLE_DIR, "02_original_target_source_families.csv"), row.names = FALSE)
target_families <- unique(discovery$source_family[discovery$C1_C2 == "C2"])
write.csv(data.frame(target = "original_C2_record_membership_retained_after_deduplication",
  C1_n = sum(discovery$C1_C2 == "C1"), C2_n = sum(discovery$C1_C2 == "C2"),
  C2_source_families = length(target_families), parent_selected_k = selected_k,
  interpretation = "fixed_original_contrast_not_a_new_k2_partition_not_independent_discovery_validation"),
  file.path(TABLE_DIR, "02_analysis_target_definition.csv"), row.names = FALSE)
log_message("Loaded", nrow(expression_raw), "genes and", ncol(expression_raw), "records.")
log_message("Fixed original groups:", paste(names(table(discovery$C1_C2)), table(discovery$C1_C2), collapse = "; "))

#Module 11. Independent reference panel and label correspondence ---------------
reference_panels <- list(
  Smooth_muscle_reference = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "LMOD1", "MYOCD", "DES"),
  Fibroblast_ECM_reference = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A1", "COL6A2"),
  T_cell_reference = c("CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK", "CD2"),
  Cytotoxic_reference = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "CTSW"),
  Myeloid_reference = c("LST1", "TYROBP", "FCER1G", "AIF1", "CSF1R", "CD68", "CTSS"),
  Endothelial_reference = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "ESAM", "CDH5"))

reference_genes <- intersect(reference_panels$Smooth_muscle_reference, rownames(expression_corrected))
reference_genes <- reference_genes[matrixStats::rowVars(expression_corrected[reference_genes, , drop = FALSE]) > 0]
if (length(reference_genes) < 5L) stop("Insufficient external smooth-muscle reference coverage.")
smooth_reference <- colMeans(t(scale(t(expression_corrected[reference_genes, , drop = FALSE]))))
reference_check <- data.frame(sampleid = discovery$sampleid, original_group = discovery$C1_C2,
  revised_k2_group = discovery$v2_k2_group, selected_cluster = discovery$selected_cluster,
  smooth_reference = smooth_reference)
write.csv(reference_check, file.path(TABLE_DIR, "03_reference_score_by_original_and_revised_groups.csv"), row.names = FALSE)
write.csv(stats::aggregate(smooth_reference ~ selected_cluster, reference_check, mean),
  file.path(TABLE_DIR, "03_selected_cluster_smooth_reference_means.csv"), row.names = FALSE)
original_reference_difference <- mean(smooth_reference[discovery$C1_C2 == "C2"]) -
  mean(smooth_reference[discovery$C1_C2 == "C1"])
# Never exchange labels because a different group has a smaller P value or higher AUC.
write.csv(data.frame(reference = "prespecified_smooth_muscle_panel",
  mean_C2_minus_C1 = original_reference_difference, target_definition_changed_by_reference = FALSE),
  file.path(TABLE_DIR, "03_reference_direction_check.csv"), row.names = FALSE)

#Module 12. Fixed original C2 contrast: descriptive DE and signature -------
# Inference uses UNCORRECTED expression + batch, never ComBat + batch again.
# Cluster-derived DE is descriptive; its P values do not validate the clustering.
de_tables <- list()
de_status <- list()
for (de_scenario in c("all", "LMS_only")) {
  selected_rows <- if (de_scenario == "all") seq_len(nrow(discovery)) else
    which(discovery$histology == "Leiomyosarcoma")
  current_metadata <- droplevels(discovery[selected_rows, , drop = FALSE])
  if (nrow(current_metadata) < 10L || nlevels(current_metadata$C1_C2) < 2L ||
      min(table(current_metadata$C1_C2)) < 3L) {
    de_status[[de_scenario]] <- data.frame(scenario = de_scenario, status = "insufficient_group_size")
    next
  }
  design <- if (nlevels(current_metadata$batch) > 1L)
    stats::model.matrix(~ batch + C1_C2, current_metadata) else stats::model.matrix(~ C1_C2, current_metadata)
  if (qr(design)$rank < ncol(design)) {
    de_status[[de_scenario]] <- data.frame(scenario = de_scenario, status = "cluster_batch_confounding_nonestimable")
    next
  }
  fit <- limma::eBayes(limma::lmFit(expression_raw[, selected_rows, drop = FALSE], design), trend = TRUE)
  de_table <- limma::topTable(fit, coef = "C1_C2C2", number = Inf, adjust.method = "BH", sort.by = "P")
  de_table$GeneSymbol <- rownames(de_table)
  de_tables[[de_scenario]] <- de_table
  write.csv(de_table, file.path(TABLE_DIR, paste0("09_DE_C2_vs_C1_", de_scenario, ".csv")), row.names = FALSE)
  de_status[[de_scenario]] <- data.frame(scenario = de_scenario, status = "completed_descriptive")
}
write.csv(do.call(rbind, de_status), file.path(TABLE_DIR, "09_DE_status.csv"), row.names = FALSE)
signature_genes <- character()
if (!is.null(de_tables$all)) {
  signature_candidates <- de_tables$all[de_tables$all$logFC > 0 & de_tables$all$adj.P.Val < 0.05, , drop = FALSE]
  signature_candidates <- signature_candidates[order(-signature_candidates$t, signature_candidates$GeneSymbol), , drop = FALSE]
  signature_genes <- head(signature_candidates$GeneSymbol, SIGNATURE_SIZE)
}
write.csv(data.frame(GeneSymbol = signature_genes, rank = seq_along(signature_genes)),
  file.path(TABLE_DIR, "09_revised_discovery_signature.csv"), row.names = FALSE)
signature_ready <- length(signature_genes) >= MIN_SIGNATURE_GENES
log_message("Original-target discovery signature:", length(signature_genes), "genes; projection eligible:", signature_ready)


signature_reference_rho <- NA_real_
if (signature_ready) {
  signature_discovery_score <- colMeans(t(scale(t(expression_corrected[signature_genes, , drop = FALSE]))))
  signature_reference_rho <- suppressWarnings(stats::cor(signature_discovery_score, smooth_reference, method = "spearman"))
}
biological_direction_supported <- is.finite(original_reference_difference) && original_reference_difference > 0 &&
  is.finite(signature_reference_rho) && signature_reference_rho > 0
signature_ready <- signature_ready && biological_direction_supported
write.csv(data.frame(signature_genes = length(signature_genes),
  original_C2_reference_difference = original_reference_difference,
  signature_reference_Spearman_rho = signature_reference_rho,
  projection_allowed = signature_ready,
  interpretation = "positive_direction_check_not_proof_of_biological_specificity"),
  file.path(TABLE_DIR, "04_signature_biological_direction.csv"), row.names = FALSE)
if (!signature_ready) log_message("Projection skipped: insufficient signature or unsupported smooth-muscle direction. Labels are not flipped.")

#Module 13. C1 biology: external gene sets and within-C1 heterogeneity -----------
# These small curated panels are exploratory reference panels, not published
# classifiers, cell fractions, MCP-counter scores or EcoTyper assignments.
pathways <- reference_panels
hallmark_status <- "disabled"
hallmark_file <- HALLMARK_FILE
hallmark_table <- NULL
if (RUN_HALLMARK) {
  hallmark_result <- try({
    if (file.exists(hallmark_file)) {
      hallmark_table <- readRDS(hallmark_file)
      input_files <- c(input_files, hallmark_file)
    } else {
      if (!requireNamespace("msigdbr", quietly = TRUE) && INSTALL_MISSING_PACKAGES) {
        install.packages("msigdbr", repos = "https://cloud.r-project.org")
      }
      if (!requireNamespace("msigdbr", quietly = TRUE)) stop("msigdbr is not installed")
      if ("collection" %in% names(formals(msigdbr::msigdbr))) {
        hallmark_table <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
      } else {
        hallmark_table <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      }
    }
    stopifnot(all(c("gs_name", "gene_symbol") %in% names(hallmark_table)))
    hallmark_sets <- split(as.character(hallmark_table$gene_symbol), hallmark_table$gs_name)
    hallmark_sets <- lapply(hallmark_sets, unique)
    pathways <- c(pathways, hallmark_sets)
    saveRDS(hallmark_table, file.path(CACHE_DIR, "hallmark_sets_used.rds"))
    hallmark_status <- "loaded"
  }, silent = TRUE)
  if (inherits(hallmark_result, "try-error")) {
    hallmark_status <- paste("unavailable:", as.character(hallmark_result))
    log_message("Hallmark unavailable; curated reference panels continue.")
  }
}
writeLines(hallmark_status, file.path(LOG_DIR, "hallmark_status.txt"))
pathway_membership <- list()
for (pathway_name in names(pathways)) pathway_membership[[pathway_name]] <- data.frame(
  pathway = pathway_name, GeneSymbol = pathways[[pathway_name]],
  type = ifelse(grepl("^HALLMARK_", pathway_name), "MSigDB_Hallmark", "curated_exploratory_reference"))
write.csv(do.call(rbind, pathway_membership), file.path(TABLE_DIR, "10_reference_gene_sets.csv"), row.names = FALSE)
pathway_scores <- list()
pathway_coverage <- list()
for (pathway_name in names(pathways)) {
  available_genes <- intersect(pathways[[pathway_name]], rownames(expression_corrected))
  available_genes <- available_genes[matrixStats::rowVars(expression_corrected[available_genes, , drop = FALSE]) > 0]
  minimum_overlap <- if (grepl("^HALLMARK_", pathway_name)) 10L else 3L
  pathway_coverage[[pathway_name]] <- data.frame(pathway = pathway_name,
    reference_size = length(pathways[[pathway_name]]), available = length(available_genes),
    used = length(available_genes) >= minimum_overlap)
  if (length(available_genes) < minimum_overlap) next
  pathway_z <- t(scale(t(expression_corrected[available_genes, , drop = FALSE])))
  pathway_scores[[pathway_name]] <- data.frame(sampleid = discovery$sampleid,
    pathway = pathway_name, score = colMeans(pathway_z), C1_C2 = discovery$C1_C2,
    selected_cluster = discovery$selected_cluster, histology = discovery$histology, gse_id = discovery$gse_id)
}
pathway_score_table <- do.call(rbind, pathway_scores)
write.csv(pathway_score_table, file.path(TABLE_DIR, "10_discovery_reference_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, pathway_coverage), file.path(TABLE_DIR, "10_reference_gene_coverage.csv"), row.names = FALSE)
focused_pathways <- c(names(reference_panels), "HALLMARK_MYOGENESIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
focused_scores <- pathway_score_table[pathway_score_table$pathway %in% focused_pathways, , drop = FALSE]
focused_scores$pathway_display <- gsub("_", " ", sub("^HALLMARK_", "", focused_scores$pathway))
for (label_field in c("C1_C2", "selected_cluster")) {
  focused_scores$label <- focused_scores[[label_field]]
  biology_plot <- ggplot2::ggplot(focused_scores, ggplot2::aes(label, score, fill = label)) +
    ggplot2::geom_boxplot(outlier.size = 0.3) + ggplot2::facet_wrap(~ pathway_display, scales = "free_y", ncol = 3, labeller = ggplot2::label_wrap_gen(width = 25)) +
    ggplot2::theme_bw() + ggplot2::guides(fill = "none") +
    ggplot2::labs(x = label_field, y = "Mean gene z-score", title = "Reference panels: original groups and revised selected clusters",
      caption = "Exploratory expression scores; not cell proportions or formal molecular classes.")
  save_figure(biology_plot, paste0("10_reference_scores_", label_field), width = 10, height = 10)
}
c1_scores <- focused_scores[focused_scores$C1_C2 == "C1", , drop = FALSE]
for (group_field in c("histology", "gse_id")) {
  c1_scores$group <- as.character(c1_scores[[group_field]])
  c1_plot <- ggplot2::ggplot(c1_scores, ggplot2::aes(group, score, fill = group)) +
    ggplot2::geom_boxplot(outlier.size = 0.3) + ggplot2::facet_wrap(~ pathway_display, scales = "free_y", ncol = 3, labeller = ggplot2::label_wrap_gen(width = 25)) +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1)) +
    ggplot2::guides(fill = "none") + ggplot2::labs(x = group_field, y = "Mean gene z-score",
      title = "Heterogeneity within retained original C1")
  save_figure(c1_plot, paste0("10_C1_heterogeneity_", group_field), width = 11, height = 11)
  score_summary <- stats::aggregate(score ~ pathway + group, c1_scores,
    function(x) c(n = length(x), median = stats::median(x), q25 = stats::quantile(x, 0.25), q75 = stats::quantile(x, 0.75)))
  score_summary <- data.frame(score_summary[, c("pathway", "group")], score_summary$score)
  write.csv(score_summary, file.path(TABLE_DIR, paste0("10_C1_summary_", group_field, ".csv")), row.names = FALSE)
}
fgsea_status <- "not_run_no_Hallmark_or_DE"
if (hallmark_status == "loaded" && !is.null(de_tables$all)) {
  enrichment_result <- try({
    if (!requireNamespace("fgsea", quietly = TRUE) && INSTALL_MISSING_PACKAGES) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("fgsea", ask = FALSE, update = FALSE)
    }
    if (!requireNamespace("fgsea", quietly = TRUE)) stop("fgsea is not installed")
    ranked_statistics <- de_tables$all$t
    names(ranked_statistics) <- de_tables$all$GeneSymbol
    ranked_statistics <- sort(ranked_statistics[is.finite(ranked_statistics)], decreasing = TRUE)
    set.seed(SEED)
    enrichment <- as.data.frame(fgsea::fgseaMultilevel(pathways = hallmark_sets,
      stats = ranked_statistics, minSize = 10L, maxSize = 500L, eps = 0, nproc = 1L))
    enrichment$leadingEdge <- vapply(enrichment$leadingEdge, paste, character(1), collapse = ";")
    enrichment$direction <- ifelse(enrichment$NES > 0, "C2", "C1")
    write.csv(enrichment, file.path(TABLE_DIR, "10_Hallmark_fgsea_descriptive.csv"), row.names = FALSE)
    fgsea_status <- "completed_descriptive_cluster_derived_ranking"
  }, silent = TRUE)
  if (inherits(enrichment_result, "try-error")) fgsea_status <- paste("unavailable:", as.character(enrichment_result))
}
writeLines(fgsea_status, file.path(LOG_DIR, "fgsea_status.txt"))

#Module 14. Recompute TCGA projection from raw RSEM and revised signature --------
tcga_dir <- file.path(PROJECT_DIR, "Input", "TCGA", "sarc_tcga_pub")
tcga_expression_file <- file.path(tcga_dir, "data_mrna_seq_v2_rsem.txt")
tcga_sample_file <- file.path(tcga_dir, "data_clinical_sample.txt")
tcga_patient_file <- file.path(tcga_dir, "data_clinical_patient.txt")
tcga_files <- c(tcga_expression_file, tcga_sample_file, tcga_patient_file)
if (!all(file.exists(tcga_files))) stop("Missing local TCGA input: ", paste(tcga_files[!file.exists(tcga_files)], collapse = "; "))
input_files <- c(input_files, tcga_files)
tcga_expression_table <- data.table::fread(tcga_expression_file, data.table = FALSE, check.names = FALSE)
tcga_samples <- data.table::fread(tcga_sample_file, skip = "PATIENT_ID", data.table = FALSE,
  na.strings = c("", "NA", "N/A", "[Not Available]", "[Not Applicable]"))
tcga_patients <- data.table::fread(tcga_patient_file, skip = "PATIENT_ID", data.table = FALSE,
  na.strings = c("", "NA", "N/A", "[Not Available]", "[Not Applicable]"))
stopifnot(all(c("SAMPLE_ID", "PATIENT_ID") %in% names(tcga_samples)),
  "PATIENT_ID" %in% names(tcga_patients), !anyNA(tcga_samples$SAMPLE_ID),
  !anyNA(tcga_samples$PATIENT_ID), !anyNA(tcga_patients$PATIENT_ID),
  !anyDuplicated(tcga_samples$SAMPLE_ID), !anyDuplicated(tcga_patients$PATIENT_ID),
  "Hugo_Symbol" %in% names(tcga_expression_table))
tcga_metadata <- merge(tcga_samples, tcga_patients, by = "PATIENT_ID", all.x = TRUE,
  sort = FALSE, suffixes = c("_sample", "_patient"))
tcga_ids <- intersect(colnames(tcga_expression_table), tcga_metadata$SAMPLE_ID)
if (length(tcga_ids) < 30L) stop("Too few exactly matched TCGA sample IDs.")
tcga_metadata <- tcga_metadata[match(tcga_ids, tcga_metadata$SAMPLE_ID), , drop = FALSE]
if (anyDuplicated(tcga_metadata$PATIENT_ID)) {
  write.csv(tcga_metadata, file.path(TABLE_DIR, "11_multiple_samples_per_patient_review.csv"), row.names = FALSE)
  stop("TCGA has multiple specimens per patient; choose specimens explicitly before survival analysis.")
}
tcga_symbols <- trimws(as.character(tcga_expression_table$Hugo_Symbol))
tcga_expression <- as.matrix(tcga_expression_table[, tcga_ids, drop = FALSE])
storage.mode(tcga_expression) <- "double"
valid_symbol <- !is.na(tcga_symbols) & nzchar(tcga_symbols)
tcga_expression <- tcga_expression[valid_symbol, , drop = FALSE]
tcga_symbols <- tcga_symbols[valid_symbol]
if (anyDuplicated(tcga_symbols)) {
  duplicate_counts <- table(tcga_symbols)
  tcga_expression <- rowsum(tcga_expression, group = tcga_symbols, reorder = TRUE)
  tcga_expression <- tcga_expression / as.numeric(duplicate_counts[rownames(tcga_expression)])
} else rownames(tcga_expression) <- tcga_symbols
if (any(tcga_expression < 0, na.rm = TRUE)) stop("TCGA input is not nonnegative RSEM abundance.")
tcga_expression <- log2(tcga_expression + 1)
tcga_complete_genes <- rowSums(!is.finite(tcga_expression)) == 0
tcga_complete_genes <- tcga_complete_genes & matrixStats::rowVars(tcga_expression, na.rm = TRUE) > 0
tcga_complete_genes[is.na(tcga_complete_genes)] <- FALSE
write.csv(data.frame(GeneSymbol = rownames(tcga_expression), used = tcga_complete_genes),
  file.path(TABLE_DIR, "11_TCGA_gene_QC.csv"), row.names = FALSE)
tcga_expression <- tcga_expression[tcga_complete_genes, , drop = FALSE]
stopifnot(identical(colnames(tcga_expression), tcga_metadata$SAMPLE_ID))
required_tcga_fields <- c("SHORT_HISTOLOGY", "AGE", "SEX", "FNCLCC_GRADE",
  "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS", "PURITY", "LEUKOCYTE_FRACTION")
if (!all(required_tcga_fields %in% names(tcga_metadata))) {
  stop("Required TCGA fields absent: ", paste(setdiff(required_tcga_fields, names(tcga_metadata)), collapse = ", "))
}
tcga_metadata$histology_group <- toupper(trimws(tcga_metadata$SHORT_HISTOLOGY))
recognized_histology <- c("STLMS", "ULMS", "DDLPS", "UPS", "MFS", "SS", "MPNST")
if (any(!tcga_metadata$histology_group %in% recognized_histology)) {
  stop("Unrecognized TCGA SHORT_HISTOLOGY; review rather than guessing histology from site.")
}
tcga_metadata$age <- suppressWarnings(as.numeric(tcga_metadata$AGE))
tcga_metadata$age[!is.finite(tcga_metadata$age) | tcga_metadata$age < 0 | tcga_metadata$age > 120] <- NA_real_
sex_values <- tolower(trimws(as.character(tcga_metadata$SEX)))
sex_values[!sex_values %in% c("male", "female")] <- NA_character_
tcga_metadata$sex <- factor(sex_values, levels = c("female", "male"))
tcga_metadata$grade <- suppressWarnings(as.numeric(tcga_metadata$FNCLCC_GRADE))
tcga_metadata$grade[!tcga_metadata$grade %in% 1:3] <- NA_real_
tcga_metadata$histology <- factor(tcga_metadata$histology_group)
tcga_metadata$myogenic_score <- NA_real_
tcga_metadata$score_z <- NA_real_
tcga_metadata$projected_axis <- NA_character_
tcga_metadata$cor_C1 <- NA_real_
tcga_metadata$cor_C2 <- NA_real_
projection_genes <- intersect(signature_genes, rownames(tcga_expression))
projection_ready <- signature_ready && length(projection_genes) >= MIN_SIGNATURE_GENES
write.csv(data.frame(GeneSymbol = signature_genes, present_in_TCGA = signature_genes %in% projection_genes),
  file.path(TABLE_DIR, "11_TCGA_signature_coverage.csv"), row.names = FALSE)
if (projection_ready) {
  # Frozen discovery signature; TCGA is not reclustered, and survival is unused.
  signature_expression <- tcga_expression[projection_genes, , drop = FALSE]
  signature_z <- t(scale(t(signature_expression)))
  tcga_metadata$myogenic_score <- colMeans(signature_z)
  score_mean <- mean(tcga_metadata$myogenic_score)
  score_sd <- stats::sd(tcga_metadata$myogenic_score)
  if (!is.finite(score_sd) || score_sd <= 0) stop("Constant projected myogenic score.")
  tcga_metadata$score_z <- (tcga_metadata$myogenic_score - score_mean) / score_sd
  centroid_c1 <- rowMeans(expression_corrected[projection_genes, discovery$C1_C2 == "C1", drop = FALSE])
  centroid_c2 <- rowMeans(expression_corrected[projection_genes, discovery$C1_C2 == "C2", drop = FALSE])
  for (sample_index in seq_len(ncol(signature_expression))) {
    tcga_metadata$cor_C1[sample_index] <- stats::cor(signature_expression[, sample_index], centroid_c1)
    tcga_metadata$cor_C2[sample_index] <- stats::cor(signature_expression[, sample_index], centroid_c2)
  }
  tcga_metadata$projected_axis <- ifelse(tcga_metadata$cor_C2 >= tcga_metadata$cor_C1, "C2_like", "C1_like")
  tcga_metadata$projected_axis[!is.finite(tcga_metadata$cor_C1) | !is.finite(tcga_metadata$cor_C2)] <- NA_character_
  write.csv(data.frame(GeneSymbol = projection_genes, C1 = centroid_c1, C2 = centroid_c2,
    TCGA_gene_mean = rowMeans(signature_expression), TCGA_gene_sd = matrixStats::rowSds(signature_expression)),
    file.path(TABLE_DIR, "11_projection_centroids_and_scaling.csv"), row.names = FALSE)
  write.csv(data.frame(reference = "all_matched_TCGA_samples_before_scenario_subsetting",
    n = nrow(tcga_metadata), score_mean = score_mean, score_sd = score_sd),
    file.path(TABLE_DIR, "11_score_scaling_reference.csv"), row.names = FALSE)
  score_plot <- ggplot2::ggplot(tcga_metadata, ggplot2::aes(histology_group, score_z, fill = histology_group)) +
    ggplot2::geom_boxplot(outlier.size = 0.7) + ggplot2::theme_bw() + ggplot2::guides(fill = "none") +
    ggplot2::labs(x = "Histology", y = "Myogenic score (SD units)", title = "Fixed original-C2-derived signature in TCGA")
  save_figure(score_plot, "11_TCGA_myogenic_score")
}
tcga_metadata$projected_axis <- factor(tcga_metadata$projected_axis, levels = c("C1_like", "C2_like"))
scenario_histologies <- list(primary_174_definition = c("STLMS", "DDLPS", "UPS", "MFS", "SS"),
  plus_ULMS = c("STLMS", "ULMS", "DDLPS", "UPS", "MFS", "SS"),
  all_histologies = recognized_histology, LMS_only = c("STLMS", "ULMS"), STLMS_only = "STLMS")
auc_rows <- list()
if (projection_ready) {
  for (scenario in names(scenario_histologies)[1:3]) {
    current_data <- tcga_metadata[tcga_metadata$histology_group %in% scenario_histologies[[scenario]], , drop = FALSE]
    positive <- current_data$histology_group %in% c("STLMS", "ULMS")
    n_positive <- sum(positive)
    n_negative <- sum(!positive)
    auc <- (sum(rank(current_data$score_z)[positive]) - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)
    auc_rows[[scenario]] <- data.frame(scenario = scenario, n = nrow(current_data),
      n_LMS = n_positive, AUC = auc, role = "descriptive_histology_separation_not_prognosis")
  }
  write.csv(do.call(rbind, auc_rows), file.path(TABLE_DIR, "11_TCGA_histology_AUC.csv"), row.names = FALSE)
}

#Module 15. TCGA OS/DFS, Cox models, KM and PH diagnostics -----------------------
# Main reporting comparison: core histology set, OS, adjusted continuous score.
# All tests remain exploratory after inspecting v2.0.0; BH includes all fitted exposure tests.
# Adjustment: age + sex + ordinal grade + stratified baseline by histology.
# Nested models use identical complete cases within each scenario/endpoint.
# No outcome-optimized threshold. DSS is not inferred from OS follow-up time.
endpoint_definitions <- data.frame(endpoint = c("OS", "DFS"),
  time_field = c("OS_MONTHS", "DFS_MONTHS"), status_field = c("OS_STATUS", "DFS_STATUS"),
  event_value = c("1:DECEASED", "1:RECURRED/PROGRESSED"),
  censor_value = c("0:LIVING", "0:DISEASEFREE"), unit = "months")
write.csv(endpoint_definitions, file.path(TABLE_DIR, "12_TCGA_endpoint_definitions.csv"), row.names = FALSE)
survival_results <- list()
survival_counts <- list()
survival_status <- list()
ph_rows <- list()
exclusion_rows <- list()
set.seed(SEED)
if (projection_ready) {
  for (scenario in names(scenario_histologies)) {
    scenario_data <- tcga_metadata[tcga_metadata$histology_group %in% scenario_histologies[[scenario]], , drop = FALSE]
    for (endpoint_index in seq_len(nrow(endpoint_definitions))) {
      endpoint <- endpoint_definitions$endpoint[endpoint_index]
      analysis_data <- scenario_data
      analysis_data$time_months <- suppressWarnings(as.numeric(analysis_data[[endpoint_definitions$time_field[endpoint_index]]]))
      raw_status <- toupper(trimws(as.character(analysis_data[[endpoint_definitions$status_field[endpoint_index]]])))
      analysis_data$event <- NA_integer_
      analysis_data$event[!is.na(raw_status) & raw_status == endpoint_definitions$event_value[endpoint_index]] <- 1L
      analysis_data$event[!is.na(raw_status) & raw_status == endpoint_definitions$censor_value[endpoint_index]] <- 0L
      valid_endpoint <- is.finite(analysis_data$time_months) & analysis_data$time_months > 0 & !is.na(analysis_data$event)
      valid_covariates <- stats::complete.cases(analysis_data[, c("score_z", "projected_axis", "age", "sex", "grade", "histology")])
      excluded <- !valid_endpoint | !valid_covariates
      exclusion_rows[[paste(scenario, endpoint)]] <- data.frame(scenario = scenario, endpoint = endpoint,
        SAMPLE_ID = analysis_data$SAMPLE_ID, endpoint_eligible = valid_endpoint,
        covariates_complete = valid_covariates, included = !excluded,
        source_time = analysis_data$time_months, source_status = raw_status)
      analysis_data <- droplevels(analysis_data[!excluded, , drop = FALSE])
      survival_counts[[paste(scenario, endpoint)]] <- data.frame(scenario = scenario, endpoint = endpoint,
        source_n = nrow(scenario_data), endpoint_eligible_n = sum(valid_endpoint),
        model_n = nrow(analysis_data), events = sum(analysis_data$event), excluded_n = sum(excluded),
        C1_like_n = sum(analysis_data$projected_axis == "C1_like"),
        C2_like_n = sum(analysis_data$projected_axis == "C2_like"))
      if (nrow(analysis_data) < 10L || sum(analysis_data$event) < 5L) {
        survival_status[[length(survival_status) + 1L]] <- data.frame(scenario = scenario, endpoint = endpoint,
          model = "all", status = "insufficient_records_or_events", detail = "n < 10 or events < 5")
        next
      }
      analysis_data$sex <- droplevels(analysis_data$sex)
      analysis_data$histology <- droplevels(analysis_data$histology)
      adjustment_terms <- character()
      if (length(unique(analysis_data$age)) > 1) adjustment_terms <- c(adjustment_terms, "age")
      if (nlevels(analysis_data$sex) > 1) adjustment_terms <- c(adjustment_terms, "sex")
      if (length(unique(analysis_data$grade)) > 1) adjustment_terms <- c(adjustment_terms, "grade")
      if (nlevels(analysis_data$histology) > 1) adjustment_terms <- c(adjustment_terms, "strata(histology)")
      models <- list(continuous_unadjusted = Surv(time_months, event) ~ score_z,
        continuous_adjusted = stats::as.formula(paste("Surv(time_months, event) ~",
          paste(c("score_z", adjustment_terms), collapse = " + "))))
      if (nlevels(analysis_data$projected_axis) == 2L) {
        models$projected_axis_unadjusted <- Surv(time_months, event) ~ projected_axis
      }
      for (model_name in names(models)) {
        model_formula <- models[[model_name]]
        model_warnings <- character()
        fitted <- try(withCallingHandlers(survival::coxph(model_formula, data = analysis_data,
          ties = "efron", x = TRUE, y = TRUE, model = TRUE, singular.ok = FALSE),
          warning = function(w) { model_warnings <<- c(model_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), silent = TRUE)
        if (inherits(fitted, "try-error")) {
          survival_status[[length(survival_status) + 1L]] <- data.frame(scenario = scenario, endpoint = endpoint,
            model = model_name, status = "fit_failed", detail = as.character(fitted))
          next
        }
        parameter_count <- length(stats::coef(fitted))
        events_per_parameter <- sum(analysis_data$event) / parameter_count
        if (model_name == "continuous_adjusted" && events_per_parameter < MIN_EVENTS_PER_PARAMETER) {
          survival_status[[length(survival_status) + 1L]] <- data.frame(scenario = scenario, endpoint = endpoint,
            model = model_name, status = "not_reported_insufficient_events", detail = paste("events/parameter", round(events_per_parameter, 2)))
          next
        }
        fit_summary <- summary(fitted)
        coefficient_table <- as.data.frame(fit_summary$coefficients)
        confidence_table <- as.data.frame(fit_summary$conf.int)
        unstable <- any(!is.finite(as.matrix(confidence_table))) || length(model_warnings) > 0L
        model_status <- if (unstable) "review_unstable_or_warning" else "completed"
        result <- data.frame(scenario = scenario, endpoint = endpoint, model = model_name,
          term = rownames(coefficient_table), n = nrow(analysis_data), events = sum(analysis_data$event),
          parameters = parameter_count, events_per_parameter = events_per_parameter,
          HR = confidence_table[["exp(coef)"]], CI_low = confidence_table[["lower .95"]],
          CI_high = confidence_table[["upper .95"]], p_value = coefficient_table[["Pr(>|z|)"]],
          status = model_status, formula = paste(deparse(model_formula), collapse = " "))
        result$main_reporting_comparison <- scenario == "primary_174_definition" & endpoint == "OS" &
          model_name == "continuous_adjusted" & result$term == "score_z"
        survival_results[[length(survival_results) + 1L]] <- result
        survival_status[[length(survival_status) + 1L]] <- data.frame(scenario = scenario, endpoint = endpoint,
          model = model_name, status = model_status, detail = paste(model_warnings, collapse = " | "))
        model_stub <- paste(scenario, endpoint, model_name, sep = "_")
        saveRDS(fitted, file.path(CACHE_DIR, paste0("cox_", model_stub, ".rds")))
        ph_test <- try(survival::cox.zph(fitted), silent = TRUE)
        if (!inherits(ph_test, "try-error")) {
          ph_rows[[length(ph_rows) + 1L]] <- data.frame(scenario = scenario, endpoint = endpoint,
            model = model_name, term = rownames(ph_test$table), as.data.frame(ph_test$table), row.names = NULL)
          grDevices::pdf(file.path(PLOT_DIR, paste0("12_PH_", model_stub, ".pdf")), width = 7, height = 5)
          try(plot(ph_test), silent = TRUE)
          grDevices::dev.off()
        } else {
          writeLines(as.character(ph_test), file.path(LOG_DIR, paste0("PH_failed_", model_stub, ".txt")))
        }
      }
      if (nlevels(analysis_data$projected_axis) == 2L) {
        km_fit <- survival::survfit(Surv(time_months, event) ~ projected_axis, data = analysis_data)
        km_summary <- summary(km_fit, censored = TRUE)
        km_data <- data.frame(time = km_summary$time, survival = km_summary$surv,
          censored = km_summary$n.censor, group = as.character(km_summary$strata))
        initial_rows <- data.frame(time = 0, survival = 1, censored = 0,
          group = names(km_fit$strata))
        km_data <- rbind(initial_rows, km_data)
        km_data <- km_data[order(km_data$group, km_data$time), , drop = FALSE]
        logrank <- survival::survdiff(Surv(time_months, event) ~ projected_axis, data = analysis_data)
        logrank_p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1L, lower.tail = FALSE)
        km_plot <- ggplot2::ggplot(km_data, ggplot2::aes(time, survival, color = group)) +
          ggplot2::geom_step(linewidth = 0.7) +
          ggplot2::geom_point(data = km_data[km_data$censored > 0, , drop = FALSE], shape = 3, size = 1.2) +
          ggplot2::theme_bw() + ggplot2::coord_cartesian(ylim = c(0, 1)) +
          ggplot2::labs(x = "Months", y = paste(endpoint, "probability"), color = "Projected axis",
            title = paste(scenario, endpoint),
            subtitle = paste0("Original-target projection; n=", nrow(analysis_data), "; events=", sum(analysis_data$event),
              "; unadjusted log-rank P=", format.pval(logrank_p, digits = 3)),
            caption = "Original-target centroid projection; exploratory reanalysis. No optimized cutpoint or CI bands.")
        save_figure(km_plot, paste0("12_KM_", scenario, "_", endpoint))
        risk_times <- pretty(c(0, max(analysis_data$time_months)), n = 6)
        risk_times <- risk_times[risk_times >= 0 & risk_times <= max(analysis_data$time_months)]
        risk_summary <- summary(km_fit, times = risk_times, extend = TRUE)
        risk_table <- data.frame(time = risk_summary$time, group = as.character(risk_summary$strata),
          n_risk = risk_summary$n.risk, n_event = risk_summary$n.event, n_censor = risk_summary$n.censor)
        write.csv(risk_table, file.path(TABLE_DIR, paste0("12_KM_risk_", scenario, "_", endpoint, ".csv")), row.names = FALSE)
        write.csv(km_data, file.path(TABLE_DIR, paste0("12_KM_coordinates_", scenario, "_", endpoint, ".csv")), row.names = FALSE)
      }
    }
  }
}
if (length(survival_results)) {
  survival_table <- do.call(rbind, survival_results)
  survival_table$exposure_test <-
    survival_table$term %in% c("score_z", "projected_axisC2_like")
  survival_table$all_exposure_family_BH <- NA_real_
  secondary_rows <- which(survival_table$exposure_test & survival_table$status == "completed")
  survival_table$all_exposure_family_BH[secondary_rows] <- stats::p.adjust(survival_table$p_value[secondary_rows], "BH")
  write.csv(survival_table, file.path(TABLE_DIR, "12_TCGA_Cox_all_results.csv"), row.names = FALSE)
  forest_data <- survival_table[survival_table$term == "score_z" & survival_table$status == "completed" &
    is.finite(survival_table$CI_low) & is.finite(survival_table$CI_high), , drop = FALSE]
  if (nrow(forest_data)) {
    forest_data$label <- paste(forest_data$scenario, forest_data$model, sep = " / ")
    forest_plot <- ggplot2::ggplot(forest_data, ggplot2::aes(HR, label)) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2) +
      ggplot2::geom_segment(ggplot2::aes(x = CI_low, xend = CI_high, yend = label)) +
      ggplot2::geom_point() + ggplot2::scale_x_log10() + ggplot2::facet_wrap(~ endpoint) +
      ggplot2::theme_bw() + ggplot2::labs(x = "HR per 1 SD myogenic score (95% CI)", y = NULL)
    save_figure(forest_plot, "12_TCGA_Cox_forest", width = 12, height = 7)
  }
}
if (length(survival_counts)) write.csv(do.call(rbind, survival_counts), file.path(TABLE_DIR, "12_TCGA_survival_counts.csv"), row.names = FALSE)
if (length(exclusion_rows)) write.csv(do.call(rbind, exclusion_rows), file.path(TABLE_DIR, "12_TCGA_survival_inclusion_audit.csv"), row.names = FALSE)
if (length(survival_status)) write.csv(do.call(rbind, survival_status), file.path(TABLE_DIR, "12_TCGA_model_status.csv"), row.names = FALSE)
if (length(ph_rows)) write.csv(do.call(rbind, ph_rows), file.path(TABLE_DIR, "12_TCGA_PH_tests.csv"), row.names = FALSE)
write.csv(tcga_metadata, file.path(TABLE_DIR, "12_TCGA_per_sample_revised.csv"), row.names = FALSE)


#Module 16. TCGA microenvironment context and optional published class labels ---
tcga_context_rows <- list()
tcga_reference_rows <- list()
for (pathway_name in intersect(focused_pathways, names(pathways))) {
  genes <- intersect(pathways[[pathway_name]], rownames(tcga_expression))
  if (length(genes) < ifelse(grepl("^HALLMARK_", pathway_name), 10L, 3L)) next
  panel_score <- colMeans(t(scale(t(tcga_expression[genes, , drop = FALSE]))))
  context_data <- data.frame(SAMPLE_ID = tcga_metadata$SAMPLE_ID, pathway = pathway_name,
    score = panel_score, histology = tcga_metadata$histology,
    PURITY = suppressWarnings(as.numeric(tcga_metadata$PURITY)),
    LEUKOCYTE_FRACTION = suppressWarnings(as.numeric(tcga_metadata$LEUKOCYTE_FRACTION)))
  tcga_reference_rows[[pathway_name]] <- context_data
  for (context_variable in c("PURITY", "LEUKOCYTE_FRACTION")) {
    context_data$value <- context_data[[context_variable]]
    observed <- droplevels(context_data[is.finite(context_data$value) & is.finite(context_data$score), , drop = FALSE])
    if (nrow(observed) < 10L || length(unique(observed$value)) < 2L) next
    context_test <- stats::cor.test(observed$score, observed$value, method = "spearman", exact = FALSE)
    context_model <- if (nlevels(observed$histology) > 1L)
      stats::lm(score ~ value + histology, observed) else stats::lm(score ~ value, observed)
    model_coefficients <- summary(context_model)$coefficients
    tcga_context_rows[[length(tcga_context_rows) + 1L]] <- data.frame(pathway = pathway_name,
      covariate = context_variable, n = nrow(observed), Spearman_rho = as.numeric(context_test$estimate),
      p_value = context_test$p.value, histology_adjusted_beta = model_coefficients["value", "Estimate"],
      histology_adjusted_p = model_coefficients["value", "Pr(>|t|)"])
    context_plot <- ggplot2::ggplot(observed, ggplot2::aes(value, score, color = histology)) +
      ggplot2::geom_point(alpha = 0.65, size = 1.3) + ggplot2::theme_bw() +
      ggplot2::labs(x = context_variable, y = "Mean gene z-score", title = pathway_name)
    save_figure(context_plot, paste0("13_", pathway_name, "_", context_variable))
  }
}
if (length(tcga_reference_rows)) write.csv(do.call(rbind, tcga_reference_rows),
  file.path(TABLE_DIR, "13_TCGA_reference_scores.csv"), row.names = FALSE)
if (length(tcga_context_rows)) {
  context_table <- do.call(rbind, tcga_context_rows)
  context_table$Spearman_BH <- stats::p.adjust(context_table$p_value, "BH")
  context_table$adjusted_model_BH <- stats::p.adjust(context_table$histology_adjusted_p, "BH")
  write.csv(context_table, file.path(TABLE_DIR, "13_TCGA_context_associations.csv"), row.names = FALSE)
}
if (RUN_HALLMARK && (hallmark_status != "loaded" || !startsWith(fgsea_status, "completed"))) {
  stop("Primary Hallmark enrichment incomplete: ", hallmark_status, "; ", fgsea_status)
}
V21_DIR <- OUTPUT_DIR

#Module 17. 03_Matched_sensitivity output paths ----
SEED <- 20261026L; set.seed(SEED)
OUTPUT_DIR <- file.path(PRIMARY_RESULTS, "03_Matched_sensitivity")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables"); PLOT_DIR <- file.path(OUTPUT_DIR, "figures")
CACHE_DIR <- file.path(OUTPUT_DIR, "cache"); LOG_DIR <- file.path(OUTPUT_DIR, "session")
for (directory in c(TABLE_DIR, PLOT_DIR, CACHE_DIR, LOG_DIR)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
input_files <- character()
#Module 18. Parent inputs and fixed original target ---------------------------
parent_files <- c(file.path(V20_DIR, "tables", "03_record_lineage_and_selection.csv"),
  file.path(V21_DIR, "tables", "02_original_target_metadata.csv"),
  file.path(V20_DIR, "cache", "discovery_deduplicated_uncorrected.rds"),
  file.path(V20_DIR, "cache", "discovery_deduplicated_corrected.rds"),
  file.path(V21_DIR, "tables", "12_TCGA_per_sample_revised.csv"),
  file.path(V21_DIR, "tables", "12_TCGA_Cox_all_results.csv"),
  file.path(V21_DIR, "tables", "12_TCGA_survival_inclusion_audit.csv"),
  file.path(V21_DIR, "tables", "10_discovery_reference_scores.csv"),
  file.path(V21_DIR, "tables", "13_TCGA_reference_scores.csv"),
  file.path(V21_DIR, "tables", "09_revised_discovery_signature.csv"))
if (!all(file.exists(parent_files))) stop("Missing input: ", paste(parent_files[!file.exists(parent_files)], collapse = "; "))
input_files <- parent_files
all_metadata <- data.table::fread(parent_files[1], data.table = FALSE)
discovery <- data.table::fread(parent_files[2], data.table = FALSE)
expression_raw <- readRDS(parent_files[3]); expression_corrected <- readRDS(parent_files[4])
tcga <- data.table::fread(parent_files[5], data.table = FALSE)
parent_cox <- data.table::fread(parent_files[6], data.table = FALSE)
parent_inclusion <- data.table::fread(parent_files[7], data.table = FALSE)
discovery_reference <- data.table::fread(parent_files[8], data.table = FALSE)
tcga_reference <- data.table::fread(parent_files[9], data.table = FALSE)
signature_genes <- data.table::fread(parent_files[10])$GeneSymbol
stopifnot(!anyNA(discovery$sampleid), !anyDuplicated(discovery$sampleid),
  !anyNA(discovery$record_group), !anyDuplicated(discovery$record_group),
  !anyNA(all_metadata$sampleid), !anyDuplicated(all_metadata$sampleid),
  !anyNA(all_metadata$record_group), all(discovery$sampleid %in% all_metadata$sampleid),
  setequal(discovery$record_group, all_metadata$record_group),
  identical(as.character(discovery$record_group), as.character(all_metadata$record_group[match(discovery$sampleid, all_metadata$sampleid)])),
  all(discovery$C1_C2 == discovery$legacy_cluster), all(discovery$C1_C2 %in% c("C1", "C2")),
  identical(dimnames(expression_raw), dimnames(expression_corrected)),
  !anyDuplicated(rownames(expression_raw)), !anyDuplicated(colnames(expression_raw)),
  setequal(discovery$sampleid, colnames(expression_raw)),
  all(is.finite(expression_raw)), all(is.finite(expression_corrected)),
  !anyDuplicated(tcga$SAMPLE_ID), !anyDuplicated(tcga$PATIENT_ID), !anyNA(tcga$score_z))
expression_raw <- expression_raw[, discovery$sampleid, drop = FALSE]
expression_corrected <- expression_corrected[, discovery$sampleid, drop = FALSE]
discovery$C1_C2 <- factor(discovery$legacy_cluster, levels = c("C1", "C2"))
stopifnot(all(table(discovery$C1_C2) > 0L))
discovery$histology <- factor(discovery$histology); discovery$batch <- factor(discovery$gse_id)
all_metadata$priority <- as.numeric(all_metadata$priority)
for (record_id in discovery$record_group) {
  old_labels <- unique(na.omit(all_metadata$legacy_cluster[all_metadata$record_group == record_id]))
  if (length(old_labels) != 1L || old_labels != as.character(discovery$C1_C2[discovery$record_group == record_id])) stop("Linked original labels disagree.")
}
write.csv(data.frame(item = c("original_records", "retained_record_groups", "original_C1", "original_C2", "TCGA_patients"),
  n = c(nrow(all_metadata), nrow(discovery), table(discovery$C1_C2), nrow(tcga))),
  file.path(TABLE_DIR, "02_fixed_sample_counts.csv"), row.names = FALSE)
log_message("Validated", nrow(all_metadata), "records;", nrow(discovery), "record groups;", nrow(tcga), "TCGA patients.")

#Module 19. Rebuild all-record matrix on one frozen gene universe --------------
# The v2.0 retained genes are fixed across conditions; no new variance filter.
cluster_status <- "disabled_by_setting"
if (RUN_CLUSTERING) {
  frozen_genes <- rownames(expression_raw); expression_list <- list()
  for (series_id in unique(all_metadata$gse_id)) {
    expression_file <- file.path(PROJECT_DIR, "Input", "GEO_prepared", series_id,
      paste0(series_id, "_ExpressionMatrix_Annotated.csv"))
    if (!file.exists(expression_file)) stop("Required prepared matrix missing: ", expression_file)
    input_files <- c(input_files, expression_file)
    source_table <- data.table::fread(expression_file, data.table = FALSE, check.names = FALSE)
    gene_column <- if ("Gene_Symbol" %in% names(source_table)) "Gene_Symbol" else names(source_table)[1]
    source_symbols <- trimws(as.character(source_table[[gene_column]]))
    series_metadata <- all_metadata[all_metadata$gse_id == series_id, , drop = FALSE]
    stopifnot(!anyNA(source_symbols), !anyDuplicated(source_symbols), all(nzchar(source_symbols)),
      all(frozen_genes %in% source_symbols), all(series_metadata$expr_col_original %in% names(source_table)))
    series_expression <- as.matrix(source_table[match(frozen_genes, source_symbols), series_metadata$expr_col_original, drop = FALSE])
    storage.mode(series_expression) <- "double"
    dimnames(series_expression) <- list(frozen_genes, series_metadata$sampleid)
    if (any(!is.finite(series_expression))) stop("Nonfinite input expression: ", series_id)
    expression_list[[series_id]] <- series_expression
    log_message("Loaded prepared matrix:", series_id, ncol(series_expression), "records")
  }
  expression_all_raw <- do.call(cbind, expression_list)[, all_metadata$sampleid, drop = FALSE]
  rm(expression_list, source_table, series_expression)
  raw_difference <- max(abs(expression_all_raw[, discovery$sampleid, drop = FALSE] - expression_raw))
  write.csv(data.frame(max_abs_difference = raw_difference, tolerance = 1e-7),
    file.path(TABLE_DIR, "03_prepared_vs_parent_matrix_check.csv"), row.names = FALSE)
  if (!is.finite(raw_difference) || raw_difference > 1e-7) stop("Prepared inputs changed since v2.0.0; resolve provenance before comparing methods.")
  all_batch <- factor(all_metadata$gse_id)
  expression_all_corrected <- if (nlevels(all_batch) > 1L) sva::ComBat(expression_all_raw,
    batch = all_batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE, mean.only = any(table(all_batch) == 1L)) else expression_all_raw
  stopifnot(all(is.finite(expression_all_corrected)), identical(dimnames(expression_all_raw), dimnames(expression_all_corrected)))
  correction_difference <- max(abs(expression_all_corrected[, discovery$sampleid, drop = FALSE] - expression_corrected))
  write.csv(data.frame(genes = length(frozen_genes), all_records = nrow(all_metadata),
    representatives = nrow(discovery), max_abs_refit_vs_fixed_correction = correction_difference),
    file.path(TABLE_DIR, "03_factorial_matrix_dimensions.csv"), row.names = FALSE)
  if (SAVE_LARGE_MATRICES) saveRDS(expression_all_corrected, file.path(CACHE_DIR, "all_records_corrected_frozen_genes.rds"))
}

#Module 20. Six matched clustering conditions and group-wise stability --------
# Selection and silhouette use the same space within each condition.
# Stability conditions on the fitted correction/embedding; it does NOT measure
# full-pipeline refitting stability, nor reproduce the original 30-resample ARI.
if (RUN_CLUSTERING) {
  condition_grid <- expand.grid(matrix_state = c("all_records", "dedup_fixed_correction", "dedup_refit_correction"),
    space = c("full_gene", "top_variable_PC"), stringsAsFactors = FALSE)
  condition_grid$condition <- paste(condition_grid$matrix_state, condition_grid$space, sep = "__")
  cluster_rows <- list(); support_rows <- list(); stability_rows <- list(); recovery_rows <- list()
  group_ids <- sort(unique(all_metadata$record_group)); group_draws <- vector("list", SUBSAMPLE_REPLICATES)
  stopifnot(length(group_ids) > max(K_GRID) + 2L)
  set.seed(SEED)
  for (b in seq_len(SUBSAMPLE_REPLICATES)) group_draws[[b]] <- sample(group_ids,
    max(max(K_GRID) + 2L, floor(length(group_ids) * SUBSAMPLE_FRACTION)), replace = FALSE)
  for (condition_index in seq_len(nrow(condition_grid))) {
    condition <- condition_grid$condition[condition_index]; state <- condition_grid$matrix_state[condition_index]
    current_metadata <- if (state == "all_records") all_metadata else discovery
    current_expression <- if (state == "all_records") expression_all_corrected else if (state == "dedup_fixed_correction")
      expression_all_corrected[, discovery$sampleid, drop = FALSE] else expression_corrected
    if (condition_grid$space[condition_index] == "full_gene") {
      current_coordinates <- t(current_expression)
    } else {
      variable_genes <- head(order(matrixStats::rowVars(current_expression), decreasing = TRUE), TOP_VARIABLE_GENES)
      current_pca <- stats::prcomp(t(current_expression[variable_genes, , drop = FALSE]), center = TRUE,
        scale. = FALSE, rank. = min(NUMBER_OF_PCS, ncol(current_expression) - 1L, length(variable_genes)))
      current_coordinates <- current_pca$x
      write.csv(data.frame(GeneSymbol = rownames(current_expression)[variable_genes]),
        file.path(TABLE_DIR, paste0("04_features_", condition, ".csv")), row.names = FALSE)
    }
    # Euclidean distances via a Gram matrix; avoids repeated high-dimensional dist().
    coordinates_centered <- sweep(current_coordinates, 2L, colMeans(current_coordinates), "-")
    gram_matrix <- tcrossprod(coordinates_centered); norms <- diag(gram_matrix)
    squared_distance <- outer(norms, norms, "+") - 2 * gram_matrix
    squared_distance[squared_distance < 0] <- 0; distance_matrix <- sqrt(squared_distance); diag(distance_matrix) <- 0
    dimnames(distance_matrix) <- list(current_metadata$sampleid, current_metadata$sampleid)
    current_distance <- stats::as.dist(distance_matrix); current_tree <- stats::hclust(current_distance, method = "ward.D2")
    full_labels <- stats::cutree(current_tree, k = K_GRID)
    representative_index <- match(discovery$sampleid, current_metadata$sampleid)
    stopifnot(!anyNA(representative_index))
    current_support <- data.frame(condition = condition, matrix_state = state, space = condition_grid$space[condition_index],
      k = K_GRID, silhouette = NA_real_, mean_conditional_ARI = NA_real_, ARI_q025 = NA_real_, ARI_q975 = NA_real_)
    current_ari <- matrix(NA_real_, nrow = SUBSAMPLE_REPLICATES, ncol = length(K_GRID))
    for (ki in seq_along(K_GRID)) {
      labels <- full_labels[, ki]
      current_support$silhouette[ki] <- mean(cluster::silhouette(labels, current_distance)[, "sil_width"])
      cluster_rows[[length(cluster_rows) + 1L]] <- data.frame(condition = condition, sampleid = discovery$sampleid,
        record_group = discovery$record_group, k = K_GRID[ki], numeric_cluster = labels[representative_index],
        original_group = discovery$C1_C2, histology = discovery$histology)
      ref_labels <- labels[representative_index]; target <- discovery$C1_C2 == "C2"
      for (cl in sort(unique(ref_labels))) {
        member <- ref_labels == cl; overlap <- sum(member & target)
        recovery_rows[[length(recovery_rows) + 1L]] <- data.frame(condition = condition, k = K_GRID[ki],
          numeric_cluster = cl, n = sum(member), original_C2_overlap = overlap,
          target_recall = overlap / sum(target), target_precision = overlap / sum(member),
          target_Jaccard = overlap / sum(member | target),
          LMS_n = sum(member & as.character(discovery$histology) == "Leiomyosarcoma"))
      }
    }
    for (b in seq_len(SUBSAMPLE_REPLICATES)) {
      chosen <- which(current_metadata$record_group %in% group_draws[[b]])
      chosen_rep <- which(discovery$record_group %in% group_draws[[b]])
      sub_tree <- stats::hclust(stats::as.dist(distance_matrix[chosen, chosen, drop = FALSE]), method = "ward.D2")
      sub_labels <- stats::cutree(sub_tree, k = K_GRID)
      sub_rep <- match(discovery$sampleid[chosen_rep], current_metadata$sampleid[chosen])
      stopifnot(!anyNA(sub_rep))
      for (ki in seq_along(K_GRID)) current_ari[b, ki] <- mclust::adjustedRandIndex(
        full_labels[representative_index[chosen_rep], ki], sub_labels[sub_rep, ki])
    }
    current_support$mean_conditional_ARI <- colMeans(current_ari)
    current_support$ARI_q025 <- apply(current_ari, 2L, stats::quantile, probs = 0.025)
    current_support$ARI_q975 <- apply(current_ari, 2L, stats::quantile, probs = 0.975)
    current_support$selected_by_silhouette <- current_support$k == current_support$k[order(-current_support$silhouette, current_support$k)][1]
    support_rows[[condition]] <- current_support
    stability_rows[[condition]] <- data.frame(condition = condition, replicate = rep(seq_len(SUBSAMPLE_REPLICATES), times = length(K_GRID)),
      k = rep(K_GRID, each = SUBSAMPLE_REPLICATES), conditional_ARI = as.vector(current_ari))
    saveRDS(current_tree, file.path(CACHE_DIR, paste0("tree_", condition, ".rds")))
    log_message("Clustering complete:", condition)
  }
  cluster_table <- do.call(rbind, cluster_rows); support_table <- do.call(rbind, support_rows)
  recovery_table <- do.call(rbind, recovery_rows); stability_table <- do.call(rbind, stability_rows)
  comparison_rows <- list()
  for (k in K_GRID) {
    lab <- reshape(cluster_table[cluster_table$k == k, c("sampleid", "condition", "numeric_cluster")],
      idvar = "sampleid", timevar = "condition", direction = "wide")
    lab <- lab[match(discovery$sampleid, lab$sampleid), , drop = FALSE]
    for (i in seq_len(nrow(condition_grid))) {
      v1 <- lab[[paste0("numeric_cluster.", condition_grid$condition[i])]]
      comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(k = k, condition_a = condition_grid$condition[i],
        condition_b = "original_C1_C2", ARI = mclust::adjustedRandIndex(v1, discovery$C1_C2), n = nrow(discovery))
      if (i < nrow(condition_grid)) for (j in seq.int(i + 1L, nrow(condition_grid))) {
        v2 <- lab[[paste0("numeric_cluster.", condition_grid$condition[j])]]
        comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(k = k, condition_a = condition_grid$condition[i],
          condition_b = condition_grid$condition[j], ARI = mclust::adjustedRandIndex(v1, v2), n = nrow(discovery))
      }
    }
  }
  write.csv(support_table, file.path(TABLE_DIR, "04_clustering_support.csv"), row.names = FALSE)
  write.csv(cluster_table, file.path(TABLE_DIR, "04_cluster_labels_on_common_records.csv"), row.names = FALSE)
  write.csv(stability_table, file.path(TABLE_DIR, "04_conditional_subsampling.csv"), row.names = FALSE)
  write.csv(recovery_table, file.path(TABLE_DIR, "04_original_C2_recovery_all_clusters.csv"), row.names = FALSE)
  write.csv(do.call(rbind, comparison_rows), file.path(TABLE_DIR, "04_partition_comparisons.csv"), row.names = FALSE)
  for (measure in c("silhouette", "mean_conditional_ARI")) {
    support_table$value <- support_table[[measure]]
    p <- ggplot2::ggplot(support_table, ggplot2::aes(k, value, color = space, group = space)) +
      ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::facet_wrap(~ matrix_state, nrow = 1) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::scale_x_continuous(breaks = K_GRID) +
      ggplot2::labs(y = measure, color = "Expression space", title = "Clustering sensitivity on a fixed gene universe",
        caption = "ARI: fixed correction/embedding, record-group subsampling, one representative per group for evaluation.\nThese ARIs are not directly comparable with earlier full-pipeline stability values.")
    save_figure(p, paste0("04_", measure), width = 12, height = 4.8)
  }
  cluster_status <- "completed_six_conditions_conditional_stability"
  rm(current_expression, current_coordinates, coordinates_centered, gram_matrix, squared_distance, distance_matrix)
}

#Module 21. Within-source biology and original-C1 heterogeneity ----------------
# Fixed labels were derived from expression. DE is descriptive, not an
# independent validation test. Source restriction cannot remove label circularity.
within_source_rows <- list(); within_source_status <- list()
for (subset_name in c("GSE71118_all", "GSE71118_LMS")) {
  keep <- discovery$gse_id == "GSE71118"
  if (subset_name == "GSE71118_LMS") keep <- keep & discovery$histology == "Leiomyosarcoma"
  d <- droplevels(discovery[keep, , drop = FALSE])
  counts <- table(factor(d$C1_C2, levels = c("C1", "C2")))
  if (any(counts < 3L)) {
    within_source_status[[subset_name]] <- data.frame(subset = subset_name, status = "insufficient_group_size", detail = paste(counts, collapse = "/"))
    next
  }
  design <- stats::model.matrix(~ C1_C2, d)
  fit <- limma::eBayes(limma::lmFit(expression_raw[, d$sampleid, drop = FALSE], design), trend = TRUE)
  de <- limma::topTable(fit, coef = "C1_C2C2", number = Inf, sort.by = "none", adjust.method = "BH")
  de$GeneSymbol <- rownames(de); de$subset <- subset_name
  de$in_fixed_signature <- de$GeneSymbol %in% signature_genes
  de$n_C1 <- counts[1]; de$n_C2 <- counts[2]
  within_source_rows[[subset_name]] <- de
  within_source_status[[subset_name]] <- data.frame(subset = subset_name, status = "completed_descriptive", detail = "Raw expression; one source; original fixed labels")
}
if (length(within_source_rows)) {
  within_de <- do.call(rbind, within_source_rows)
  write.csv(within_de, file.path(TABLE_DIR, "05_within_source_DE_all_genes.csv"), row.names = FALSE)
  write.csv(within_de[within_de$in_fixed_signature, ], file.path(TABLE_DIR, "05_within_source_fixed_signature.csv"), row.names = FALSE)
}
write.csv(do.call(rbind, within_source_status), file.path(TABLE_DIR, "05_within_source_status.csv"), row.names = FALSE)
reference_panels <- list(
  Smooth_muscle_reference = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "LMOD1", "MYOCD", "DES"),
  Fibroblast_ECM_reference = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A1", "COL6A2"),
  T_cell_reference = c("CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK", "CD2"),
  Cytotoxic_reference = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "CTSW"),
  Myeloid_reference = c("LST1", "TYROBP", "FCER1G", "AIF1", "CSF1R", "CD68", "CTSS"),
  Endothelial_reference = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "ESAM", "CDH5"))
reference_summary <- list(); correlation_rows <- list(); panel_coverage <- list()
z_expression <- sweep(expression_corrected, 1L, rowMeans(expression_corrected), "-")
gene_sd <- sqrt(matrixStats::rowVars(expression_corrected)); gene_sd[gene_sd == 0] <- NA_real_
z_expression <- sweep(z_expression, 1L, gene_sd, "/")
fixed_score <- colMeans(z_expression[intersect(signature_genes, rownames(z_expression)), , drop = FALSE], na.rm = TRUE)
for (pathway in names(reference_panels)) {
  genes <- intersect(reference_panels[[pathway]], rownames(z_expression))
  independent_genes <- setdiff(genes, signature_genes)
  panel_coverage[[pathway]] <- data.frame(pathway = pathway, panel_size = length(reference_panels[[pathway]]),
    observed_n = length(genes), signature_overlap_n = sum(genes %in% signature_genes),
    nonoverlap_n = length(independent_genes), observed_genes = paste(genes, collapse = ";"),
    nonoverlap_genes = paste(independent_genes, collapse = ";"))
  for (scope in c("all_records", "within_LMS", "within_GSE71118")) {
    keep <- if (scope == "all_records") rep(TRUE, nrow(discovery)) else if (scope == "within_LMS")
      discovery$histology == "Leiomyosarcoma" else discovery$gse_id == "GSE71118"
    for (version in c("full_panel", "signature_overlap_removed")) {
      current_genes <- if (version == "full_panel") genes else independent_genes
      if (length(current_genes) < 3L) next
      panel_score <- colMeans(z_expression[current_genes, , drop = FALSE], na.rm = TRUE)
      ok <- keep & is.finite(fixed_score) & is.finite(panel_score)
      if (sum(ok) < 10L || stats::sd(panel_score[ok]) == 0 || stats::sd(fixed_score[ok]) == 0) next
      test <- stats::cor.test(fixed_score[ok], panel_score[ok], method = "spearman", exact = FALSE)
      correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(pathway = pathway, scope = scope,
        version = version, n = sum(ok), genes_n = length(current_genes), rho = unname(test$estimate), p_value = test$p.value)
    }
  }
}
write.csv(do.call(rbind, panel_coverage), file.path(TABLE_DIR, "05_reference_signature_overlap.csv"), row.names = FALSE)
if (length(correlation_rows)) {
  correlations <- do.call(rbind, correlation_rows)
  correlations$BH_within_module <- p.adjust(correlations$p_value, "BH")
  write.csv(correlations, file.path(TABLE_DIR, "05_reference_correlations_overlap_sensitivity.csv"), row.names = FALSE)
}
# Existing Hallmark/reference scores are summarized; no ecotype labels are inferred.
c1_reference <- discovery_reference[discovery_reference$C1_C2 == "C1", , drop = FALSE]
for (group_variable in c("histology", "gse_id")) {
  for (key in unique(c1_reference[[group_variable]])) {
    for (pathway in unique(c1_reference$pathway)) {
      x <- c1_reference$score[c1_reference[[group_variable]] == key & c1_reference$pathway == pathway]
      x <- x[is.finite(x)]; if (!length(x)) next
      reference_summary[[length(reference_summary) + 1L]] <- data.frame(group_variable = group_variable,
        group = key, pathway = pathway, n = length(x), median = median(x),
        q25 = unname(quantile(x, .25)), q75 = unname(quantile(x, .75)))
    }
  }
}
c1_summary <- do.call(rbind, reference_summary)
write.csv(c1_summary, file.path(TABLE_DIR, "05_C1_reference_by_histology_and_source.csv"), row.names = FALSE)
for (group_variable in c("histology", "gse_id")) {
  plotted <- c1_summary[c1_summary$group_variable == group_variable &
    c1_summary$pathway %in% names(reference_panels), , drop = FALSE]
  if (!nrow(plotted)) next
  plotted$pathway <- sub("_reference$", "", plotted$pathway)
  p <- ggplot2::ggplot(plotted, ggplot2::aes(pathway, group, fill = median)) +
    ggplot2::geom_tile(color = "white") + ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n)), size = 2.6) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Median score", title = "Original C1: descriptive expression programs",
      caption = "Reference expression scores are not measured cell fractions or validated ecotype assignments.")
  save_figure(p, paste0("05_C1_reference_", group_variable), 10, 6)
}
rm(z_expression)
log_message("Completed within-source fixed-target DE and C1 descriptive context.")

#Module 22. Metadata completeness and linked-record disagreements -------------
# Report source values; do not resolve contradictory clinical records by majority.
metadata_fields <- intersect(c("age", "sex", "gender", "grade_h", "stage_h", "treatment_h", "time", "metastasis"), names(all_metadata))
missing_tokens <- c("", "na", "n/a", "nan", "unknown", "[unknown]", "not available", "not reported")
completeness_rows <- list(); consensus_rows <- list(); conflict_rows <- list()
for (field in metadata_fields) {
  values <- trimws(as.character(all_metadata[[field]]))
  values[is.na(values) | tolower(values) %in% missing_tokens] <- NA_character_
  if (field %in% c("sex", "gender", "metastasis")) values <- tolower(values)
  for (series_id in unique(all_metadata$gse_id)) {
    keep <- all_metadata$gse_id == series_id
    completeness_rows[[length(completeness_rows) + 1L]] <- data.frame(gse_id = series_id, field = field,
      n = sum(keep), nonmissing_n = sum(!is.na(values[keep])), missing_n = sum(is.na(values[keep])))
  }
  for (record_id in discovery$record_group) {
    ii <- which(all_metadata$record_group == record_id); observed <- unique(na.omit(values[ii]))
    consensus_rows[[length(consensus_rows) + 1L]] <- data.frame(record_group = record_id, field = field,
      consensus_value = if (length(observed) == 1L) observed else NA_character_,
      status = if (!length(observed)) "missing" else if (length(observed) == 1L) "agreement_or_single_value" else "conflict_unresolved",
      source_values = paste(paste(all_metadata$gse_id[ii], values[ii], sep = ":"), collapse = " | "))
    if (length(observed) > 1L) conflict_rows[[length(conflict_rows) + 1L]] <- data.frame(
      record_group = record_id, field = field, sampleid = all_metadata$sampleid[ii],
      gse_id = all_metadata$gse_id[ii], value = values[ii])
  }
}
write.csv(do.call(rbind, completeness_rows), file.path(TABLE_DIR, "06_metadata_completeness_by_source.csv"), row.names = FALSE)
write.csv(do.call(rbind, consensus_rows), file.path(TABLE_DIR, "06_linked_metadata_consensus_audit.csv"), row.names = FALSE)
if (length(conflict_rows)) write.csv(do.call(rbind, conflict_rows), file.path(TABLE_DIR, "06_linked_metadata_conflicts.csv"), row.names = FALSE)
confounding_design <- model.matrix(~ batch + histology, discovery)
write.csv(data.frame(n = nrow(discovery), columns = ncol(confounding_design), rank = qr(confounding_design)$rank,
  full_rank = qr(confounding_design)$rank == ncol(confounding_design)),
  file.path(TABLE_DIR, "06_batch_histology_identifiability.csv"), row.names = FALSE)
write.csv(as.data.frame(table(discovery$gse_id, discovery$histology, discovery$C1_C2)),
  file.path(TABLE_DIR, "06_source_histology_original_group_counts.csv"), row.names = FALSE)

#Module 23. GEO observed metastasis and gated time-to-event analysis -----------
# Each source is analyzed separately. GSE21050 shares records with GSE71118;
# these are sensitivity annotations, not independent validation cohorts.
# Observed yes/no status is NOT a risk estimate at a fixed time or a survival endpoint.
geo_sources <- intersect(c("GSE71118", "GSE21050"), unique(all_metadata$gse_id))
geo <- all_metadata[all_metadata$gse_id %in% geo_sources, , drop = FALSE]
geo$original_group <- discovery$C1_C2[match(geo$record_group, discovery$record_group)]
geo$resolved_histology <- as.character(discovery$histology[match(geo$record_group, discovery$record_group)])
stopifnot(!anyNA(geo$original_group), !anyDuplicated(paste(geo$gse_id, geo$record_group)))
geo$metastasis_observed <- tolower(trimws(as.character(geo$metastasis)))
geo$metastasis_binary <- ifelse(geo$metastasis_observed == "yes", 1L,
  ifelse(geo$metastasis_observed == "no", 0L, NA_integer_))
geo$time_numeric <- suppressWarnings(as.numeric(geo$time))
write.csv(geo[, c("sampleid", "record_group", "gse_id", "original_group", "resolved_histology", "metastasis", "time", "metastasis_binary", "time_numeric")],
  file.path(TABLE_DIR, "07_GEO_observed_endpoint_audit.csv"), row.names = FALSE)
zero_geo <- geo[!is.na(geo$time_numeric) & geo$time_numeric <= 0, c("sampleid", "record_group", "gse_id", "original_group", "metastasis", "time")]
write.csv(zero_geo, file.path(TABLE_DIR, "07_GEO_nonpositive_times_for_review.csv"), row.names = FALSE)
binary_rows <- list(); rate_rows <- list()
for (series_id in geo_sources) for (subset_name in c("all", "LMS_only")) {
  keep <- geo$gse_id == series_id
  if (subset_name == "LMS_only") keep <- keep & geo$resolved_histology == "Leiomyosarcoma"
  source_data <- geo[keep, , drop = FALSE]
  d <- source_data[!is.na(source_data$metastasis_binary), , drop = FALSE]
  # Rows: C2,C1; columns: yes,no -> odds ratio for C2 vs C1.
  tab <- table(factor(d$original_group, levels = c("C2", "C1")), factor(d$metastasis_binary, levels = c(1, 0)))
  for (group in c("C1", "C2")) {
    n <- sum(tab[group, ]); y <- tab[group, "1"]
    ci <- if (n > 0) stats::binom.test(y, n)$conf.int else c(NA_real_, NA_real_)
    rate_rows[[length(rate_rows) + 1L]] <- data.frame(gse_id = series_id, subset = subset_name,
      original_group = group, source_n = sum(source_data$original_group == group), observed_n = n, yes_n = y,
      observed_proportion = if (n) y / n else NA_real_, CI_low = ci[1], CI_high = ci[2])
  }
  if (any(rowSums(tab) == 0L) || any(colSums(tab) == 0L)) {
    binary_rows[[length(binary_rows) + 1L]] <- data.frame(gse_id = series_id, subset = subset_name,
      n = nrow(d), OR_C2_vs_C1 = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_, status = "insufficient_variation")
  } else {
    test <- stats::fisher.test(tab)
    binary_rows[[length(binary_rows) + 1L]] <- data.frame(gse_id = series_id, subset = subset_name,
      n = nrow(d), OR_C2_vs_C1 = unname(test$estimate), CI_low = test$conf.int[1], CI_high = test$conf.int[2],
      p_value = test$p.value, status = "descriptive_observed_status_only")
  }
}
if (length(binary_rows)) {
  binary <- do.call(rbind, binary_rows); binary$BH_binary_family <- p.adjust(binary$p_value, "BH")
  write.csv(binary, file.path(TABLE_DIR, "07_GEO_observed_metastasis_Fisher.csv"), row.names = FALSE)
  write.csv(do.call(rbind, rate_rows), file.path(TABLE_DIR, "07_GEO_observed_metastasis_proportions.csv"), row.names = FALSE)
}
# A verified dictionary is optional. Never infer years from decimal values alone.
# Document why time zero means baseline disease; exclude only with verified policy.
dictionary_file <- file.path(CONFIG_DIR, "GEO_verified_endpoint_dictionary.csv")
dictionary_template <- data.frame(gse_id = c("GSE71118", "GSE21050"), endpoint = "MFS_candidate",
  time_field = "time", event_field = "metastasis", event_value = "yes", censor_value = "no",
  unit = "UNVERIFIED", origin = "UNVERIFIED", censoring_definition = "UNVERIFIED",
  zero_time_policy = "BLOCK", zero_time_evidence = "UNVERIFIED", field_mapping_evidence = "UNVERIFIED",
  source_doi = c("10.1016/j.ejca.2015.12.027", "10.1038/nm.2174"), verified = FALSE)
write.csv(dictionary_template, file.path(TABLE_DIR, "07_GEO_endpoint_dictionary_TEMPLATE.csv"), row.names = FALSE)
dictionary <- dictionary_template
if (file.exists(dictionary_file)) {
  dictionary <- data.table::fread(dictionary_file, data.table = FALSE)
  input_files <- c(input_files, dictionary_file)
  stopifnot(all(names(dictionary_template) %in% names(dictionary)), !anyDuplicated(dictionary$gse_id))
}
write.csv(dictionary, file.path(TABLE_DIR, "07_GEO_endpoint_dictionary_used.csv"), row.names = FALSE)
geo_status_rows <- list(); geo_survival_rows <- list(); geo_survival_audit <- list()
for (series_id in geo_sources) {
  spec <- dictionary[dictionary$gse_id == series_id, , drop = FALSE]
  ready <- nrow(spec) == 1L && isTRUE(as.logical(spec$verified))
  if (ready) ready <- isTRUE(all(!is.na(unlist(spec[c("origin", "censoring_definition", "field_mapping_evidence")])) &
    !unlist(spec[c("origin", "censoring_definition", "field_mapping_evidence")]) %in% c("", "UNVERIFIED")) &&
    spec$endpoint == "MFS" && spec$unit %in% c("days", "months", "years") &&
    all(c(spec$time_field, spec$event_field) %in% names(geo)) &&
    !is.na(spec$event_value) && !is.na(spec$censor_value) &&
    tolower(spec$event_value) != tolower(spec$censor_value) &&
    spec$zero_time_policy %in% c("block_if_present", "exclude_verified_baseline_disease"))
  if (!ready) {
    geo_status_rows[[length(geo_status_rows) + 1L]] <- data.frame(gse_id = series_id, subset = "all",
      status = "BLOCKED_endpoint_definition", detail = "Confirm field mapping, unit, origin, censoring and zero-time meaning with source evidence.")
    next
  }
  d0 <- geo[geo$gse_id == series_id, , drop = FALSE]
  d0$time_months <- suppressWarnings(as.numeric(d0[[spec$time_field]])) * switch(spec$unit, days = 12 / 365.25, months = 1, years = 12)
  event_text <- tolower(trimws(as.character(d0[[spec$event_field]])))
  d0$event <- ifelse(event_text == tolower(spec$event_value), 1L, ifelse(event_text == tolower(spec$censor_value), 0L, NA_integer_))
  zero_rows <- is.finite(d0$time_months) & d0$time_months == 0
  zero_evidence_ok <- !is.na(spec$zero_time_evidence) && !spec$zero_time_evidence %in% c("", "UNVERIFIED")
  # Exclusion applies only if every time-zero record is a verified baseline event.
  if (any(d0$time_months < 0, na.rm = TRUE) || (any(zero_rows) &&
    (spec$zero_time_policy != "exclude_verified_baseline_disease" || !zero_evidence_ok || any(is.na(d0$event[zero_rows])) || any(d0$event[zero_rows] != 1L)))) {
    geo_status_rows[[length(geo_status_rows) + 1L]] <- data.frame(gse_id = series_id, subset = "all",
      status = "BLOCKED_nonpositive_time", detail = "No epsilon replacement; baseline disease exclusion requires documented meaning and event=1.")
    next
  }
  d0$included <- is.finite(d0$time_months) & d0$time_months > 0 & !is.na(d0$event)
  geo_survival_audit[[series_id]] <- d0[, c("sampleid", "record_group", "gse_id", "original_group", "time_months", "event", "included")]
  for (subset_name in c("all", "LMS_only")) {
    keep <- d0$included
    if (subset_name == "LMS_only") keep <- keep & d0$resolved_histology == "Leiomyosarcoma"
    d <- droplevels(d0[keep, , drop = FALSE]); d$original_group <- factor(d$original_group, levels = c("C1", "C2"))
    group_n <- table(d$original_group); group_event <- tapply(d$event, d$original_group, sum)
    if (any(group_n < 5L) || anyNA(group_event) || any(group_event < 3L) || sum(d$event) < MIN_EVENTS_PER_PARAMETER) {
      geo_status_rows[[length(geo_status_rows) + 1L]] <- data.frame(gse_id = series_id, subset = subset_name,
        status = "insufficient_group_size_or_events", detail = "Require >=5 records and >=3 events per group, >=10 events overall.")
      next
    }
    warnings_fit <- character()
    fit <- try(withCallingHandlers(survival::coxph(Surv(time_months, event) ~ original_group, data = d, ties = "efron", x = TRUE),
      warning = function(w) { warnings_fit <<- c(warnings_fit, conditionMessage(w)); invokeRestart("muffleWarning") }), silent = TRUE)
    if (inherits(fit, "try-error") || length(warnings_fit) || any(!is.finite(coef(fit))) || any(!is.finite(vcov(fit)))) {
      geo_status_rows[[length(geo_status_rows) + 1L]] <- data.frame(gse_id = series_id, subset = subset_name,
        status = "fit_failed_or_unstable", detail = if (length(warnings_fit)) paste(warnings_fit, collapse = " | ") else "Inspect model inputs")
      next
    }
    ci <- exp(stats::confint(fit)); ph <- try(survival::cox.zph(fit), silent = TRUE)
    if (any(!is.finite(ci))) {
      geo_status_rows[[length(geo_status_rows) + 1L]] <- data.frame(gse_id = series_id, subset = subset_name,
        status = "unbounded_confidence_interval", detail = "Not reported as a stable estimate")
      next
    }
    geo_survival_rows[[length(geo_survival_rows) + 1L]] <- data.frame(gse_id = series_id, subset = subset_name,
      endpoint = spec$endpoint, n = nrow(d), events = sum(d$event), HR_C2_vs_C1 = exp(coef(fit)[1]),
      CI_low = ci[1, 1], CI_high = ci[1, 2], p_value = summary(fit)$coefficients[1, "Pr(>|z|)"],
      PH_p = if (inherits(ph, "try-error")) NA_real_ else ph$table[1, "p"],
      estimand = if (any(zero_rows)) "conditional_on_no_verified_baseline_metastasis" else "source_defined_time_to_event")
    geo_status_rows[[length(geo_status_rows) + 1L]] <- data.frame(gse_id = series_id, subset = subset_name,
      status = "completed_exploratory_source_specific", detail = "Overlapping source cohorts; no independent replication claim")
    saveRDS(fit, file.path(CACHE_DIR, paste0("07_GEO_", series_id, "_", subset_name, ".rds")))
  }
}
if (length(geo_survival_rows)) {
  geo_survival <- do.call(rbind, geo_survival_rows); geo_survival$BH_GEO_time_family <- p.adjust(geo_survival$p_value, "BH")
  write.csv(geo_survival, file.path(TABLE_DIR, "07_GEO_verified_survival.csv"), row.names = FALSE)
}
if (length(geo_survival_audit)) write.csv(do.call(rbind, geo_survival_audit), file.path(TABLE_DIR, "07_GEO_survival_inclusion_audit.csv"), row.names = FALSE)
geo_status <- if (length(geo_status_rows)) do.call(rbind, geo_status_rows) else data.frame(gse_id = NA, subset = "all", status = "no_candidate_sources", detail = "")
write.csv(geo_status, file.path(TABLE_DIR, "07_GEO_survival_status.csv"), row.names = FALSE)
log_message("GEO time-to-event status:", paste(unique(geo_status$status), collapse = "; "))

#Module 24. TCGA AUC uncertainty ----------------
# Fixed score and fixed sign. Bootstrap histology strata without refitting the
# signature; this CI conditions on the discovered signature and cohort composition.
scenario_histologies <- list(primary_174_definition = c("STLMS", "DDLPS", "UPS", "MFS", "SS"),
  plus_ULMS = c("STLMS", "ULMS", "DDLPS", "UPS", "MFS", "SS"), all_histologies = unique(tcga$histology_group),
  LMS_only = c("STLMS", "ULMS"), STLMS_only = "STLMS")
auc_rows <- list(); auc_draws <- list(); set.seed(SEED + 8L)
for (scenario in names(scenario_histologies)[1:3]) {
  d <- tcga[tcga$histology_group %in% scenario_histologies[[scenario]], , drop = FALSE]
  positive <- d$histology_group %in% c("STLMS", "ULMS"); np <- sum(positive); nn <- sum(!positive)
  if (np < 2L || nn < 2L) next
  estimate <- (sum(rank(d$score_z)[positive]) - np * (np + 1) / 2) / (np * nn)
  draws <- numeric(AUC_BOOTSTRAPS); strata_indices <- split(seq_len(nrow(d)), d$histology_group)
  for (b in seq_len(AUC_BOOTSTRAPS)) {
    sampled <- integer()
    for (ii in strata_indices) sampled <- c(sampled, ii[sample.int(length(ii), length(ii), replace = TRUE)])
    pos <- positive[sampled]
    draws[b] <- (sum(rank(d$score_z[sampled])[pos]) - np * (np + 1) / 2) / (np * nn)
  }
  ci <- quantile(draws, c(.025, .975), names = FALSE)
  auc_rows[[scenario]] <- data.frame(scenario = scenario, n = nrow(d), n_LMS = np, AUC = estimate,
    CI_low = ci[1], CI_high = ci[2], bootstrap_replicates = AUC_BOOTSTRAPS,
    role = "fixed_signature_histology_discrimination_not_prognosis_or_independent_cluster_replication")
  auc_draws[[scenario]] <- data.frame(scenario = scenario, replicate = seq_along(draws), AUC = draws)
}
if (length(auc_rows)) {
  write.csv(do.call(rbind, auc_rows), file.path(TABLE_DIR, "08_TCGA_AUC_bootstrap_CI.csv"), row.names = FALSE)
  write.csv(do.call(rbind, auc_draws), file.path(TABLE_DIR, "08_TCGA_AUC_bootstrap_draws.csv"), row.names = FALSE)
}
#Module 25. Before/after PCA, provenance and analysis readiness ----------------
variable_genes <- head(order(matrixStats::rowVars(expression_raw), decreasing = TRUE), TOP_VARIABLE_GENES)
pca_rows <- list(); pca_variance_rows <- list()
for (stage in c("before_ComBat", "after_ComBat")) {
  x <- if (stage == "before_ComBat") expression_raw else expression_corrected
  pca <- prcomp(t(x[variable_genes, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 2)
  explained <- pca$sdev^2 / sum(pca$sdev^2)
  pca_rows[[stage]] <- data.frame(sampleid = discovery$sampleid, stage = stage, PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    gse_id = discovery$gse_id, histology = discovery$histology, original_group = discovery$C1_C2)
  if ("platform_id" %in% names(discovery)) pca_rows[[stage]]$platform_id <- as.character(discovery$platform_id)
  pca_variance_rows[[stage]] <- data.frame(stage = stage, PC1_variance = explained[1], PC2_variance = explained[2])
}
pca_table <- do.call(rbind, pca_rows)
write.csv(pca_table, file.path(TABLE_DIR, "10_PCA_coordinates.csv"), row.names = FALSE)
write.csv(do.call(rbind, pca_variance_rows), file.path(TABLE_DIR, "10_PCA_explained_variance.csv"), row.names = FALSE)
for (color_variable in intersect(c("gse_id", "platform_id", "histology", "original_group"), names(pca_table))) {
  pca_table$color_group <- pca_table[[color_variable]]
  p <- ggplot2::ggplot(pca_table, ggplot2::aes(PC1, PC2, color = color_group)) +
    ggplot2::geom_point(alpha = .7, size = 1.5) + ggplot2::facet_wrap(~ stage, scales = "free") +
    ggplot2::theme_bw(base_size = 10) + ggplot2::labs(color = color_variable,
      title = "Representative records: before and after ComBat",
      caption = "Same genes selected by raw variance; separate PCA fits. Axes are not directly aligned.\nVisual mixing does not establish preservation of biology when source and histology are confounded.")
  save_figure(p, paste0("10_PCA_", color_variable), 11, 5.2)
}
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo.txt"))

#Module 26. 04_Representative_policies output paths ----
SEED <- 20261102L; set.seed(SEED)
OUTPUT_DIR <- file.path(PRIMARY_RESULTS, "04_Representative_policies")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables"); PLOT_DIR <- file.path(OUTPUT_DIR, "figures")
CACHE_DIR <- file.path(OUTPUT_DIR, "cache"); LOG_DIR <- file.path(OUTPUT_DIR, "session")
for (directory in c(TABLE_DIR, PLOT_DIR, CACHE_DIR, LOG_DIR)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
input_files <- character()

RECORD_GROUPS_FILE <- file.path(CODE_DIR, "Inputs", "Primary", "Record_groups.csv")
FIXED_TARGETS_FILE <- file.path(CODE_DIR, "Inputs", "Primary", "Fixed_original_targets.csv")
SIGNATURE_FILE <- file.path(CODE_DIR, "Inputs", "Primary", "Revised_signature_genes.csv")
PARENT_MATRIX_FILE <- file.path(V20_DIR, "cache", "discovery_deduplicated_uncorrected.rds")
#Module 27. Reconstruct identical input records and frozen genes --------------
signature_table <- read.csv(SIGNATURE_FILE, check.names = FALSE)
all_records <- read.csv(RECORD_GROUPS_FILE, check.names = FALSE)
targets <- read.csv(FIXED_TARGETS_FILE, check.names = FALSE)
record_columns <- c("sampleid", "gse_id", "geo_accession", "expr_col_original", "diagnosis_group", "legacy_cluster", "record_group", "priority", "representative", "selection_reason")
target_columns <- c("sampleid", "record_group", "C1_C2", "diagnosis_group", "gse_id")
stopifnot("GeneSymbol" %in% names(signature_table), all(record_columns %in% names(all_records)),
  all(target_columns %in% names(targets)))
all_records <- all_records[, record_columns, drop = FALSE]
targets <- targets[, target_columns, drop = FALSE]
revised_genes <- trimws(as.character(signature_table$GeneSymbol))
stopifnot(length(revised_genes) == 30L, !anyNA(revised_genes), !anyDuplicated(revised_genes), all(nzchar(revised_genes)))
prepared_series <- unique(as.character(all_records$gse_id))
prepared_files <- file.path(PROJECT_DIR, "Input", "GEO_prepared", prepared_series,
  paste0(prepared_series, "_ExpressionMatrix_Annotated.csv"))
if (!all(file.exists(prepared_files))) stop("Missing prepared expression matrices:\n", paste(prepared_files[!file.exists(prepared_files)], collapse = "\n"))
parent_path <- PARENT_MATRIX_FILE
parent_raw <- readRDS(parent_path)
stopifnot(is.matrix(parent_raw), is.numeric(parent_raw), !anyDuplicated(rownames(parent_raw)),
  !anyDuplicated(colnames(parent_raw)), all(is.finite(parent_raw)), nrow(targets) == 629L,
  !anyDuplicated(targets$record_group), !anyDuplicated(all_records$sampleid),
  !anyNA(all_records$record_group), setequal(all_records$record_group, targets$record_group),
  setequal(targets$sampleid, colnames(parent_raw)), all(revised_genes %in% rownames(parent_raw)))
stopifnot(nrow(all_records) == 1115L, nrow(parent_raw) == 20287L,
  !anyNA(targets[c("sampleid", "record_group", "C1_C2", "gse_id")]), !anyNA(all_records[c("sampleid", "gse_id", "expr_col_original", "legacy_cluster", "record_group")]),
  all(targets$C1_C2 %in% c("C1", "C2")), sum(targets$C1_C2 == "C2") == 37L,
  length(unique(all_records$gse_id)) == 9L)
targets <- targets[order(targets$record_group), , drop = FALSE]
for (g in targets$record_group) {
  labels <- unique(all_records$legacy_cluster[all_records$record_group == g])
  if (length(labels) != 1L || labels != targets$C1_C2[targets$record_group == g]) stop("Original labels disagree within linked group: ", g)
}
frozen_genes <- rownames(parent_raw); matrix_list <- list()
for (series_id in unique(all_records$gse_id)) {
  matrix_path <- file.path(PROJECT_DIR, "Input", "GEO_prepared", series_id, paste0(series_id, "_ExpressionMatrix_Annotated.csv"))
  if (!file.exists(matrix_path)) stop("Missing prepared matrix: ", matrix_path)
  source_data <- data.table::fread(matrix_path, data.table = FALSE, check.names = FALSE)
  gene_column <- if ("Gene_Symbol" %in% names(source_data)) "Gene_Symbol" else names(source_data)[1]
  symbols <- trimws(as.character(source_data[[gene_column]]))
  sm <- all_records[all_records$gse_id == series_id, , drop = FALSE]
  stopifnot(!anyNA(symbols), !anyDuplicated(symbols), all(nzchar(symbols)),
    all(frozen_genes %in% symbols), all(sm$expr_col_original %in% names(source_data)))
  m <- as.matrix(source_data[match(frozen_genes, symbols), sm$expr_col_original, drop = FALSE])
  storage.mode(m) <- "double"; dimnames(m) <- list(frozen_genes, sm$sampleid)
  stopifnot(all(is.finite(m))); matrix_list[[series_id]] <- m
}
all_raw <- do.call(cbind, matrix_list)[, all_records$sampleid, drop = FALSE]
parent_difference <- max(abs(all_raw[, targets$sampleid, drop = FALSE] - parent_raw[, targets$sampleid, drop = FALSE]))
if (!is.finite(parent_difference) || parent_difference > 1e-7) stop("Prepared matrices differ from parent inputs; do not compare changed input versions.")
rm(matrix_list, source_data, m)

#Module 28. Refit correction and compare matched groups under three policies --
policies <- c("current_representatives", "prefer_GSE21050", "prefer_GSE16382")
selection_rows <- list(); support_rows <- list(); label_rows <- list(); recovery_rows <- list(); score_rows <- list()
for (policy in policies) {
  if (policy == "current_representatives") {
    chosen <- all_records[match(targets$sampleid, all_records$sampleid), , drop = FALSE]
  } else {
    first_source <- sub("prefer_", "", policy)
    source_order <- unique(c(first_source, "GSE71118", "GSE21050", "GSE16382", sort(unique(all_records$gse_id))))
    ranked <- all_records[order(all_records$record_group, match(all_records$gse_id, source_order), all_records$sampleid), , drop = FALSE]
    chosen <- ranked[!duplicated(ranked$record_group), , drop = FALSE]
    chosen <- chosen[match(targets$record_group, chosen$record_group), , drop = FALSE]
  }
  stopifnot(nrow(chosen) == 629L, identical(as.character(chosen$record_group), as.character(targets$record_group)))
  chosen$fixed_target <- targets$C1_C2
  selection_rows[[policy]] <- data.frame(policy = policy, chosen)
  expr <- all_raw[, chosen$sampleid, drop = FALSE]
  batch <- factor(chosen$gse_id)
  corrected <- sva::ComBat(expr, batch = batch, mod = NULL, par.prior = TRUE,
    prior.plots = FALSE, mean.only = any(table(batch) == 1L))
  stopifnot(identical(dimnames(corrected), dimnames(expr)), all(is.finite(corrected)))
  z <- t(scale(t(corrected[revised_genes, , drop = FALSE])))
  stopifnot(all(is.finite(z)))
  score_rows[[policy]] <- data.frame(policy = policy, record_group = chosen$record_group,
    sampleid = chosen$sampleid, fixed_target = chosen$fixed_target, myogenic_score = colMeans(z))
  for (space in c("full_gene", "top_variable_PC")) {
    if (space == "full_gene") coords <- t(corrected) else {
      genes <- head(order(matrixStats::rowVars(corrected), decreasing = TRUE), 3000L)
      pc <- prcomp(t(corrected[genes, , drop = FALSE]), center = TRUE, scale. = FALSE,
        rank. = min(30L, ncol(corrected) - 1L, length(genes)))
      coords <- pc$x
      write.csv(data.frame(GeneSymbol = rownames(corrected)[genes]), file.path(TABLE_DIR, paste0("07_features_", policy, ".csv")), row.names = FALSE)
    }
    centered <- sweep(coords, 2L, colMeans(coords), "-")
    gram <- tcrossprod(centered); norms <- diag(gram)
    squared <- outer(norms, norms, "+") - 2 * gram; squared[squared < 0] <- 0
    distances <- sqrt(squared); diag(distances) <- 0
    dimnames(distances) <- list(chosen$record_group, chosen$record_group)
    distance <- as.dist(distances); tree <- hclust(distance, method = "ward.D2")
    selected_support <- data.frame(policy = policy, space = space, k = 2:6, silhouette = NA_real_)
    for (k in 2:6) {
      labels <- cutree(tree, k = k)
      selected_support$silhouette[selected_support$k == k] <- mean(cluster::silhouette(labels, distance)[, "sil_width"])
      label_rows[[length(label_rows) + 1L]] <- data.frame(policy = policy, space = space, k = k,
        record_group = chosen$record_group, cluster = labels, fixed_target = chosen$fixed_target)
      target <- chosen$fixed_target == "C2"
      for (cl in sort(unique(labels))) {
        member <- labels == cl
        recovery_rows[[length(recovery_rows) + 1L]] <- data.frame(policy = policy, space = space, k = k, cluster = cl,
          n = sum(member), C2_overlap = sum(member & target), target_recall = sum(member & target) / sum(target),
          target_precision = sum(member & target) / sum(member))
      }
    }
    selected_support$selected_by_silhouette <- selected_support$k == selected_support$k[which.max(selected_support$silhouette)]
    support_rows[[length(support_rows) + 1L]] <- selected_support
  }
}
all_labels <- do.call(rbind, label_rows); all_scores <- do.call(rbind, score_rows)
comparison_rows <- list(); score_agreement_rows <- list()
for (policy in policies[-1]) {
  a <- all_scores[all_scores$policy == policies[1], , drop = FALSE]
  b <- all_scores[all_scores$policy == policy, , drop = FALSE]; b <- b[match(a$record_group, b$record_group), , drop = FALSE]
  stopifnot(identical(as.character(a$record_group), as.character(b$record_group)), nrow(a) == 629L)
  score_agreement_rows[[policy]] <- data.frame(policy = policy, n = nrow(a),
    Pearson = cor(a$myogenic_score, b$myogenic_score), Spearman = cor(a$myogenic_score, b$myogenic_score, method = "spearman"))
  for (space in c("full_gene", "top_variable_PC")) for (k in 2:6) {
    a <- all_labels[all_labels$policy == policies[1] & all_labels$space == space & all_labels$k == k, , drop = FALSE]
    b <- all_labels[all_labels$policy == policy & all_labels$space == space & all_labels$k == k, , drop = FALSE]
    b <- b[match(a$record_group, b$record_group), , drop = FALSE]
    stopifnot(identical(as.character(a$record_group), as.character(b$record_group)))
    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(policy = policy, space = space, k = k,
      ARI_vs_current_representatives = mclust::adjustedRandIndex(a$cluster, b$cluster))
  }
}
write.csv(do.call(rbind, selection_rows), file.path(TABLE_DIR, "07_representative_selections.csv"), row.names = FALSE)
write.csv(do.call(rbind, support_rows), file.path(TABLE_DIR, "07_clustering_support.csv"), row.names = FALSE)
write.csv(all_labels, file.path(TABLE_DIR, "07_group_aligned_cluster_labels.csv"), row.names = FALSE)
write.csv(do.call(rbind, recovery_rows), file.path(TABLE_DIR, "07_fixed_C2_recovery.csv"), row.names = FALSE)
write.csv(all_scores, file.path(TABLE_DIR, "07_fixed_signature_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, comparison_rows), file.path(TABLE_DIR, "07_partition_agreement.csv"), row.names = FALSE)
write.csv(do.call(rbind, score_agreement_rows), file.path(TABLE_DIR, "07_score_agreement.csv"), row.names = FALSE)
representative_status <- "completed"

#Module 29. Summarize policy changes and save returnable outputs ---------------
selection <- do.call(rbind, selection_rows)
selection_summary <- list()
for (policy in policies) {
  chosen <- selection[selection$policy == policy, , drop = FALSE]
  chosen <- chosen[match(targets$record_group, chosen$record_group), , drop = FALSE]
  stopifnot(identical(as.character(chosen$record_group), as.character(targets$record_group)))
  selection_summary[[policy]] <- data.frame(policy = policy, connected_groups = nrow(chosen),
    changed_representatives = sum(chosen$sampleid != targets$sampleid),
    changed_C2_representatives = sum(chosen$sampleid != targets$sampleid & targets$C1_C2 == "C2"))
}
write.csv(do.call(rbind, selection_summary), file.path(TABLE_DIR, "07_selection_summary.csv"), row.names = FALSE)
support <- do.call(rbind, support_rows)
write.csv(support[support$selected_by_silhouette, , drop = FALSE], file.path(TABLE_DIR, "07_selected_k.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "R_sessionInfo.txt"))

#Module 30. Completion ----
writeLines(capture.output(sessionInfo()), file.path(PRIMARY_RESULTS, "R_sessionInfo.txt"))
message("Primary expression analysis completed: ", PRIMARY_RESULTS)
