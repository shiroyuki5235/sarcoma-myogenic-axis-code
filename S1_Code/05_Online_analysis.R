# Sarcoma online discovery, enrichment and TCGA validation
# Run 04_Online_inputs.R first. This script completes the analysis in one R session.
# 처음에는 전 자료를 온라인으로 다운로드합니다. 이후에는 이 버전의 검증 입력만 이어서 사용합니다.

#Module 00. 새 온라인 프로젝트 ----
# 04와 05에서 같은 PROJECT_DIR을 사용합니다. 빈 폴더에서도 시작할 수 있습니다.
# 기본 경로는 Windows R의 사용자 Documents 아래입니다. D:/Sarcoma_Online 등으로 변경 가능합니다.
PROJECT_DIR <- Sys.getenv("SARCOMA_ONLINE_DIR", unset = file.path(path.expand("~"), "Sarcoma_Online_v2_1_0"))
PIPELINE_VERSION <- "2.1.0"
options(stringsAsFactors = FALSE, warn = 1); set.seed(123)
dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
# 기존 Input/GEO_prepared와 이전 Output/cache를 읽지 않습니다.
INPUT_DIR <- file.path(PROJECT_DIR, "Input_Online_v2_1_0")
STATE_FILE <- file.path(INPUT_DIR, "sarcoma_online_run.rds")
gse_list <- c("GSE71118", "GSE21050", "GSE16382", "GSE34620", "GSE142162", "GSE66533", "GSE20196", "GSE34800", "GSE94321")
INPUT_GEO_PREP_DIR <- file.path(INPUT_DIR, "GEO_prepared")
INPUT_GEO_RAW_DIR <- file.path(INPUT_DIR, "GEO_raw")
INPUT_TCGA_DIR <- file.path(INPUT_DIR, "TCGA")
GEO_CACHE <- file.path(INPUT_DIR, "GEO_metadata")
READY_DIR <- file.path(INPUT_DIR, "GEO_ready_v2")
GENESET_DIR <- file.path(INPUT_DIR, "GeneSets")

if (!file.exists(STATE_FILE)) stop("04_Online_inputs.R을 먼저 실행하십시오: ", STATE_FILE)
run_state <- readRDS(STATE_FILE)
if (!identical(run_state$pipeline_version, PIPELINE_VERSION) || !identical(run_state$project_dir, PROJECT_DIR)) stop("실행 기록의 버전/프로젝트 경로가 다릅니다. 이 프로젝트에서 04를 다시 실행하십시오.")
OUTPUT_DIR <- run_state$output_dir
if (!dir.exists(OUTPUT_DIR)) stop("기록된 출력 폴더가 없습니다. 04를 다시 실행하십시오.")

#Module 01. 04에서 준비한 입력과 실행 상태 확인 ----
if (!isTRUE(run_state$inputs_ready) || length(run_state$ready_files) != length(gse_list) || !all(file.exists(run_state$ready_files))) stop("GEO 9개 입력이 준비되지 않았습니다. 04의 pending_inputs.csv를 확인하십시오.")
if (!identical(unname(tools::md5sum(run_state$ready_files)), run_state$ready_md5)) stop("준비 입력이 바뀌었습니다. 04를 다시 실행하십시오.")
# Each analysis uses a fresh result folder.
{
  previous_output <- OUTPUT_DIR
  OUTPUT_DIR <- file.path(PROJECT_DIR, paste0("Output_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_Sarcoma_OnlineAnalysis"))
  if (dir.exists(OUTPUT_DIR)) stop("동일 시각 출력 폴더가 있습니다. 잠시 뒤 다시 실행하십시오.")
  dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE)
  diagnostic_files <- file.path(previous_output, "tables", c("input_readiness.csv", "pending_inputs.csv", "input_download_manifest.csv", "network_attempts.csv", "tcga_input_status.csv"))
  for (f in diagnostic_files[file.exists(diagnostic_files)]) file.copy(f, file.path(OUTPUT_DIR, "tables", basename(f)))
  run_state$output_dir <- OUTPUT_DIR
  saveRDS(run_state, STATE_FILE)
}
saveRDS(list(started = as.character(Sys.time()), version = PIPELINE_VERSION), file.path(OUTPUT_DIR, "05_started.rds"))
exprs_list <- list(); meta_list <- list()
for (gse_id in gse_list) {
  ready <- readRDS(file.path(READY_DIR, paste0(gse_id, "_ready.rds")))
  if (!identical(ready$gse_id, gse_id) || !identical(ready$pipeline_version, PIPELINE_VERSION)) stop("입력 버전/GSE 불일치: ", gse_id)
  if (!all(file.exists(ready$input_paths)) || !identical(unname(tools::md5sum(ready$input_paths)), ready$input_md5)) stop("CSV/metadata가 바뀌었습니다. 04 재실행 필요: ", gse_id)
  x <- ready$expression; m <- ready$metadata
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) || nrow(x) < 2L || ncol(x) < 2L) stop("유효하지 않은 발현행렬: ", gse_id)
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) || !identical(colnames(x), m$sampleid) || anyNA(m$geo_accession)) stop("표본/유전자 ID 검사 실패: ", gse_id)
  exprs_list[[gse_id]] <- x; meta_list[[gse_id]] <- m
}
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
LOG_DIR <- file.path(OUTPUT_DIR, "logs")
ADDON_TAB_DIR <- file.path(OUTPUT_DIR, "signature")
TCGA_TAB_DIR <- file.path(OUTPUT_DIR, "TCGA")
for (p in c(TABLE_DIR, LOG_DIR, ADDON_TAB_DIR, TCGA_TAB_DIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
#Module 02. Required analysis packages and settings ----
cran_packages <- c("data.table", "matrixStats", "cluster", "dplyr", "tidyr")
bioc_packages <- c("limma", "sva", "fgsea", "BiocParallel")
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos = "https://cloud.r-project.org")
if (length(missing_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
stopifnot(all(vapply(c(cran_packages, bioc_packages), requireNamespace, logical(1), quietly = TRUE)))
data.table::setDTthreads(2L)
LOG_FILE <- file.path(LOG_DIR, "05_online_analysis.log")
K_GRID <- 2:6; K_FINAL_OVERRIDE <- NA_integer_
BOOT_B <- 30L; BOOT_PROP <- 0.80; KSUPPORT_TOP_VAR <- 3000L
KSUPPORT_N_PC <- 30L; KSUPPORT_LOG_EVERY <- 5L
LOW_VAR_QUANTILE <- 0.05; DEG_FDR_CUTOFF <- 0.05; DEG_LOGFC_CUTOFF <- 1
MIN_SAMPLES_PER_GSE <- 5L; TOP_MARKER_N_FOR_SCORE <- 30L
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


#Module 03. 공개 코드의 표본 선별 규칙 적용 ----
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
if (any(!nzchar(meta_df$sampleid)) || anyDuplicated(meta_df$sampleid)) stop("Metadata sampleid 중복 또는 누락")

repeated_gsm <- duplicated(meta_df$geo_accession) | duplicated(meta_df$geo_accession, fromLast = TRUE)
write.csv(meta_df[repeated_gsm, , drop = FALSE], file.path(TABLE_DIR, "repeated_GSM_across_GSE_REVIEW.csv"), row.names = FALSE)
# GSE_GSM은 record ID입니다. 서로 다른 GSM이 동일 환자인지는 별도 linkage 자료가 필요합니다.
all_sample_ids <- unlist(lapply(exprs_list, colnames), use.names = FALSE)

sample_master <- data.frame(
  sampleid = all_sample_ids,
  gse_id   = sub("_.*$", "", all_sample_ids),
  stringsAsFactors = FALSE
)
sample_master <- dplyr::left_join(sample_master, meta_df, by = c("sampleid", "gse_id"))

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
    "title", "source_name_ch1", "diagnosis_h", "histology_h", "subtype_h", "fusion", "pathology",
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
# Fusion-positive BCOR-CCNB3 sarcoma is distinct from Ewing; "Ewing negative" is not a positive diagnosis.
bcor_gsm <- paste0("GSM", 855583:855592)
sample_master$diagnosis_group[grepl("bcor[^a-z0-9]*ccnb3", dx_txt) |
  sample_master$geo_accession %in% bcor_gsm] <- "BCOR-CCNB3 sarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("leiomyosarcoma|leiomyosarcome|\\blms\\b", dx_txt)] <- "Leiomyosarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("ewing", dx_txt) &
  !grepl("ewing[ -]*(negative|negative sample)|non[ -]*ewing|not[ -]*ewing", dx_txt)] <- "Ewing sarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("synovial", dx_txt)] <- "Synovial sarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("rhabdomyosarcoma", dx_txt)] <- "Rhabdomyosarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("liposarcoma.*dediff|dediff.*liposarcoma|dedifferentiated liposarcoma", dx_txt)] <- "Dedifferentiated liposarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("liposarcoma", dx_txt)] <- "Liposarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("undifferentiated sarcoma|pleomorphic sarcoma|\\bups\\b", dx_txt)] <- "Undifferentiated/pleomorphic sarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("mpnst|malignant peripheral nerve sheath", dx_txt)] <- "MPNST"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("osteosarcoma", dx_txt)] <- "Osteosarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("chondrosarcoma", dx_txt)] <- "Chondrosarcoma"
sample_master$diagnosis_group[is.na(sample_master$diagnosis_group) & grepl("sarcoma|rhabdo|synovial|lipo|mpnst|osteosarcoma|chondrosarcoma", dx_txt)] <- "Other sarcoma"
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

# 공개 코드의 unknown-species 포함 정책을 유지하고 별도 표시합니다.
sample_master$organism_unknown <- is.na(sample_master$organism2) | sample_master$organism2 == ""
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

gse_post <- sample_master
gse_post <- dplyr::filter(gse_post, gse_id %in% analysis_gse_tbl$gse_id)
gse_post <- dplyr::group_by(gse_post, gse_id)
gse_post <- dplyr::summarise(gse_post, n_total = dplyr::n(),
    n_keep_initial = sum(sample_decision == "KEEP", na.rm = TRUE),
    .groups = "drop")

small_gse <- gse_post$gse_id[gse_post$n_keep_initial > 0 & gse_post$n_keep_initial < MIN_SAMPLES_PER_GSE]
if (length(small_gse) > 0) {
  sample_master$sample_decision[sample_master$gse_id %in% small_gse & sample_master$sample_decision == "KEEP"] <- "DROP"
  sample_master$sample_reason[sample_master$gse_id %in% small_gse & grepl("^passed_conservative_sample_filter$", sample_master$sample_reason)] <-
    paste0("dropped_posthoc_small_gse_lt_", MIN_SAMPLES_PER_GSE)
}

gse_post2 <- sample_master
gse_post2 <- dplyr::filter(gse_post2, gse_id %in% analysis_gse_tbl$gse_id)
gse_post2 <- dplyr::group_by(gse_post2, gse_id, cohort_role)
gse_post2 <- dplyr::summarise(gse_post2, n_total = dplyr::n(),
    n_keep = sum(sample_decision == "KEEP", na.rm = TRUE),
    n_drop = sum(sample_decision == "DROP", na.rm = TRUE),
    .groups = "drop")
gse_post2 <- dplyr::arrange(gse_post2, dplyr::desc(n_keep), gse_id)

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
#Module 04. 공통 유전자를 같은 순서로 정렬한 뒤 결합 ----
common_genes <- Reduce(intersect, lapply(exprs_keep_list, rownames))
if (length(common_genes) < 1000L) stop("공통 유전자가 1000개 미만입니다. gene ID/annotation을 확인하십시오.")
gene_overlap <- data.frame(gse_id = names(exprs_keep_list), n_input_genes = vapply(exprs_keep_list, nrow, integer(1)), n_common_genes = length(common_genes))
write.csv(gene_overlap, file.path(TABLE_DIR, "common_gene_alignment.csv"), row.names = FALSE)
write.csv(data.frame(GeneSymbol = common_genes), file.path(TABLE_DIR, "common_genes_before_filter.csv"), row.names = FALSE)
for (g in names(exprs_keep_list)) exprs_keep_list[[g]] <- exprs_keep_list[[g]][common_genes, , drop = FALSE]
stopifnot(all(vapply(exprs_keep_list, function(x) identical(rownames(x), common_genes), logical(1))))
exprs_all <- do.call(cbind, exprs_keep_list)
if (anyDuplicated(colnames(exprs_all)) || any(!is.finite(exprs_all))) stop("통합행렬 표본 중복/비유한값")
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

#Module 05. 명시적인 ComBat 보정 ----
if (length(unique(batch_info)) >= 2) {
  exprs_combat <- sva::ComBat(dat = exprs_all, batch = batch_info, par.prior = TRUE, prior.plots = FALSE)
  batch_method <- "sva::ComBat"
} else {
  exprs_combat <- exprs_all
  batch_method <- "single_batch_no_correction"
}
if (any(!is.finite(exprs_combat)) || !identical(dimnames(exprs_all), dimnames(exprs_combat))) stop("ComBat 산출물 검사 실패")
write.csv(data.frame(method = batch_method), file.path(TABLE_DIR, "batch_correction_method.csv"), row.names = FALSE)

storage.mode(exprs_combat) <- "double"
saveRDS(exprs_all, file.path(OUTPUT_DIR, "exprs_all_curated_prefiltered.rds"))
saveRDS(exprs_combat, file.path(OUTPUT_DIR, "exprs_combat_curated.rds"))

msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("Preparing clustering objects"))
cat(msg, "\n")
cat(msg, "\n", file = LOG_FILE, append = TRUE)

#Module 06. 군집 수 평가와 최종 군집에 동일한 PCA 공간 사용 ----

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
    metric = c("n_ksupport_genes", "n_ksupport_pcs", "subsampling_B", "subsampling_prop"),
    value  = c(n_k_features, n_pc_use, BOOT_B, BOOT_PROP),
    stringsAsFactors = FALSE
  ),
  file.path(TABLE_DIR, "cluster_k_support_runinfo.csv"),
  row.names = FALSE
)

dist_mat_k <- dist(t(exprs_k_eval))
hc_k <- hclust(dist_mat_k, method = "ward.D2")
hc_full <- hc_k
saveRDS(list(genes = rownames(exprs_k), pcs = exprs_k_eval, hclust = hc_k), file.path(OUTPUT_DIR, "clustering_space.rds"))
if (any(K_GRID < 2L) || max(K_GRID) >= floor(ncol(exprs_k_eval) * BOOT_PROP)) stop("K_GRID/표본 수/부분추출 비율 확인 필요")
if (!is.na(K_FINAL_OVERRIDE) && !K_FINAL_OVERRIDE %in% K_GRID) stop("K_FINAL_OVERRIDE는 평가한 K_GRID에 포함되어야 합니다.")

k_metrics <- list()
for (k in K_GRID) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("k-support start: k=%d", k))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
  cl_k <- cutree(hc_k, k = k)
  sil <- cluster::silhouette(cl_k, dist_mat_k)
  sil_avg <- if (!is.null(sil)) mean(sil[, "sil_width"], na.rm = TRUE) else NA_real_

  set.seed(123 + k)
  n_sub <- max(5L, floor(ncol(exprs_k_eval) * BOOT_PROP))
  ari_vec <- rep(NA_real_, BOOT_B)

  for (b in seq_len(BOOT_B)) {
    idx_b <- sort(sample.int(ncol(exprs_k_eval), size = n_sub, replace = FALSE))
    sub_expr <- exprs_k_eval[, idx_b, drop = FALSE]
    sub_hc <- hclust(dist(t(sub_expr)), method = "ward.D2")

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
        } else {
          ari_vec[b] <- 1
        }
      }
    }

    if (b %% max(1L, KSUPPORT_LOG_EVERY) == 0L || b == BOOT_B) {
      msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf("k-support subsampling [%s]: %d / %d complete", paste0("k=", k), b, BOOT_B))
      cat(msg, "\n")
      cat(msg, "\n", file = LOG_FILE, append = TRUE)
    }
  }

  stab_df <- data.frame(
    k = k,
    subsampling_n = BOOT_B,
    subsampling_prop = BOOT_PROP,
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
if (any(!is.finite(k_metrics_df$avg_silhouette)) || any(!is.finite(k_metrics_df$ari_mean))) stop("군집 수 평가에 유효하지 않은 지표가 있습니다.")
k_metrics_df$rank_silhouette <- rank(-k_metrics_df$avg_silhouette, ties.method = "min")
k_metrics_df$rank_ari_mean   <- rank(-k_metrics_df$ari_mean, ties.method = "min")
k_metrics_df$rank_sum        <- k_metrics_df$rank_silhouette + k_metrics_df$rank_ari_mean

if (is.na(K_FINAL_OVERRIDE)) {
  k_pick <- k_metrics_df
  k_pick <- dplyr::arrange(k_pick, rank_sum, dplyr::desc(ari_mean), dplyr::desc(avg_silhouette), k)
  k_pick <- dplyr::slice(k_pick, 1)
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

cluster_assignments <- cutree(hc_full, k = K_FINAL)
cluster_df <- data.frame(
  SampleID = names(cluster_assignments),
  Cluster  = factor(cluster_assignments, levels = seq_len(K_FINAL), labels = paste0("C", seq_len(K_FINAL))),
  stringsAsFactors = FALSE
)
write.csv(cluster_df, file.path(TABLE_DIR, paste0("cluster_assignments_k", K_FINAL, ".csv")), row.names = FALSE)

sample_master2 <- sample_master
sample_master2 <- dplyr::filter(sample_master2, sample_decision == "KEEP")
sample_master2 <- dplyr::left_join(sample_master2, cluster_df, by = c("sampleid" = "SampleID"))

write.csv(sample_master2, file.path(TABLE_DIR, "sample_master_metadata_curated_with_cluster.csv"), row.names = FALSE)

cluster_sizes <- as.data.frame(table(cluster_df$Cluster), stringsAsFactors = FALSE)
names(cluster_sizes) <- c("Cluster", "N")
write.csv(cluster_sizes, file.path(TABLE_DIR, paste0("cluster_sizes_k", K_FINAL, ".csv")), row.names = FALSE)

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

#Module 07. 군집별 DEG: 군집 생성 자료를 이용한 기술적 비교 ----
if (any(table(groups) < 2L)) stop("2개 미만 표본의 군집이 있어 DEG를 진행할 수 없습니다.")
contr <- limma::makeContrasts(contrasts = unname(contrast_strings), levels = design)
colnames(contr) <- names(contrast_strings)
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
  rank_lists[[coef_name]] <- sort(setNames(tt$t[is.finite(tt$t)], tt$GeneSymbol[is.finite(tt$t)]), decreasing = TRUE)

  write.csv(tt, file.path(TABLE_DIR, paste0("DEG_", coef_name, ".csv")), row.names = FALSE)

  sig_lists[[coef_name]] <- tt$GeneSymbol[is.finite(tt$adj.P.Val) & is.finite(tt$logFC) & tt$adj.P.Val < DEG_FDR_CUTOFF & abs(tt$logFC) > DEG_LOGFC_CUTOFF]

}

#Module 08. Freeze the discovery-only marker signature ----
stopifnot(identical(dimnames(exprs_all), dimnames(exprs_combat)),
  !anyDuplicated(sample_master2$sampleid), setequal(colnames(exprs_combat), sample_master2$sampleid))
sample_master2 <- sample_master2[match(colnames(exprs_combat), sample_master2$sampleid), , drop = FALSE]
stopifnot(identical(sample_master2$sampleid, colnames(exprs_combat)), !anyNA(sample_master2$Cluster))
minor_order <- order(cluster_sizes$N, cluster_sizes$Cluster)
minor_cluster <- as.character(cluster_sizes$Cluster[minor_order[1L]])
major_cluster <- as.character(cluster_sizes$Cluster[order(-cluster_sizes$N, cluster_sizes$Cluster)[1L]])
if (minor_cluster == major_cluster) stop("Distinct minor and major clusters are required.")
# Minor denotes cluster size; it does not establish a molecular subtype.
write.csv(data.frame(minor_cluster, major_cluster), file.path(ADDON_TAB_DIR, "minor_major_definition.csv"), row.names = FALSE)
minor_contrast <- paste0(minor_cluster, "_vsRest")
minor_markers_df <- deg_tables[[minor_contrast]]
minor_markers_df <- minor_markers_df[is.finite(minor_markers_df$adj.P.Val) &
  minor_markers_df$adj.P.Val < DEG_FDR_CUTOFF & minor_markers_df$logFC > DEG_LOGFC_CUTOFF, , drop = FALSE]
# Preserve limma topTable(sort.by = "P") order used in the reported online run.
minor_markers <- head(unique(minor_markers_df$GeneSymbol), TOP_MARKER_N_FOR_SCORE)
marker_sd <- apply(exprs_combat[minor_markers, , drop = FALSE], 1L, stats::sd)
minor_markers <- minor_markers[is.finite(marker_sd) & marker_sd > 0]
if (length(minor_markers) < 5L) stop("Fewer than five usable discovery markers.")
frozen_signature <- list(version = PIPELINE_VERSION, genes = minor_markers,
  minor_cluster = minor_cluster, major_cluster = major_cluster,
  minor_centroid = rowMeans(exprs_combat[minor_markers, sample_master2$Cluster == minor_cluster, drop = FALSE]),
  major_centroid = rowMeans(exprs_combat[minor_markers, sample_master2$Cluster == major_cluster, drop = FALSE]),
  selection = "top significant discovery up-markers; at most 30; nonconstant genes")
saveRDS(frozen_signature, file.path(OUTPUT_DIR, "frozen_minor_signature.rds"))
write.csv(data.frame(GeneSymbol = minor_markers), file.path(ADDON_TAB_DIR, "signature_genes.csv"), row.names = FALSE)
discovery_signature_z <- t(scale(t(exprs_combat[minor_markers, , drop = FALSE])))
stopifnot(all(is.finite(discovery_signature_z)))
write.csv(data.frame(sampleid = colnames(exprs_combat), Cluster = sample_master2$Cluster,
  minor_cluster_score = colMeans(discovery_signature_z)), file.path(ADDON_TAB_DIR, "discovery_signature_scores.csv"), row.names = FALSE)

#Module 09. Hallmark and Reactome GSEA: sequential execution ----
# A fixed seed per collection makes results independent of removed exploratory plots/tests.
# Multilevel is preferred; permutation fallback is labelled and never silently substituted.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
gsea_results <- list(); gsea_status <- list()
for (collection in c("Hallmark", "Reactome")) {
  geneset_file <- file.path(GENESET_DIR, paste0(collection, ".rds"))
  if (!file.exists(geneset_file)) stop("Missing gene sets: ", geneset_file)
  pathways <- readRDS(geneset_file); stats_vec <- rank_lists[[minor_contrast]]
  if (!is.list(pathways) || !length(pathways) || is.null(names(pathways))) stop("Invalid gene-set collection.")
  if (anyDuplicated(names(stats_vec)) || any(!is.finite(stats_vec))) stop("Invalid GSEA ranks.")
  set.seed(123L + match(collection, c("Hallmark", "Reactome")))
  method <- "fgseaMultilevel"; fallback_reason <- ""
  message("GSEA ", collection, ": ", minor_contrast, " (SerialParam)")
  out <- tryCatch(fgsea::fgseaMultilevel(pathways = pathways, stats = stats_vec, minSize = 10L, maxSize = 500L,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)), error = function(e) e)
  if (inherits(out, "error")) {
    fallback_reason <- conditionMessage(out); message("Multilevel failed: ", fallback_reason)
    method <- "fgseaSimple_20000_permutations"
    set.seed(123L + match(collection, c("Hallmark", "Reactome")))
    out <- fgsea::fgseaSimple(pathways = pathways, stats = stats_vec, minSize = 10L, maxSize = 500L,
      nperm = 20000L, BPPARAM = BiocParallel::SerialParam(progressbar = FALSE))
  }
  out <- as.data.frame(out)
  if (!nrow(out) || all(!is.finite(out$NES))) stop("GSEA returned no usable pathways: ", collection)
  out$leadingEdge <- vapply(out$leadingEdge, paste, character(1), collapse = ";")
  out <- out[order(out$padj, -abs(out$NES)), , drop = FALSE]
  write.csv(out, file.path(TABLE_DIR, paste0("GSEA_", collection, "_", minor_contrast, ".csv")), row.names = FALSE)
  gsea_results[[collection]] <- out
  gsea_status[[collection]] <- data.frame(collection, contrast = minor_contrast, method, fallback_reason,
    n_pathways = nrow(out), nonfinite_NES = sum(!is.finite(out$NES)), nonfinite_padj = sum(!is.finite(out$padj)))
  if (collection == "Hallmark" && (!"HALLMARK_MYOGENESIS" %in% out$pathway ||
    any(!is.finite(out$padj[out$pathway == "HALLMARK_MYOGENESIS"])))) stop("Myogenesis enrichment did not complete.")
}
write.csv(dplyr::bind_rows(gsea_status), file.path(TABLE_DIR, "GSEA_execution_summary.csv"), row.names = FALSE)
#Module 10. Load TCGA expression and clinical inputs ----
TCGA_DIR <- file.path(INPUT_TCGA_DIR, "sarc_tcga_pub")
expr_file <- file.path(TCGA_DIR, "data_mrna_seq_v2_rsem.txt")
clin_sample_file <- file.path(TCGA_DIR, "data_clinical_sample.txt")
clin_patient_file <- file.path(TCGA_DIR, "data_clinical_patient.txt")
tcga_files <- c(expr_file, clin_sample_file, clin_patient_file)
if (!all(file.exists(tcga_files))) stop("TCGA 입력 누락. 04를 실행하거나 input_download_manifest/README의 위치에 파일을 저장하십시오.")
write.csv(data.frame(file = tcga_files, md5 = unname(tools::md5sum(tcga_files))), file.path(TCGA_TAB_DIR, "tcga_source_files.csv"), row.names = FALSE)
data.table::setDTthreads(2L)
expr_dt <- data.table::fread(expr_file, data.table = FALSE, check.names = FALSE)

# 임상 파일의 설명 주석 줄 수가 달라져도 실제 header를 찾아 읽습니다.
clinical <- list()
for (f in c(clin_sample_file, clin_patient_file)) {
  lines <- readLines(f, n = 100L, warn = FALSE)
  header_idx <- which(!grepl("^#", lines) & nzchar(trimws(lines)))[1L]
  if (is.na(header_idx)) stop("Clinical header를 찾지 못했습니다: ", f)
  clinical[[basename(f)]] <- data.table::fread(f, skip = header_idx - 1L, data.table = FALSE, check.names = FALSE,
                                             na.strings = c("NA", "[Not Available]", "[Not Applicable]", ""))
}
clin_sample <- clinical[[basename(clin_sample_file)]]
clin_patient <- clinical[[basename(clin_patient_file)]]
if (!all(c("SAMPLE_ID", "PATIENT_ID") %in% names(clin_sample)) || !"PATIENT_ID" %in% names(clin_patient)) stop("Clinical ID 열 누락")
if (anyNA(clin_sample$SAMPLE_ID) || anyNA(clin_sample$PATIENT_ID) || anyNA(clin_patient$PATIENT_ID) ||
    anyDuplicated(clin_sample$SAMPLE_ID) || anyDuplicated(clin_patient$PATIENT_ID)) stop("Clinical ID 누락/중복. 임의로 첫 행을 선택하지 않습니다.")
if (!"Hugo_Symbol" %in% names(expr_dt)) stop("Hugo_Symbol 없는 TCGA expression 형식")
expr_dt$Hugo_Symbol <- trimws(as.character(expr_dt$Hugo_Symbol))
expr_dt <- expr_dt[!is.na(expr_dt$Hugo_Symbol) & nzchar(expr_dt$Hugo_Symbol), , drop = FALSE]
# 원본의 첫 symbol 행 유지 정책을 명시하고 중복 symbol을 기록합니다.
duplicate_symbols <- unique(expr_dt$Hugo_Symbol[duplicated(expr_dt$Hugo_Symbol)])
write.csv(data.frame(Hugo_Symbol = duplicate_symbols), file.path(TCGA_TAB_DIR, "tcga_duplicate_symbols_first_row_retained.csv"), row.names = FALSE)
expr_dt <- expr_dt[!duplicated(expr_dt$Hugo_Symbol), , drop = FALSE]
sample_cols <- grep("^TCGA-", names(expr_dt), value = TRUE)
if (length(sample_cols) < 2L || anyDuplicated(sample_cols)) stop("TCGA expression sample 열을 확인하십시오.")
expr_mat <- as.matrix(expr_dt[, sample_cols, drop = FALSE])
rownames(expr_mat) <- expr_dt$Hugo_Symbol
if (!is.numeric(expr_mat) || any(!is.finite(expr_mat)) || any(expr_mat < 0)) stop("RSEM 값은 비음수의 유한한 숫자여야 합니다. 다른 데이터 종류로 대체하지 마십시오.")
storage.mode(expr_mat) <- "double"
expr_mat <- log2(expr_mat + 1)

#Module 11. Clinical 결합과 histology 분류 ----
tcga_meta <- merge(clin_sample, clin_patient, by = "PATIENT_ID", all.x = TRUE, sort = FALSE, suffixes = c("_sample", "_patient"))
if (anyDuplicated(tcga_meta$SAMPLE_ID)) stop("Clinical 결합으로 SAMPLE_ID 중복 발생")
need_meta <- c("CANCER_TYPE_DETAILED", "ONCOTREE_CODE", "SHORT_HISTOLOGY", "TUMOR_TISSUE_SITE")
for (nm in need_meta) {
  value <- rep(NA_character_, nrow(tcga_meta))
  for (source_nm in intersect(c(nm, paste0(nm, "_sample"), paste0(nm, "_patient")), names(tcga_meta))) {
    source_value <- trimws(as.character(tcga_meta[[source_nm]]))
    take <- (is.na(value) | !nzchar(value)) & !is.na(source_value) & nzchar(source_value)
    value[take] <- source_value[take]
  }
  tcga_meta[[nm]] <- value
}
tcga_meta$histology_label <- tcga_meta$CANCER_TYPE_DETAILED
for (nm in c("SHORT_HISTOLOGY", "ONCOTREE_CODE")) {
  take <- is.na(tcga_meta$histology_label) | !nzchar(tcga_meta$histology_label)
  tcga_meta$histology_label[take] <- tcga_meta[[nm]][take]
}
# 세부 표기와 tissue site를 함께 사용합니다. uterine LMS를 일반 LMS보다 먼저 판정합니다.
txt <- tolower(paste(ifelse(is.na(tcga_meta$histology_label), "", tcga_meta$histology_label),
                     ifelse(is.na(tcga_meta$SHORT_HISTOLOGY), "", tcga_meta$SHORT_HISTOLOGY),
                     ifelse(is.na(tcga_meta$ONCOTREE_CODE), "", tcga_meta$ONCOTREE_CODE)))
site_txt <- tolower(ifelse(is.na(tcga_meta$TUMOR_TISSUE_SITE), "", tcga_meta$TUMOR_TISSUE_SITE))
is_lms <- grepl("leiomyosarcoma|leiomyosarcome|\\blms\\b|\\bstlms\\b|\\bulms\\b", txt)
is_uterine <- grepl("uterine|uterus|\\bulms\\b", txt) | grepl("uter", site_txt)
tcga_meta$histology_group <- NA_character_
tcga_meta$histology_group[grepl("dedifferentiated liposarcoma|\\bddls\\b|\\bddlps\\b", txt)] <- "DDLPS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("undifferentiated pleomorphic sarcoma|malignant fibrous histiocytoma|high-grade spindle|\\bups\\b", txt)] <- "UPS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("myxofibrosarcoma|\\bmfs\\b", txt)] <- "MFS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("synovial sarcoma|\\bss\\b|\\bsyns\\b", txt)] <- "SS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("malignant peripheral nerve sheath|mpnst", txt)] <- "MPNST"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & is_lms & is_uterine] <- "ULMS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & is_lms & !is_uterine] <- "STLMS"
write.csv(tcga_meta, file.path(TCGA_TAB_DIR, "tcga_histology_mapping_review.csv"), row.names = FALSE)
write.csv(tcga_meta[is.na(tcga_meta$histology_group), , drop = FALSE], file.path(TCGA_TAB_DIR, "tcga_histology_unclassified.csv"), row.names = FALSE)

write.csv(data.frame(SAMPLE_ID = setdiff(colnames(expr_mat), tcga_meta$SAMPLE_ID)), file.path(TCGA_TAB_DIR, "expression_samples_without_clinical.csv"), row.names = FALSE)
keep_expr_ids <- intersect(colnames(expr_mat), tcga_meta$SAMPLE_ID)
if (length(keep_expr_ids) < 10L) stop("TCGA expression/clinical 공통 표본이 너무 적습니다.")
expr_mat <- expr_mat[, keep_expr_ids, drop = FALSE]
tcga_meta <- tcga_meta[match(colnames(expr_mat), tcga_meta$SAMPLE_ID), , drop = FALSE]
stopifnot(identical(colnames(expr_mat), tcga_meta$SAMPLE_ID), !anyNA(tcga_meta$PATIENT_ID))
if (anyDuplicated(tcga_meta$PATIENT_ID)) {
  duplicated_patients <- unique(tcga_meta$PATIENT_ID[duplicated(tcga_meta$PATIENT_ID)])
  write.csv(tcga_meta[tcga_meta$PATIENT_ID %in% duplicated_patients, , drop = FALSE], file.path(TCGA_TAB_DIR, "repeated_patient_samples_REVIEW.csv"), row.names = FALSE)
  stop("환자별 복수 표본이 있습니다. 대표 표본 선정 규칙을 정한 뒤 진행해야 합니다.")
}

#Module 12. 고정 gene list의 외부 플랫폼 매핑 ----
minor_cluster <- frozen_signature$minor_cluster
major_cluster <- frozen_signature$major_cluster
fixed_genes <- frozen_signature$genes
if (anyDuplicated(fixed_genes) || length(fixed_genes) < 5L) stop("고정 signature 형식 오류")
gene_mapping <- data.frame(GeneSymbol = fixed_genes, present_in_TCGA = fixed_genes %in% rownames(expr_mat), stringsAsFactors = FALSE)
minor_genes <- fixed_genes[fixed_genes %in% rownames(expr_mat)]
# 빠진 상위 marker를 다음 순위 marker로 채우지 않습니다.
tcga_sd <- apply(expr_mat[minor_genes, , drop = FALSE], 1L, stats::sd)
constant_genes <- minor_genes[!is.finite(tcga_sd) | tcga_sd <= 0]
common_centroid_genes <- setdiff(minor_genes, constant_genes)
gene_mapping$nonconstant_in_TCGA <- gene_mapping$GeneSymbol %in% common_centroid_genes
write.csv(gene_mapping, file.path(TCGA_TAB_DIR, "frozen_signature_gene_mapping.csv"), row.names = FALSE)
if (length(common_centroid_genes) < 5L) stop("고정 signature의 유효 TCGA 공통 유전자가 5개 미만입니다.")
minor_centroid <- frozen_signature$minor_centroid[common_centroid_genes]
major_centroid <- frozen_signature$major_centroid[common_centroid_genes]
if (any(!is.finite(c(minor_centroid, major_centroid))) || stats::sd(minor_centroid) == 0 || stats::sd(major_centroid) == 0) stop("Discovery centroid 비유한값 또는 분산 0")
expr_tcga_sig <- expr_mat[common_centroid_genes, , drop = FALSE]
# 원본 정책: TCGA 전체 공통 표본 안에서 유전자별 z-score를 한 번 계산합니다.
# 이 점수는 cohort-relative이며 단일 환자 적용이 검증된 고정 임상 점수가 아닙니다.
expr_tcga_z <- t(scale(t(expr_tcga_sig)))
if (any(!is.finite(expr_tcga_z))) stop("TCGA z-score 비유한값")
minor_score <- colMeans(expr_tcga_z)
minor_cor <- major_cor <- rep(NA_real_, ncol(expr_tcga_sig))
for (i in seq_len(ncol(expr_tcga_sig))) {
  if (stats::sd(expr_tcga_sig[, i]) > 0) {
    minor_cor[i] <- cor(expr_tcga_sig[, i], minor_centroid, method = "pearson")
    major_cor[i] <- cor(expr_tcga_sig[, i], major_centroid, method = "pearson")
  }
}

#Module 13. 검증 표본별 결과 ----
val_df <- data.frame(SAMPLE_ID = colnames(expr_tcga_sig), minor_cluster_score = as.numeric(minor_score),
  cor_minor_centroid = minor_cor, cor_major_centroid = major_cor,
  predicted_axis = ifelse(minor_cor >= major_cor, "minor_like", "major_like"), stringsAsFactors = FALSE)
meta_add <- tcga_meta[match(val_df$SAMPLE_ID, tcga_meta$SAMPLE_ID), , drop = FALSE]
meta_add <- meta_add[, setdiff(names(meta_add), names(val_df)), drop = FALSE]
val_df <- cbind(val_df, meta_add)
val_df$is_STLMS <- ifelse(val_df$histology_group == "STLMS", 1L, 0L)
val_df$is_ULMS <- ifelse(val_df$histology_group == "ULMS", 1L, 0L)
val_df$is_MPNST <- ifelse(val_df$histology_group == "MPNST", 1L, 0L)
write.csv(val_df, file.path(TCGA_TAB_DIR, "tcga_validation_per_sample_full.csv"), row.names = FALSE)
write.csv(data.frame(GeneSymbol = common_centroid_genes), file.path(TCGA_TAB_DIR, "tcga_validation_signature_genes_used.csv"), row.names = FALSE)

#Module 14. 원본의 세 검증 시나리오 ----
scenario_list <- list(
  primary_non_uterine_non_mpnst = c("STLMS", "DDLPS", "UPS", "MFS", "SS"),
  sensitivity_plus_ulms = c("STLMS", "ULMS", "DDLPS", "UPS", "MFS", "SS"),
  sensitivity_plus_ulms_plus_mpnst = c("STLMS", "ULMS", "DDLPS", "UPS", "MFS", "SS", "MPNST")
)

scenario_summary <- list()
for (scenario_nm in names(scenario_list)) {
  allowed <- scenario_list[[scenario_nm]]
  sub_df <- val_df[val_df$histology_group %in% allowed, , drop = FALSE]
  if (nrow(sub_df) == 0) next
  sub_df$is_target_lms <- ifelse(sub_df$histology_group == "STLMS", 1L, 0L)
  if (scenario_nm != "primary_non_uterine_non_mpnst") sub_df$is_target_lms <- ifelse(sub_df$histology_group %in% c("STLMS", "ULMS"), 1L, 0L)

  pos <- sub_df$minor_cluster_score[sub_df$is_target_lms == 1]
  neg <- sub_df$minor_cluster_score[sub_df$is_target_lms == 0]

  auc_value <- NA_real_
  if (length(pos) > 0 && length(neg) > 0) {
    all_scores <- c(pos, neg)
    all_labels <- c(rep(1, length(pos)), rep(0, length(neg)))
    r <- rank(all_scores, ties.method = "average")
    u <- sum(r[all_labels == 1]) - length(pos) * (length(pos) + 1) / 2
    auc_value <- u / (length(pos) * length(neg))
  }


  scenario_summary[[scenario_nm]] <- data.frame(
    scenario = scenario_nm,
    n_samples = nrow(sub_df),
    n_centroid_classified = sum(!is.na(sub_df$predicted_axis)),
    n_target_lms = sum(sub_df$is_target_lms == 1),
    n_other = sum(sub_df$is_target_lms == 0),
    mean_score_target = if (length(pos) > 0) mean(pos, na.rm = TRUE) else NA_real_,
    mean_score_other = if (length(neg) > 0) mean(neg, na.rm = TRUE) else NA_real_,
    auc = auc_value,
    stringsAsFactors = FALSE
  )

  hist_sum <- sub_df
  hist_sum <- dplyr::group_by(hist_sum, histology_group)
  hist_sum <- dplyr::summarise(hist_sum, n = dplyr::n(), mean_score = mean(minor_cluster_score, na.rm = TRUE), sd_score = sd(minor_cluster_score, na.rm = TRUE), mean_cor_minor = mean(cor_minor_centroid, na.rm = TRUE), prop_minor_like = mean(predicted_axis == "minor_like", na.rm = TRUE), .groups = "drop")
  hist_sum <- dplyr::arrange(hist_sum, dplyr::desc(mean_score), dplyr::desc(n))
  write.csv(hist_sum, file.path(TCGA_TAB_DIR, paste0("tcga_histology_summary_", scenario_nm, ".csv")), row.names = FALSE)


}

if (length(scenario_summary) > 0) write.csv(dplyr::bind_rows(scenario_summary), file.path(TCGA_TAB_DIR, "tcga_validation_scenario_summary.csv"), row.names = FALSE)


#Module 15. Completion ----
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "05_R_sessionInfo.txt"))
print(cluster_sizes); print(dplyr::bind_rows(scenario_summary)); print(dplyr::bind_rows(gsea_status))
saveRDS(list(version = PIPELINE_VERSION, finished = as.character(Sys.time()),
  signature_genes = minor_markers, ready_md5 = run_state$ready_md5), file.path(OUTPUT_DIR, "05_complete.rds"))
message("Online analysis completed: ", OUTPUT_DIR)
