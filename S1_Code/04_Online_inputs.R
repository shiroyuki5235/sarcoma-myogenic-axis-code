# Sarcoma online input preparation
# Run 04_Online_inputs.R, then 05_Online_analysis.R.
# 처음에는 전 자료를 온라인으로 다운로드합니다. 이후에는 이 버전의 검증 입력만 이어서 사용합니다.

#Module 00. 새 온라인 프로젝트 ----
# 04와 05에서 같은 PROJECT_DIR을 사용합니다. 빈 폴더에서도 시작할 수 있습니다.
# 기본 경로는 Windows R의 사용자 Documents 아래입니다. D:/Sarcoma_Online 등으로 변경 가능합니다.
CODE_DIR <- normalizePath(Sys.getenv("SARCOMA_CODE_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
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
FROZEN_GENESET_DIR <- file.path(CODE_DIR, "Inputs", "Primary_project", "Input", "GeneSets")

#Module 01. 실행 폴더와 패키지 ----
DOWNLOAD_TIMEOUT <- 7200L
INSTALL_MISSING_PACKAGES <- TRUE
PREPARE_GENESETS <- TRUE
options(timeout = max(DOWNLOAD_TIMEOUT, getOption("timeout", 60)))
for (p in c(INPUT_DIR, INPUT_GEO_PREP_DIR, INPUT_GEO_RAW_DIR, INPUT_TCGA_DIR, GEO_CACHE, READY_DIR, GENESET_DIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
if (file.exists(STATE_FILE)) {
  run_state <- readRDS(STATE_FILE)
  if (!identical(run_state$pipeline_version, PIPELINE_VERSION)) stop("Incompatible input-state version; choose a new project folder.")
  OUTPUT_DIR <- run_state$output_dir
  if (!identical(run_state$project_dir, PROJECT_DIR) || !dir.exists(OUTPUT_DIR)) {
    OUTPUT_DIR <- file.path(PROJECT_DIR, paste0("Output_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_Sarcoma_Online"))
    run_state$project_dir <- PROJECT_DIR; run_state$output_dir <- OUTPUT_DIR
  }
} else {
  OUTPUT_DIR <- file.path(PROJECT_DIR, paste0("Output_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_Sarcoma_Online"))
  run_state <- list(pipeline_version = PIPELINE_VERSION, project_dir = PROJECT_DIR, output_dir = OUTPUT_DIR,
                    inputs_ready = FALSE, ready_files = character(), ready_md5 = character())
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
TABLE_DIR <- file.path(OUTPUT_DIR, "tables"); LOG_DIR <- file.path(OUTPUT_DIR, "logs")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(LOG_DIR, "01_prepare_online.log")
run_state$inputs_ready <- FALSE
saveRDS(run_state, STATE_FILE)
cat("새 온라인 프로젝트: ", PROJECT_DIR, "\n결과 폴더: ", OUTPUT_DIR, "\n", sep = "")

cran_core <- c("curl", "data.table", "dplyr", "R.utils", "digest")
bioc_core <- c("Biobase", "GEOquery", "affy", "affyio", "org.Hs.eg.db", "hgu133plus2cdf", "hgu133plus2.db", "AnnotationDbi")
cran_install <- cran_core; bioc_install <- bioc_core
if (INSTALL_MISSING_PACKAGES) {
  absent <- cran_install[!vapply(cran_install, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) tryCatch(install.packages(absent, repos = "https://cloud.r-project.org"), error = function(e) message(conditionMessage(e)))
  absent <- bioc_install[!vapply(bioc_install, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    tryCatch(BiocManager::install(absent, ask = FALSE, update = FALSE), error = function(e) message(conditionMessage(e)))
  }
}
package_status <- data.frame(package = c(cran_core, bioc_core), role = "required")
package_status$available <- vapply(package_status$package, requireNamespace, logical(1), quietly = TRUE)
write.csv(package_status, file.path(LOG_DIR, "package_status.csv"), row.names = FALSE)
absent <- package_status$package[package_status$role == "required" & !package_status$available]
if (length(absent)) stop("필수 패키지 설치를 완료하지 못했습니다: ", paste(absent, collapse = ", "))

manifest <- data.frame(gse_id = gse_list, stringsAsFactors = FALSE)
manifest$raw_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", sub("[0-9]{1,3}$", "nnn", gse_list), "/", gse_list, "/suppl/", gse_list, "_RAW.tar")
manifest$metadata_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse_list, "&targ=gsm&view=brief&form=text")
manifest$series_matrix_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", sub("[0-9]{1,3}$", "nnn", gse_list), "/", gse_list, "/matrix/", gse_list, "_series_matrix.txt.gz")
manifest$raw_destination <- file.path(INPUT_GEO_RAW_DIR, gse_list, paste0(gse_list, "_RAW.tar"))
manifest$metadata_destination <- file.path(GEO_CACHE, paste0(gse_list, "_GSM_metadata.soft"))
manifest$expression_csv <- file.path(INPUT_GEO_PREP_DIR, gse_list, paste0(gse_list, "_ExpressionMatrix_Annotated.csv"))
manifest$metadata_status <- "not_attempted"; manifest$raw_status <- "not_attempted"
manifest$metadata_error <- ""; manifest$raw_error <- ""
manifest$n_metadata_samples <- 0L; manifest$n_CEL <- 0L
manifest$metadata_source_file <- ""
write.csv(manifest, file.path(TABLE_DIR, "input_download_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "01_sessionInfo_start.txt"))

#Module 02. 온라인 다운로드: 파일별 성공 확인과 제한된 재시도 ----
# Windows의 기본 curl.exe를 우선 사용합니다. 서버의 403을 재시도 사유로 쓰지 않습니다.
# 403은 해당 파일에만 적용하고 다른 코호트 요청을 임의로 차단하지 않습니다.
# 429/인증/프록시 인증 문제는 해당 호스트의 이번 실행 요청을 보류합니다.
CURL_EXE <- unname(Sys.which(if (.Platform$OS.type == "windows") "curl.exe" else "curl"))
if (.Platform$OS.type == "windows") {
  windows_curl <- file.path(Sys.getenv("SystemRoot", unset = "C:/Windows"), "System32", "curl.exe")
  if (file.exists(windows_curl)) CURL_EXE <- windows_curl
}
download_state <- new.env(parent = emptyenv())
download_state$paused_hosts <- character(); download_state$attempts <- list()

# 공통 I/O만 함수로 묶고, 전처리/분석은 각 Module 본문에서 실행합니다.
fetch_online <- function(url, destination) {
  if (!grepl("^https://", url)) stop("HTTPS 주소만 허용합니다: ", url)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  receipt_file <- paste0(destination, ".download.rds")
  if (file.exists(destination) && file.exists(receipt_file)) {
    receipt <- tryCatch(readRDS(receipt_file), error = function(e) NULL)
    if (is.list(receipt) && identical(receipt$version, PIPELINE_VERSION) && identical(receipt$url, url) &&
        identical(receipt$md5, unname(tools::md5sum(destination)))) return(invisible(destination))
  }
  # 이 프로젝트에서 완료된 다운로드임을 확인할 수 없는 파일은 별도로 보존합니다.
  if (file.exists(destination)) {
    archived <- paste0(destination, ".unverified_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    if (!file.rename(destination, archived)) stop("이전 파일 보존 실패: ", destination)
  }
  host <- sub("^https://([^/]+).*$", "\\1", url)
  if (host %in% download_state$paused_hosts) stop("해당 서버의 인증/요청 제한으로 이번 실행 요청 보류: ", url)
  part <- tempfile("transfer_", tmpdir = dirname(destination), fileext = ".part")
  on.exit(if (file.exists(part)) unlink(part), add = TRUE)
  label <- paste0(gsub("[^A-Za-z0-9_.-]", "_", basename(destination)), "_", length(download_state$attempts) + 1L)
  header_file <- file.path(LOG_DIR, paste0(label, "_headers.txt"))
  error_file <- file.path(LOG_DIR, paste0(label, "_transport.txt"))
  http_code <- NA_integer_; transport_code <- NA_integer_; detail <- ""
  started <- Sys.time(); transport <- if (nzchar(CURL_EXE)) "system_curl" else "R_curl"
  cat(format(started), " Download: ", url, "\n", sep = "")
  if (nzchar(CURL_EXE)) {
    # TLS 검증을 유지합니다. 전달 인자는 shQuote로 경로/URL을 인용합니다.
    args <- c("--location", "--fail", "--silent", "--show-error", "--globoff",
      "--proto", shQuote("=https"), "--proto-redir", shQuote("=https"),
      "--connect-timeout", "60", "--max-time", as.character(DOWNLOAD_TIMEOUT),
      "--retry", "2", "--retry-delay", "5", "--retry-max-time", "180",
      "--dump-header", shQuote(header_file), "--output", shQuote(part), "--url", shQuote(url))
    transport_code <- tryCatch(suppressWarnings(system2(CURL_EXE, args = args, stdout = error_file, stderr = error_file)),
                               error = function(e) { detail <<- conditionMessage(e); 127L })
    if (file.exists(header_file)) {
      headers <- readLines(header_file, warn = FALSE)
      status_lines <- grep("^HTTP/[0-9.]+[[:space:]]+[0-9]{3}", headers, value = TRUE)
      if (length(status_lines)) http_code <- as.integer(sub("^HTTP/[0-9.]+[[:space:]]+([0-9]{3}).*$", "\\1", tail(status_lines, 1L)))
    }
    if (file.exists(error_file)) detail <- paste(c(detail, readLines(error_file, warn = FALSE)), collapse = " | ")
  }
  # curl 프로그램이 없거나 전송 자체가 성립하지 않은 경우에만 R의 curl을 사용합니다.
  # 받은 HTTP 403/404/429를 같은 주소의 다른 전송 방식으로 반복 요청하지 않습니다.
  if (!nzchar(CURL_EXE) || (is.na(http_code) && transport_code != 0L)) {
    transport <- if (nzchar(CURL_EXE)) "system_curl_then_R_curl" else "R_curl"
    for (try_i in seq_len(3L)) {
      http_code <- NA_integer_
      if (file.exists(part)) unlink(part)
      h <- curl::new_handle(followlocation = TRUE, connecttimeout = 60, timeout = DOWNLOAD_TIMEOUT,
                            ssl_verifypeer = TRUE, ssl_verifyhost = 2L, protocols = 2L, redir_protocols = 2L)
      response <- tryCatch(curl::curl_fetch_disk(url, part, handle = h), error = function(e) {
        detail <<- paste(detail, conditionMessage(e), sep = " | "); NULL
      })
      if (!is.null(response)) {
        http_code <- as.integer(response$status_code); transport_code <- 0L
        writeBin(response$headers, header_file)
      } else { transport_code <- 1L }
      if (!is.na(http_code) && !http_code %in% c(408L, 500L, 502L, 503L, 504L)) break
      if (try_i < 3L) Sys.sleep(5 * try_i)
    }
  }
  ok <- identical(as.integer(transport_code), 0L) && !is.na(http_code) && http_code >= 200L && http_code < 300L &&
    file.exists(part) && isTRUE(file.info(part)$size > 0)
  if (!is.na(http_code) && http_code %in% c(401L, 407L, 429L)) download_state$paused_hosts <- unique(c(download_state$paused_hosts, host))
  download_state$attempts[[length(download_state$attempts) + 1L]] <- data.frame(time = as.character(started),
    url = url, destination = destination, transport = transport, transport_code = transport_code, http_code = http_code,
    transfer_ok = ok, seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), detail = detail)
  write.csv(dplyr::bind_rows(download_state$attempts), file.path(TABLE_DIR, "network_attempts.csv"), row.names = FALSE)
  if (!ok) stop("다운로드 실패 [HTTP ", if (is.na(http_code)) "unknown" else http_code, "]: ", url, " | ", detail)
  if (!file.rename(part, destination)) stop("다운로드 저장 실패: ", destination)
  saveRDS(list(version = PIPELINE_VERSION, url = url, md5 = unname(tools::md5sum(destination)),
               bytes = file.info(destination)$size, fetched_at = as.character(Sys.time())), receipt_file)
  cat("  Download complete: ", round(file.info(destination)$size / 1024^2, 1), " MiB\n", sep = "")
  invisible(destination)
}

# GEO의 documented brief SOFT에서 각 GSM의 metadata만 읽습니다.
# https://www.ncbi.nlm.nih.gov/geo/info/download.html
read_geo_brief <- function(filename) {
  lines <- readLines(filename, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE[[:space:]]*=[[:space:]]*GSM[0-9]+", lines)
  if (!length(starts)) stop("GSM SOFT metadata가 아닙니다: ", filename)
  ends <- c(starts[-1L] - 1L, length(lines)); rows <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    gsm_id <- trimws(sub("^\\^SAMPLE[[:space:]]*=[[:space:]]*", "", lines[starts[i]]))
    block <- lines[seq.int(starts[i] + 1L, ends[i])]
    block <- block[grepl("^!Sample_[^=]+ = ", block)]
    if (!length(block)) stop("빈 GSM metadata: ", gsm_id)
    sep_pos <- regexpr(" = ", block, fixed = TRUE)
    keys <- substring(block, 9L, sep_pos - 1L)
    values <- substring(block, sep_pos + 3L)
    keys <- make.unique(keys, sep = "__")
    record <- as.list(values); names(record) <- keys
    record$geo_accession <- gsm_id
    rows[[i]] <- as.data.frame(record, stringsAsFactors = FALSE, check.names = FALSE)
  }
  pheno <- dplyr::bind_rows(rows)
  if (anyDuplicated(pheno$geo_accession) || anyNA(pheno$geo_accession)) stop("GSM SOFT 중복/누락")
  if (!all(c("title", "organism_ch1", "platform_id") %in% names(pheno))) stop("GSM SOFT 필수 속성 누락")
  pheno
}

#Module 03. 모든 코호트의 metadata 온라인 수집 ----
# 작은 metadata를 먼저 받고, RAW는 별도 Module에서 모두 시도합니다.
geo_pheno_list <- list()
for (gse_id in gse_list) {
  i_gse <- match(gse_id, manifest$gse_id)
  cat("\nMetadata: ", gse_id, "\n", sep = "")
  pheno <- NULL; meta_errors <- character(); source_file <- ""
  brief_file <- manifest$metadata_destination[i_gse]
  metadata_error <- tryCatch({
    fetch_online(manifest$metadata_url[i_gse], brief_file)
    pheno <- read_geo_brief(brief_file)
    source_file <- brief_file
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(metadata_error)) {
    meta_errors <- c(meta_errors, paste("brief SOFT:", metadata_error))
    if (file.exists(brief_file)) {
      if (!file.rename(brief_file, paste0(brief_file, ".invalid_", format(Sys.time(), "%Y%m%d_%H%M%S")))) stop("Invalid metadata 보존 실패")
      if (file.exists(paste0(brief_file, ".download.rds"))) unlink(paste0(brief_file, ".download.rds"))
    }
  }
  if (is.null(pheno)) {
    matrix_file <- file.path(GEO_CACHE, paste0(gse_id, "_series_matrix.txt.gz"))
    metadata_error <- tryCatch({
      fetch_online(manifest$series_matrix_url[i_gse], matrix_file)
      geo_args <- list(filename = matrix_file, getGPL = FALSE, parseCharacteristics = FALSE)
      if ("returnType" %in% names(formals(GEOquery::getGEO))) geo_args$returnType <- "ExpressionSet"
      geo_obj <- do.call(GEOquery::getGEO, geo_args)
      if (!methods::is(geo_obj, "ExpressionSet")) stop("Series Matrix가 ExpressionSet으로 읽히지 않았습니다.")
      pheno <- as.data.frame(Biobase::pData(geo_obj)); source_file <- matrix_file
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(metadata_error)) {
      meta_errors <- c(meta_errors, paste("Series Matrix:", metadata_error))
      if (file.exists(matrix_file)) {
        if (!file.rename(matrix_file, paste0(matrix_file, ".invalid_", format(Sys.time(), "%Y%m%d_%H%M%S")))) stop("Invalid Series Matrix 보존 실패")
        if (file.exists(paste0(matrix_file, ".download.rds"))) unlink(paste0(matrix_file, ".download.rds"))
      }
    }
  }
  metadata_error <- tryCatch({
    if (is.null(pheno) || !is.data.frame(pheno) || !nrow(pheno)) stop("온라인 metadata 수집 실패")
    if (!"geo_accession" %in% names(pheno)) stop("metadata의 geo_accession 열 누락")
    for (cc in names(pheno)) pheno[[cc]] <- as.character(pheno[[cc]])
    # 반복 characteristics를 펼치고 원문 열도 보존합니다.
    char_cols <- grep("^characteristics", names(pheno), value = TRUE, ignore.case = TRUE)
    for (cc in char_cols) for (i in seq_len(nrow(pheno))) {
      value <- pheno[[cc]][i]
      if (is.na(value) || !nzchar(value)) next
      for (piece in trimws(strsplit(value, ";", fixed = TRUE)[[1]])) {
        colon <- regexpr(":", piece, fixed = TRUE)[1L]
        if (colon < 2L) next
        key <- tolower(gsub("[^A-Za-z0-9]+", "_", trimws(substr(piece, 1L, colon - 1L))))
        val <- trimws(substring(piece, colon + 1L))
        if (!nzchar(key) || !nzchar(val)) next
        if (!key %in% names(pheno)) pheno[[key]] <- NA_character_
        if (is.na(pheno[[key]][i]) || !nzchar(pheno[[key]][i])) pheno[[key]][i] <- val
      }
    }
    names(pheno) <- make.unique(gsub("^_|_$", "", gsub("_+", "_", gsub("[^a-z0-9_]+", "_", tolower(names(pheno))))), sep = "__dup")
    pheno$geo_accession <- toupper(trimws(pheno$geo_accession))
    if (anyNA(pheno$geo_accession) || anyDuplicated(pheno$geo_accession) || any(!grepl("^GSM[0-9]+$", pheno$geo_accession))) stop("metadata GSM 누락/중복/형식 오류")
    if ("series_id" %in% names(pheno)) {
      series_cols <- grep("^series_id($|_)", names(pheno), value = TRUE)
      membership <- rep(FALSE, nrow(pheno))
      for (cc in series_cols) membership <- membership | (!is.na(pheno[[cc]]) & pheno[[cc]] == gse_id)
      if (!all(membership)) stop("요청 GSE와 GSM metadata의 Series 연결 불일치")
    }
    if ("platform_id" %in% names(pheno) && any(!pheno$platform_id %in% c("GPL570", "GPL16311", "GPL17810"))) stop("예상하지 않은 microarray platform; 임의 CDF 처리를 중단합니다.")
    geo_pheno_list[[gse_id]] <- pheno
    saveRDS(pheno, file.path(GEO_CACHE, paste0(gse_id, "_geo_pheno.rds")))
    manifest$metadata_status[i_gse] <- "ready"
    manifest$metadata_source_file[i_gse] <- source_file
    manifest$n_metadata_samples[i_gse] <- nrow(pheno)
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(metadata_error)) {
    manifest$metadata_status[i_gse] <- "failed"
    meta_errors <- c(meta_errors, metadata_error)
    if (nzchar(source_file) && file.exists(source_file)) {
      if (!file.rename(source_file, paste0(source_file, ".invalid_", format(Sys.time(), "%Y%m%d_%H%M%S")))) stop("Invalid metadata 보존 실패")
      if (file.exists(paste0(source_file, ".download.rds"))) unlink(paste0(source_file, ".download.rds"))
    }
  }
  manifest$metadata_error[i_gse] <- paste(meta_errors, collapse = " | ")
  cat(gse_id, ": ", manifest$metadata_status[i_gse], " (", manifest$n_metadata_samples[i_gse], " samples)\n", sep = "")
  write.csv(manifest, file.path(TABLE_DIR, "input_download_manifest.csv"), row.names = FALSE)
  Sys.sleep(1)
}

#Module 04. 모든 코호트의 RAW CEL archive 온라인 수집 ----
raw_members_list <- list()
for (gse_id in gse_list) {
  i_gse <- match(gse_id, manifest$gse_id)
  cat("\nRAW CEL archive: ", gse_id, "\n", sep = "")
  raw_tar <- manifest$raw_destination[i_gse]
  raw_error <- tryCatch({
    fetch_online(manifest$raw_url[i_gse], raw_tar)
    members <- withCallingHandlers(utils::untar(raw_tar, list = TRUE), warning = function(w) stop(conditionMessage(w)))
    cel_members <- members[grepl("\\.cel(\\.gz)?$", members, ignore.case = TRUE)]
    if (!length(cel_members)) stop("RAW archive에 CEL 파일이 없습니다.")
    if (any(grepl("(^/|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))", cel_members))) stop("안전하지 않은 archive 경로")
    cel_members <- cel_members[order(!grepl("\\.gz$", cel_members, ignore.case = TRUE))]
    cel_members <- cel_members[!duplicated(tolower(sub("\\.gz$", "", cel_members, ignore.case = TRUE)))]
    raw_members_list[[gse_id]] <- cel_members
    manifest$raw_status[i_gse] <- "ready"; manifest$n_CEL[i_gse] <- length(cel_members)
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(raw_error)) {
    manifest$raw_status[i_gse] <- "failed"; manifest$raw_error[i_gse] <- raw_error
    if (file.exists(raw_tar)) {
      if (!file.rename(raw_tar, paste0(raw_tar, ".invalid_", format(Sys.time(), "%Y%m%d_%H%M%S")))) stop("Invalid RAW archive 보존 실패")
      if (file.exists(paste0(raw_tar, ".download.rds"))) unlink(paste0(raw_tar, ".download.rds"))
    }
    message(gse_id, ": ", raw_error)
  }
  write.csv(manifest, file.path(TABLE_DIR, "input_download_manifest.csv"), row.names = FALSE)
  Sys.sleep(1)
}

#Module 05. TCGA와 gene set 온라인 수집 ----
# Use the same fixed cBioPortal files as the historical primary analysis.
tcga_sources <- read.csv(file.path(CODE_DIR, "Inputs", "Required_input_files.csv"),
  check.names = FALSE, colClasses = "character", na.strings = NULL)
tcga_sources <- tcga_sources[tcga_sources$source == "cbioportal_fixed_commit", , drop = FALSE]
stopifnot(nrow(tcga_sources) == 3L, !anyDuplicated(tcga_sources$relative_path))
TCGA_DIR <- file.path(INPUT_TCGA_DIR, "sarc_tcga_pub")
tcga_files <- file.path(TCGA_DIR, basename(tcga_sources$relative_path))
tcga_error <- tryCatch({
  for (input_index in seq_len(nrow(tcga_sources))) {
    if (!file.exists(tcga_files[input_index])) fetch_online(tcga_sources$url[input_index], tcga_files[input_index])
    if (!identical(digest::digest(file = tcga_files[input_index], algo = "sha256"), tcga_sources$sha256[input_index]))
      stop("TCGA input differs from the frozen analysis input: ", basename(tcga_files[input_index]))
  }
  if (!any(grepl("Hugo_Symbol", readLines(tcga_files[1L], n = 5L, warn = FALSE), fixed = TRUE))) stop("TCGA RSEM header 오류")
  NULL
}, error = function(e) conditionMessage(e))
write.csv(data.frame(url = tcga_sources$url, file = tcga_files, exists = file.exists(tcga_files),
                     status = if (is.null(tcga_error)) "ready" else "failed", detail = if (is.null(tcga_error)) "" else tcga_error),
          file.path(TABLE_DIR, "tcga_input_status.csv"), row.names = FALSE)
geneset_rows <- list()
# Reported enrichment uses archived gene sets; no current-database lookup is needed.
for (gs_kind in c("Hallmark", "Reactome")) {
  gs_error <- tryCatch({
    if (!PREPARE_GENESETS) stop("PREPARE_GENESETS=FALSE")
    for (suffix in c(".rds", "_source_table.rds")) {
      frozen_file <- file.path(FROZEN_GENESET_DIR, paste0(gs_kind, suffix))
      if (!file.exists(frozen_file)) stop("Missing archived MSigDB file: ", frozen_file)
      if (!file.copy(frozen_file, file.path(GENESET_DIR, basename(frozen_file)), overwrite = TRUE)) stop("Could not copy archived gene sets.")
    }
    frozen_source <- readRDS(file.path(GENESET_DIR, paste0(gs_kind, "_source_table.rds")))
    stopifnot(identical(unique(as.character(frozen_source$db_version)), "2026.1.Hs"))
    gs_check <- readRDS(file.path(GENESET_DIR, paste0(gs_kind, ".rds")))
    if (!is.list(gs_check) || !length(gs_check) || is.null(names(gs_check))) stop("Invalid archived gene-set list")
    NULL
  }, error = function(e) conditionMessage(e))
  geneset_rows[[gs_kind]] <- data.frame(collection = gs_kind, status = if (is.null(gs_error)) "ready" else "failed",
                                        detail = if (is.null(gs_error)) "" else gs_error)
}
write.csv(dplyr::bind_rows(geneset_rows), file.path(TABLE_DIR, "geneset_download_status.csv"), row.names = FALSE)

#Module 06. 새로 받은 CEL에서 RMA와 gene-symbol 행렬 생성 ----
# 기존 CSV를 입력으로 삼지 않습니다. 이 버전이 원자료 지문과 함께 만든 checkpoint만 재사용합니다.
cohort_rows <- list()
annotation_db <- get("hgu133plus2.db", envir = asNamespace("hgu133plus2.db"))
for (gse_id in gse_list) {
  i_gse <- match(gse_id, manifest$gse_id)
  prep_dir <- file.path(INPUT_GEO_PREP_DIR, gse_id); raw_dir <- file.path(INPUT_GEO_RAW_DIR, gse_id)
  dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)
  ready_file <- file.path(READY_DIR, paste0(gse_id, "_ready.rds"))
  checkpoint_file <- file.path(prep_dir, "RMA_checkpoint_v2_1_0.rds")
  n_genes <- 0L; n_samples <- 0L; stage <- "download_completeness"
  cohort_error <- tryCatch({
    if (manifest$metadata_status[i_gse] != "ready" || manifest$raw_status[i_gse] != "ready") {
      stop("온라인 입력 미완료. metadata=", manifest$metadata_status[i_gse], "; RAW=", manifest$raw_status[i_gse])
    }
    stage <- "CEL_extraction"
    raw_tar <- manifest$raw_destination[i_gse]; raw_md5 <- unname(tools::md5sum(raw_tar))
    pheno <- geo_pheno_list[[gse_id]]; cel_members <- raw_members_list[[gse_id]]
    expected_gsm <- sort(pheno$geo_accession)
    if (any(!grepl("GSM[0-9]+", toupper(basename(cel_members))))) stop("CEL 파일명에 GSM이 없습니다. 잘못된 archive인지 확인하십시오.")
    cel_gsm <- regmatches(toupper(basename(cel_members)), regexpr("GSM[0-9]+", toupper(basename(cel_members))))
    if (anyDuplicated(cel_gsm)) stop("archive에 중복 GSM의 CEL이 있습니다. 임의로 대표 파일을 선택하지 않습니다.")
    if (!setequal(cel_gsm, expected_gsm)) {
      write.csv(data.frame(geo_accession = setdiff(expected_gsm, cel_gsm)), file.path(prep_dir, "metadata_samples_without_CEL.csv"), row.names = FALSE)
      write.csv(data.frame(geo_accession = setdiff(cel_gsm, expected_gsm)), file.path(prep_dir, "CEL_samples_without_metadata.csv"), row.names = FALSE)
      stop("GEO metadata와 RAW CEL 표본 집합이 다릅니다.")
    }
    mat <- NULL
    if (file.exists(checkpoint_file)) {
      checkpoint <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
      if (is.list(checkpoint) && identical(checkpoint$version, PIPELINE_VERSION) && identical(checkpoint$raw_md5, raw_md5) &&
          identical(checkpoint$annotation_version, as.character(utils::packageVersion("hgu133plus2.db"))) &&
          identical(checkpoint$affy_version, as.character(utils::packageVersion("affy"))) &&
          identical(checkpoint$orgdb_version, as.character(utils::packageVersion("org.Hs.eg.db"))) &&
          identical(checkpoint$cdf_version, as.character(utils::packageVersion("hgu133plus2cdf")))) {
        mat <- checkpoint$expression
        cat(gse_id, ": 원자료가 일치하는 이번 온라인 실행의 RMA checkpoint 재사용\n")
      }
    }
    if (is.null(mat)) {
      extract_dir <- file.path(raw_dir, "CEL_extracted")
      dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
      if (utils::untar(raw_tar, files = cel_members, exdir = extract_dir) != 0L) stop("CEL 압축 해제 실패")
      cel_paths <- file.path(extract_dir, cel_members)
      if (!all(file.exists(cel_paths)) || any(file.info(cel_paths)$size <= 0)) stop("CEL 추출 파일 누락/빈 파일")
      cels <- sub("\\.gz$", "", cel_paths, ignore.case = TRUE)
      for (i in seq_along(cel_paths)) if (grepl("\\.gz$", cel_paths[i], ignore.case = TRUE)) {
        # R.utils can leave this temporary destination after an interrupted run.
        stale_cel_tmp <- paste0(cels[i], ".tmp")
        if (file.exists(stale_cel_tmp) && unlink(stale_cel_tmp) != 0L) stop("Cannot clear interrupted CEL temporary file: ", stale_cel_tmp)
        R.utils::gunzip(cel_paths[i], destname = cels[i], remove = FALSE, overwrite = TRUE)
      }
      # GSM 정렬을 고정하여 재실행 때 표본 순서가 바뀌지 않도록 합니다.
      order_cel <- order(cel_gsm); cels <- cels[order_cel]; cel_gsm <- cel_gsm[order_cel]
      stage <- "RMA"; cat("RMA: ", gse_id, " / ", length(cels), " CEL files\n", sep = "")
      chip_names <- vapply(cels, function(f) affyio::read.celfile.header(f)$cdfName, character(1))
      chip_names <- tolower(gsub("[^A-Za-z0-9]", "", chip_names))
      if (any(chip_names != "hgu133plus2")) stop("Unexpected physical chip; standard hgu133plus2 CDF is required.")
      # Same background correction, quantile normalization and median polish as affy::rma.
      # justRMA reads PM intensities directly and reduces peak memory use.
      cel_relative <- substring(cels, nchar(extract_dir) + 2L)
      eset <- affy::justRMA(filenames = cel_relative, celfile.path = extract_dir, verbose = TRUE)
      ex <- Biobase::exprs(eset)
      # RMA가 보존한 파일명으로 expression의 실제 열과 GSM 순서를 다시 확인합니다.
      expression_gsm <- regmatches(toupper(colnames(ex)), regexpr("GSM[0-9]+", toupper(colnames(ex))))
      if (length(expression_gsm) != ncol(ex) || anyDuplicated(expression_gsm) || !setequal(expression_gsm, cel_gsm)) stop("RMA 표본명/GSM 불일치")
      ex <- ex[, match(cel_gsm, expression_gsm), drop = FALSE]
      probe_ids <- rownames(ex)
      anno <- AnnotationDbi::select(annotation_db, keys = probe_ids, columns = "SYMBOL", keytype = "PROBEID")
      symbols <- anno$SYMBOL[match(probe_ids, anno$PROBEID)]
      keep <- !is.na(symbols) & nzchar(symbols)
      # 공개 코드의 첫 SYMBOL 매핑 + 동일 symbol probe 평균 정책을 유지합니다.
      agg <- stats::aggregate(ex[keep, , drop = FALSE], by = list(Gene_Symbol = symbols[keep]), FUN = mean)
      mat <- as.matrix(agg[, -1L, drop = FALSE]); rownames(mat) <- agg$Gene_Symbol
      colnames(mat) <- cel_gsm; storage.mode(mat) <- "double"
      saveRDS(list(version = PIPELINE_VERSION, raw_md5 = raw_md5, expression = mat,
        orgdb_version = as.character(utils::packageVersion("org.Hs.eg.db")),
        affy_version = as.character(utils::packageVersion("affy")), annotation_version = as.character(utils::packageVersion("hgu133plus2.db")),
        cdf_version = as.character(utils::packageVersion("hgu133plus2cdf"))), checkpoint_file, compress = FALSE)
      write.csv(anno, file.path(prep_dir, "probe_symbol_mapping.csv"), row.names = FALSE)
      write.csv(data.frame(GSE_ID = gse_id, method = "affy::justRMA (RMA); standard hgu133plus2; first SYMBOL per probe; mean across probes",
        raw_url = manifest$raw_url[i_gse], raw_md5 = raw_md5, n_probes = nrow(ex), n_mapped_probes = sum(keep),
        n_genes = nrow(mat), n_samples = ncol(mat),
        orgdb_version = as.character(utils::packageVersion("org.Hs.eg.db")),
        affy_version = as.character(utils::packageVersion("affy")), annotation_version = as.character(utils::packageVersion("hgu133plus2.db")),
        cdf_version = as.character(utils::packageVersion("hgu133plus2cdf"))), file.path(prep_dir, "preprocessing_provenance.csv"), row.names = FALSE)
      rm(eset, ex, agg); gc()
    }
    stage <- "matrix_and_ID_validation"
    if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) < 1000L || any(!is.finite(mat)) ||
        anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat)) || !setequal(colnames(mat), expected_gsm)) stop("RMA 행렬 크기/숫자/NA/중복/GSM 검사 실패")
    mat <- mat[, expected_gsm, drop = FALSE]
    n_genes <- nrow(mat); n_samples <- ncol(mat)
    meta <- pheno[match(colnames(mat), pheno$geo_accession), , drop = FALSE]
    meta$sampleid <- paste(gse_id, meta$geo_accession, sep = "_"); meta$gse_id <- gse_id
    meta$expr_col_original <- colnames(mat); meta$id_match_method <- "fresh_CEL_GSM_exact"
    colnames(mat) <- meta$sampleid
    stopifnot(identical(colnames(mat), meta$sampleid), identical(meta$geo_accession, expected_gsm), !anyNA(meta$sampleid))
    write.csv(mat, manifest$expression_csv[i_gse])
    write.csv(meta, file.path(prep_dir, paste0(gse_id, "_SampleMetadata.csv")), row.names = FALSE)
    input_paths <- c(raw_tar, manifest$metadata_source_file[i_gse], manifest$expression_csv[i_gse])
    ready <- list(pipeline_version = PIPELINE_VERSION, gse_id = gse_id, expression = mat, metadata = meta,
      input_paths = input_paths, input_md5 = unname(tools::md5sum(input_paths)), loaded_from = "online_RAW_CEL_RMA")
    reuse_ready <- FALSE
    if (file.exists(ready_file)) {
      old_ready <- tryCatch(readRDS(ready_file), error = function(e) NULL)
      reuse_ready <- identical(old_ready, ready)
    }
    if (!reuse_ready) {
      temp_rds <- paste0(ready_file, ".part"); saveRDS(ready, temp_rds, compress = FALSE)
      if (!file.copy(temp_rds, ready_file, overwrite = TRUE)) stop("준비 입력 저장 실패")
      unlink(temp_rds)
    }
    NULL
  }, error = function(e) conditionMessage(e))
  cohort_rows[[gse_id]] <- data.frame(gse_id = gse_id, status = if (is.null(cohort_error)) "ready" else "pending",
    stage = stage, n_genes = n_genes, n_samples = n_samples, ready_file = ready_file,
    detail = if (is.null(cohort_error)) "" else cohort_error)
  cat(gse_id, ": ", cohort_rows[[gse_id]]$status, " | ", n_genes, " genes x ", n_samples, " samples\n", sep = "")
  if (!is.null(cohort_error)) message(cohort_error)
  cat(format(Sys.time()), gse_id, cohort_rows[[gse_id]]$status, stage, cohort_rows[[gse_id]]$detail, "\n", file = LOG_FILE, append = TRUE)
  write.csv(dplyr::bind_rows(cohort_rows), file.path(TABLE_DIR, "input_readiness.csv"), row.names = FALSE)
}

#Module 07. 준비 상태와 다음 실행 단계 ----
readiness <- dplyr::bind_rows(cohort_rows)
pending <- readiness[readiness$status != "ready", , drop = FALSE]
pending <- dplyr::left_join(pending, manifest[, c("gse_id", "metadata_error", "raw_error", "metadata_url", "raw_url")], by = "gse_id")
write.csv(pending, file.path(TABLE_DIR, "pending_inputs.csv"), row.names = FALSE)
run_state$inputs_ready <- nrow(pending) == 0L && nrow(readiness) == length(gse_list) &&
  is.null(tcga_error) && all(vapply(geneset_rows, function(x) identical(x$status, "ready"), logical(1)))
run_state$ready_files <- readiness$ready_file
run_state$ready_md5 <- if (run_state$inputs_ready) unname(tools::md5sum(readiness$ready_file)) else character()
run_state$tcga_files <- tcga_files
run_state$tcga_md5 <- if (is.null(tcga_error)) unname(tools::md5sum(tcga_files)) else character()
run_state$input_checked_at <- as.character(Sys.time())
saveRDS(run_state, STATE_FILE)
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "01_sessionInfo.txt"))
if (!run_state$inputs_ready) {
  if (!is.null(tcga_error)) message("TCGA: ", tcga_error)
  print(dplyr::bind_rows(geneset_rows))
  print(pending[, c("gse_id", "stage", "detail", "metadata_error", "raw_error")], row.names = FALSE)
  stop("온라인 입력 준비가 끝나지 않았습니다. 위의 실제 실패 원인과 tables/pending_inputs.csv를 확인하십시오. 이 버전에서 받은 정상 파일과 완료된 RMA는 다음 실행에 이어서 사용합니다.", call. = FALSE)
}
cat("GEO 9개, TCGA, gene set 준비 완료. 05_Online_analysis.R을 실행하십시오.\n")
