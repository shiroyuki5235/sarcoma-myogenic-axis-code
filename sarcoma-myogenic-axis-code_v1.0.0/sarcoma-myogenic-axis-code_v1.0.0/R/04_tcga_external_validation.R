rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(123)

PROJECT_DIR <- getwd()
INPUT_DIR <- file.path(PROJECT_DIR, "Input")
INPUT_TCGA_DIR <- file.path(INPUT_DIR, "TCGA")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Output")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
TCGA_TAB_DIR <- file.path(OUTPUT_DIR, "tables_validation")
TCGA_PLOT_PNG_DIR <- file.path(OUTPUT_DIR, "plots_png_validation")
TCGA_PLOT_PDF_DIR <- file.path(OUTPUT_DIR, "plots_pdf_validation")
LOG_DIR <- file.path(OUTPUT_DIR, "logs")

dir.create(INPUT_TCGA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TCGA_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TCGA_PLOT_PNG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TCGA_PLOT_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "tcga_external_validation.log")
PLOT_DPI <- 600L

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

data.table::setDTthreads(0L)

TCGA_TAR <- file.path(INPUT_TCGA_DIR, "sarc_tcga_pub.tar.gz")
TCGA_URL <- "https://datahub.assets.cbioportal.org/sarc_tcga_pub.tar.gz"
TCGA_DIR <- file.path(INPUT_TCGA_DIR, "sarc_tcga_pub")

if (!file.exists(TCGA_TAR)) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Downloading TCGA-SARC tar.gz")
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
  download.file(TCGA_URL, TCGA_TAR, mode = "wb", quiet = FALSE)
}

if (!dir.exists(TCGA_DIR)) {
  dir.create(TCGA_DIR, recursive = TRUE, showWarnings = FALSE)
  untar(TCGA_TAR, exdir = INPUT_TCGA_DIR)
}

expr_file <- file.path(TCGA_DIR, "data_mrna_seq_v2_rsem.txt")
clin_sample_file <- file.path(TCGA_DIR, "data_clinical_sample.txt")
clin_patient_file <- file.path(TCGA_DIR, "data_clinical_patient.txt")
if (!file.exists(expr_file)) stop("TCGA expression file not found after untar")
if (!file.exists(clin_sample_file)) stop("TCGA sample clinical file not found after untar")
if (!file.exists(clin_patient_file)) stop("TCGA patient clinical file not found after untar")

expr_dt <- data.table::fread(expr_file, data.table = FALSE)
clin_sample <- data.table::fread(clin_sample_file, skip = 4L, data.table = FALSE)
clin_patient <- data.table::fread(clin_patient_file, skip = 4L, data.table = FALSE)

if (!("Hugo_Symbol" %in% names(expr_dt))) names(expr_dt)[1] <- "Hugo_Symbol"
expr_dt$Hugo_Symbol <- as.character(expr_dt$Hugo_Symbol)
expr_dt <- expr_dt[!is.na(expr_dt$Hugo_Symbol) & expr_dt$Hugo_Symbol != "", , drop = FALSE]
expr_dt <- expr_dt[!duplicated(expr_dt$Hugo_Symbol), , drop = FALSE]

sample_cols <- setdiff(names(expr_dt), c("Hugo_Symbol", "Entrez_Gene_Id"))
expr_mat <- as.matrix(expr_dt[, sample_cols, drop = FALSE])
rownames(expr_mat) <- expr_dt$Hugo_Symbol
storage.mode(expr_mat) <- "double"
expr_mat <- log2(expr_mat + 1)

if (!("PATIENT_ID" %in% names(clin_sample))) stop("PATIENT_ID missing in TCGA sample clinical")
if (!("SAMPLE_ID" %in% names(clin_sample))) stop("SAMPLE_ID missing in TCGA sample clinical")
if (!("PATIENT_ID" %in% names(clin_patient))) stop("PATIENT_ID missing in TCGA patient clinical")

tcga_meta <- merge(clin_sample, clin_patient, by = "PATIENT_ID", all.x = TRUE, suffixes = c("_sample", "_patient"))
names(tcga_meta) <- make.unique(names(tcga_meta), sep = "__dup")

need_meta <- c("CANCER_TYPE_DETAILED", "ONCOTREE_CODE", "SHORT_HISTOLOGY", "TUMOR_TISSUE_SITE")
for (nm in need_meta) if (!(nm %in% names(tcga_meta))) tcga_meta[[nm]] <- NA_character_

tcga_meta$histology_label <- tcga_meta$CANCER_TYPE_DETAILED
tcga_meta$histology_label[is.na(tcga_meta$histology_label) | tcga_meta$histology_label == ""] <- tcga_meta$SHORT_HISTOLOGY[is.na(tcga_meta$histology_label) | tcga_meta$histology_label == ""]
tcga_meta$histology_label[is.na(tcga_meta$histology_label) | tcga_meta$histology_label == ""] <- tcga_meta$ONCOTREE_CODE[is.na(tcga_meta$histology_label) | tcga_meta$histology_label == ""]

tcga_meta$histology_group <- NA_character_
txt <- tolower(tcga_meta$histology_label)
site_txt <- tolower(tcga_meta$TUMOR_TISSUE_SITE)
tcga_meta$histology_group[grepl("dedifferentiated liposarcoma|\\bddls\\b", txt)] <- "DDLPS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("undifferentiated pleomorphic sarcoma|malignant fibrous histiocytoma|high-grade spindle", txt)] <- "UPS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("myxofibrosarcoma", txt)] <- "MFS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("synovial sarcoma", txt)] <- "SS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("malignant peripheral nerve sheath|mpnst", txt)] <- "MPNST"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("uterine leiomyosarcoma|uterine leiomyosarcome", txt)] <- "ULMS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("leiomyosarcoma|leiomyosarcome", txt)] <- "STLMS"
tcga_meta$histology_group[is.na(tcga_meta$histology_group) & grepl("uter", site_txt) & grepl("leiomyo|\\blms\\b", txt)] <- "ULMS"

keep_expr_ids <- intersect(colnames(expr_mat), tcga_meta$SAMPLE_ID)
expr_mat <- expr_mat[, keep_expr_ids, drop = FALSE]
tcga_meta <- tcga_meta[match(colnames(expr_mat), tcga_meta$SAMPLE_ID), , drop = FALSE]

sample_master2 <- read.csv(file.path(TABLE_DIR, "sample_master_metadata_curated_with_cluster.csv"), stringsAsFactors = FALSE, check.names = FALSE)
cluster_final <- read.csv(file.path(TABLE_DIR, "cluster_final_k_selected.csv"), stringsAsFactors = FALSE, check.names = FALSE)
K_FINAL <- cluster_final$K_FINAL[1]
cluster_sizes <- read.csv(file.path(TABLE_DIR, paste0("cluster_sizes_k", K_FINAL, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
exprs_combat <- readRDS(file.path(OUTPUT_DIR, "exprs_combat_curated.rds"))

minor_cluster <- cluster_sizes$Cluster[which.min(cluster_sizes$N)]
major_cluster <- cluster_sizes$Cluster[which.max(cluster_sizes$N)]

sig_file1 <- file.path(OUTPUT_DIR, "tables_addon", "additional_minor_cluster_myogenic_signature_genes.csv")
sig_file2 <- file.path(TABLE_DIR, paste0("DEG_top_up_", minor_cluster, "_vsRest.csv"))

if (file.exists(sig_file1)) {
  sig_df <- read.csv(sig_file1, stringsAsFactors = FALSE, check.names = FALSE)
  minor_genes <- unique(sig_df$GeneSymbol)
} else {
  sig_df <- read.csv(sig_file2, stringsAsFactors = FALSE, check.names = FALSE)
  minor_genes <- unique(sig_df$GeneSymbol)
}
minor_genes <- minor_genes[!is.na(minor_genes) & minor_genes != ""]
minor_genes <- intersect(minor_genes, rownames(exprs_combat))
minor_genes <- intersect(minor_genes, rownames(expr_mat))
minor_genes <- head(minor_genes, 30)
if (length(minor_genes) < 5) stop("Too few overlapping minor signature genes between discovery and TCGA")

minor_samples <- sample_master2$sampleid[sample_master2$Cluster == minor_cluster]
major_samples <- sample_master2$sampleid[sample_master2$Cluster == major_cluster]

common_centroid_genes <- intersect(minor_genes, rownames(exprs_combat))
common_centroid_genes <- intersect(common_centroid_genes, rownames(expr_mat))
minor_centroid <- rowMeans(exprs_combat[common_centroid_genes, minor_samples, drop = FALSE])
major_centroid <- rowMeans(exprs_combat[common_centroid_genes, major_samples, drop = FALSE])

expr_tcga_sig <- expr_mat[common_centroid_genes, , drop = FALSE]
expr_tcga_z <- t(scale(t(expr_tcga_sig)))
expr_tcga_z[!is.finite(expr_tcga_z)] <- 0
minor_score <- colMeans(expr_tcga_z, na.rm = TRUE)
minor_cor <- rep(NA_real_, ncol(expr_tcga_sig))
major_cor <- rep(NA_real_, ncol(expr_tcga_sig))
for (i in seq_len(ncol(expr_tcga_sig))) {
  minor_cor[i] <- suppressWarnings(cor(as.numeric(expr_tcga_sig[, i]), as.numeric(minor_centroid), method = "pearson", use = "pairwise.complete.obs"))
  major_cor[i] <- suppressWarnings(cor(as.numeric(expr_tcga_sig[, i]), as.numeric(major_centroid), method = "pearson", use = "pairwise.complete.obs"))
}

val_df <- data.frame(
  SAMPLE_ID = colnames(expr_tcga_sig),
  minor_cluster_score = as.numeric(minor_score),
  cor_minor_centroid = minor_cor,
  cor_major_centroid = major_cor,
  predicted_axis = ifelse(minor_cor >= major_cor, "minor_like", "major_like"),
  stringsAsFactors = FALSE
)
meta_add <- tcga_meta[match(val_df$SAMPLE_ID, tcga_meta$SAMPLE_ID), , drop = FALSE]
meta_add <- meta_add[, setdiff(names(meta_add), names(val_df)), drop = FALSE]
val_df <- cbind(val_df, meta_add)
val_df$is_STLMS <- ifelse(val_df$histology_group == "STLMS", 1L, 0L)
val_df$is_ULMS <- ifelse(val_df$histology_group == "ULMS", 1L, 0L)
val_df$is_MPNST <- ifelse(val_df$histology_group == "MPNST", 1L, 0L)

write.csv(val_df, file.path(TCGA_TAB_DIR, "tcga_validation_per_sample_full.csv"), row.names = FALSE)
write.csv(data.frame(GeneSymbol = common_centroid_genes), file.path(TCGA_TAB_DIR, "tcga_validation_signature_genes_used.csv"), row.names = FALSE)

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
  wt <- try(wilcox.test(pos, neg), silent = TRUE)
  if (inherits(wt, "try-error")) wt <- NULL

  auc_value <- NA_real_
  if (length(pos) > 0 && length(neg) > 0) {
    all_scores <- c(pos, neg)
    all_labels <- c(rep(1, length(pos)), rep(0, length(neg)))
    r <- rank(all_scores, ties.method = "average")
    u <- sum(r[all_labels == 1]) - length(pos) * (length(pos) + 1) / 2
    auc_value <- u / (length(pos) * length(neg))
  }

  pred_tab <- table(factor(sub_df$is_target_lms, levels = c(0, 1)), factor(sub_df$predicted_axis == "minor_like", levels = c(FALSE, TRUE)))
  fisher_res <- try(fisher.test(pred_tab), silent = TRUE)
  if (inherits(fisher_res, "try-error")) fisher_res <- NULL

  scenario_summary[[scenario_nm]] <- data.frame(
    scenario = scenario_nm,
    n_samples = nrow(sub_df),
    n_target_lms = sum(sub_df$is_target_lms == 1),
    n_other = sum(sub_df$is_target_lms == 0),
    mean_score_target = if (length(pos) > 0) mean(pos, na.rm = TRUE) else NA_real_,
    mean_score_other = if (length(neg) > 0) mean(neg, na.rm = TRUE) else NA_real_,
    wilcox_p = if (is.null(wt)) NA_real_ else wt$p.value,
    auc = auc_value,
    fisher_or = if (is.null(fisher_res)) NA_real_ else unname(fisher_res$estimate),
    fisher_p = if (is.null(fisher_res)) NA_real_ else fisher_res$p.value,
    stringsAsFactors = FALSE
  )

  hist_sum <- sub_df %>% group_by(histology_group) %>% summarise(n = n(), mean_score = mean(minor_cluster_score, na.rm = TRUE), sd_score = sd(minor_cluster_score, na.rm = TRUE), mean_cor_minor = mean(cor_minor_centroid, na.rm = TRUE), prop_minor_like = mean(predicted_axis == "minor_like", na.rm = TRUE), .groups = "drop") %>% arrange(desc(mean_score), desc(n))
  write.csv(hist_sum, file.path(TCGA_TAB_DIR, paste0("tcga_histology_summary_", scenario_nm, ".csv")), row.names = FALSE)

  p1 <- ggplot(sub_df, aes(histology_group, minor_cluster_score, fill = histology_group)) + geom_boxplot(outlier.size = 0.3) + coord_flip() + theme_bw() + guides(fill = "none") + labs(title = paste0("TCGA validation score by histology: ", scenario_nm), x = NULL, y = "Minor-cluster myogenic signature score")
  ggsave(file.path(TCGA_PLOT_PNG_DIR, paste0("tcga_score_by_histology_", scenario_nm, ".png")), plot = p1, width = 8, height = 5, dpi = PLOT_DPI)
  ggsave(file.path(TCGA_PLOT_PDF_DIR, paste0("tcga_score_by_histology_", scenario_nm, ".pdf")), plot = p1, width = 8, height = 5, device = "pdf")

  p2 <- ggplot(sub_df, aes(histology_group, cor_minor_centroid, fill = histology_group)) + geom_boxplot(outlier.size = 0.3) + coord_flip() + theme_bw() + guides(fill = "none") + labs(title = paste0("TCGA centroid correlation by histology: ", scenario_nm), x = NULL, y = "Correlation to discovery minor centroid")
  ggsave(file.path(TCGA_PLOT_PNG_DIR, paste0("tcga_cor_minor_by_histology_", scenario_nm, ".png")), plot = p2, width = 8, height = 5, dpi = PLOT_DPI)
  ggsave(file.path(TCGA_PLOT_PDF_DIR, paste0("tcga_cor_minor_by_histology_", scenario_nm, ".pdf")), plot = p2, width = 8, height = 5, device = "pdf")
}

if (length(scenario_summary) > 0) write.csv(bind_rows(scenario_summary), file.path(TCGA_TAB_DIR, "tcga_validation_scenario_summary.csv"), row.names = FALSE)

marker_genes <- intersect(c("MYH11", "LMOD1", "MYOCD", "MYLK", "CNN1"), rownames(expr_mat))
if (length(marker_genes) > 0) {
  marker_df <- as.data.frame(t(expr_mat[marker_genes, , drop = FALSE]))
  marker_df$SAMPLE_ID <- rownames(marker_df)
  marker_df <- marker_df %>% left_join(val_df[, c("SAMPLE_ID", "histology_group")], by = "SAMPLE_ID")
  marker_long <- pivot_longer(marker_df, cols = all_of(marker_genes), names_to = "GeneSymbol", values_to = "Expression")
  write.csv(marker_long, file.path(TCGA_TAB_DIR, "tcga_validation_marker_long.csv"), row.names = FALSE)
  p_marker <- ggplot(marker_long[marker_long$histology_group %in% c("STLMS", "ULMS", "DDLPS", "UPS", "MFS", "SS", "MPNST"), , drop = FALSE], aes(histology_group, Expression, fill = histology_group)) + geom_boxplot(outlier.size = 0.2) + coord_flip() + facet_wrap(~ GeneSymbol, scales = "free_x") + theme_bw() + guides(fill = "none") + labs(title = "Representative myogenic markers in TCGA-SARC", x = NULL, y = "log2(RSEM + 1)")
  ggsave(file.path(TCGA_PLOT_PNG_DIR, "tcga_representative_myogenic_markers.png"), plot = p_marker, width = 11, height = 7, dpi = PLOT_DPI)
  ggsave(file.path(TCGA_PLOT_PDF_DIR, "tcga_representative_myogenic_markers.pdf"), plot = p_marker, width = 11, height = 7, device = "pdf")
}

msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "03_tcga_external_validation finished")
cat(msg, "\n")
cat(msg, "\n", file = LOG_FILE, append = TRUE)
