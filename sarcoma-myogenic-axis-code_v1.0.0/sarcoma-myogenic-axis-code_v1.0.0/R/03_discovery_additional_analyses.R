rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(123)

PROJECT_DIR <- getwd()
OUTPUT_DIR <- file.path(PROJECT_DIR, "Output")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
ADDON_TAB_DIR <- file.path(OUTPUT_DIR, "tables_addon")
ADDON_PLOT_PNG_DIR <- file.path(OUTPUT_DIR, "plots_png_addon")
ADDON_PLOT_PDF_DIR <- file.path(OUTPUT_DIR, "plots_pdf_addon")
LOG_DIR <- file.path(OUTPUT_DIR, "logs")

dir.create(ADDON_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ADDON_PLOT_PNG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ADDON_PLOT_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "discovery_additional_analyses.log")
PLOT_DPI <- 600L

library(matrixStats)
library(cluster)
library(ggplot2)
library(dplyr)
library(tidyr)
library(limma)
library(sva)
library(igraph)
library(ggraph)

RUN_NETWORK <- TRUE
NET_OVERALL_THR <- 0.70
NET_RESIDUAL_THR <- 0.30
K_GRID_ADDON <- 2:6
BOOT_B_ADDON <- 10L
BOOT_PROP_ADDON <- 0.80
KSUPPORT_TOP_VAR_ADDON <- 3000L
KSUPPORT_N_PC_ADDON <- 30L
LEAVE_ONE_OUT_MIN_SAMPLES <- 100L
TOP_MARKER_N_FOR_SCORE <- 30L
PLOT_DPI <- 600L

exprs_all <- readRDS(file.path(OUTPUT_DIR, "exprs_all_curated_prefiltered.rds"))
exprs_combat <- readRDS(file.path(OUTPUT_DIR, "exprs_combat_curated.rds"))
sample_master2 <- read.csv(file.path(TABLE_DIR, "sample_master_metadata_curated_with_cluster.csv"), stringsAsFactors = FALSE, check.names = FALSE)
cluster_final <- read.csv(file.path(TABLE_DIR, "cluster_final_k_selected.csv"), stringsAsFactors = FALSE, check.names = FALSE)
K_FINAL <- cluster_final$K_FINAL[1]
cluster_df <- read.csv(file.path(TABLE_DIR, paste0("cluster_assignments_k", K_FINAL, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
cluster_sizes <- read.csv(file.path(TABLE_DIR, paste0("cluster_sizes_k", K_FINAL, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)

minor_cluster <- cluster_sizes$Cluster[which.min(cluster_sizes$N)]
major_cluster <- cluster_sizes$Cluster[which.max(cluster_sizes$N)]

sample_master2$is_minor_cluster <- ifelse(sample_master2$Cluster == minor_cluster, 1L, 0L)
sample_master2$is_lms <- ifelse(sample_master2$diagnosis_group == "Leiomyosarcoma", 1L, 0L)

cohort_summary <- sample_master2 %>% count(gse_id, diagnosis_group, name = "n") %>% arrange(desc(n), gse_id, diagnosis_group)
write.csv(cohort_summary, file.path(ADDON_TAB_DIR, "additional_cohort_summary_by_gse_and_diagnosis.csv"), row.names = FALSE)

cluster_diag_summary <- sample_master2 %>% count(Cluster, diagnosis_group, name = "n") %>% group_by(Cluster) %>% mutate(prop = n / sum(n)) %>% ungroup() %>% arrange(Cluster, desc(prop), desc(n))
write.csv(cluster_diag_summary, file.path(ADDON_TAB_DIR, "additional_cluster_diagnosis_summary.csv"), row.names = FALSE)

comp_effects <- list()
for (v in c("gse_id", "diagnosis_group", "tissue_group")) {
  tmp <- sample_master2[, c("Cluster", v), drop = FALSE]
  tmp <- tmp[!is.na(tmp[[v]]) & tmp[[v]] != "" & !is.na(tmp$Cluster), , drop = FALSE]
  if (nrow(tmp) > 0 && length(unique(tmp[[v]])) >= 2) {
    tab <- table(tmp$Cluster, tmp[[v]])
    chisq_obj <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 10000))
    x_cv <- suppressWarnings(chisq.test(tab, correct = FALSE))
    n_cv <- sum(tab)
    k_cv <- min(nrow(tab) - 1, ncol(tab) - 1)
    cramers_v_value <- if (n_cv <= 0 || k_cv <= 0) NA_real_ else sqrt(as.numeric(x_cv$statistic) / (n_cv * k_cv))
    comp_effects[[v]] <- data.frame(variable = v, n_samples_used = nrow(tmp), n_categories = ncol(tab), p_value = chisq_obj$p.value, cramers_v = cramers_v_value, stringsAsFactors = FALSE)
  }
}
if (length(comp_effects) > 0) write.csv(bind_rows(comp_effects), file.path(ADDON_TAB_DIR, "additional_composition_effect_sizes.csv"), row.names = FALSE)

mh_df <- sample_master2 %>% filter(!is.na(gse_id), gse_id != "", !is.na(is_lms), !is.na(is_minor_cluster)) %>% mutate(gse_id = as.character(gse_id))
mh_keep_gse <- mh_df %>% count(gse_id, is_lms, is_minor_cluster, name = "n") %>% group_by(gse_id) %>% summarise(n_total = sum(n), .groups = "drop") %>% filter(n_total >= 20)
mh_df2 <- mh_df %>% filter(gse_id %in% mh_keep_gse$gse_id)
if (nrow(mh_df2) > 0 && length(unique(mh_df2$gse_id)) >= 2) {
  strata <- sort(unique(mh_df2$gse_id))
  arr <- array(0, dim = c(2, 2, length(strata)), dimnames = list(LMS = c("No", "Yes"), MinorCluster = c("No", "Yes"), GSE = strata))
  for (g in strata) {
    sub <- mh_df2[mh_df2$gse_id == g, , drop = FALSE]
    tab <- table(factor(sub$is_lms, levels = c(0, 1)), factor(sub$is_minor_cluster, levels = c(0, 1)))
    arr[, , g] <- tab
  }
  mh <- mantelhaen.test(arr)
  mh_result_df <- data.frame(test = "Mantel-Haenszel", strata_n = length(strata), common_odds_ratio = unname(mh$estimate), p_value = mh$p.value, stringsAsFactors = FALSE)
  write.csv(mh_result_df, file.path(ADDON_TAB_DIR, "additional_mantel_haenszel_lms_vs_minor_cluster.csv"), row.names = FALSE)
}

minor_deg_file <- file.path(TABLE_DIR, paste0("DEG_top_up_", minor_cluster, "_vsRest.csv"))
major_deg_file <- file.path(TABLE_DIR, paste0("DEG_top_up_", major_cluster, "_vsRest.csv"))
if (file.exists(minor_deg_file)) {
  minor_markers_df <- read.csv(minor_deg_file, stringsAsFactors = FALSE, check.names = FALSE)
  minor_markers <- unique(minor_markers_df$GeneSymbol)
  minor_markers <- minor_markers[!is.na(minor_markers) & minor_markers != ""]
  minor_markers <- intersect(minor_markers, rownames(exprs_combat))
  minor_markers <- head(minor_markers, TOP_MARKER_N_FOR_SCORE)
  if (length(minor_markers) > 0) {
    expr_sig2 <- exprs_combat[minor_markers, , drop = FALSE]
    expr_sig2_z <- t(scale(t(expr_sig2)))
    expr_sig2_z[!is.finite(expr_sig2_z)] <- 0
    minor_score <- colMeans(expr_sig2_z, na.rm = TRUE)
    score_df <- data.frame(sampleid = colnames(exprs_combat), minor_cluster_score = as.numeric(minor_score), stringsAsFactors = FALSE) %>% left_join(as.data.frame(sample_master2), by = "sampleid")
    write.csv(score_df, file.path(ADDON_TAB_DIR, "additional_minor_cluster_myogenic_signature_score_per_sample.csv"), row.names = FALSE)
    write.csv(data.frame(GeneSymbol = minor_markers), file.path(ADDON_TAB_DIR, "additional_minor_cluster_myogenic_signature_genes.csv"), row.names = FALSE)
    score_summary_cluster <- score_df %>% group_by(Cluster) %>% summarise(n = n(), mean_score = mean(minor_cluster_score, na.rm = TRUE), sd_score = sd(minor_cluster_score, na.rm = TRUE), .groups = "drop")
    write.csv(score_summary_cluster, file.path(ADDON_TAB_DIR, "additional_minor_cluster_myogenic_score_summary_by_cluster.csv"), row.names = FALSE)
    score_summary_dx <- score_df %>% filter(!is.na(diagnosis_group), diagnosis_group != "") %>% group_by(diagnosis_group) %>% summarise(n = n(), mean_score = mean(minor_cluster_score, na.rm = TRUE), sd_score = sd(minor_cluster_score, na.rm = TRUE), .groups = "drop") %>% arrange(desc(mean_score), desc(n))
    write.csv(score_summary_dx, file.path(ADDON_TAB_DIR, "additional_minor_cluster_myogenic_score_summary_by_diagnosis.csv"), row.names = FALSE)
    lm_df <- score_df %>% filter(!is.na(gse_id), gse_id != "", !is.na(diagnosis_group), diagnosis_group != "")
    lm_df$is_lms_factor <- factor(ifelse(lm_df$diagnosis_group == "Leiomyosarcoma", "LMS", "Other"))
    lm_df$gse_id <- factor(lm_df$gse_id)
    lm_fit <- stats::lm(minor_cluster_score ~ is_lms_factor + gse_id, data = lm_df)
    lm_coef <- as.data.frame(summary(lm_fit)$coefficients)
    lm_coef$term <- rownames(lm_coef)
    rownames(lm_coef) <- NULL
    names(lm_coef) <- c("estimate", "std_error", "t_value", "p_value", "term")
    write.csv(lm_coef, file.path(ADDON_TAB_DIR, "additional_minor_cluster_myogenic_score_lm_isLMS_plus_gse.csv"), row.names = FALSE)
    p_score_cluster <- ggplot(score_df, aes(Cluster, minor_cluster_score, fill = Cluster)) + geom_boxplot(outlier.size = 0.3) + theme_bw() + guides(fill = "none") + labs(title = "Minor-cluster myogenic signature score by cluster", x = "Cluster", y = "Mean z-score across minor-cluster markers")
    ggsave(file.path(ADDON_PLOT_PNG_DIR, "additional_minor_cluster_myogenic_score_by_cluster.png"), plot = p_score_cluster, width = 6, height = 4, dpi = PLOT_DPI)
    ggsave(file.path(ADDON_PLOT_PDF_DIR, "additional_minor_cluster_myogenic_score_by_cluster.pdf"), plot = p_score_cluster, width = 6, height = 4, device = "pdf")
    plot_dx_df <- score_df %>% filter(!is.na(diagnosis_group), diagnosis_group != "") %>% count(diagnosis_group, name = "n_dx") %>% filter(n_dx >= 10) %>% right_join(score_df, by = "diagnosis_group") %>% filter(!is.na(diagnosis_group), diagnosis_group != "")
    if (nrow(plot_dx_df) > 0) {
      dx_order <- plot_dx_df %>% group_by(diagnosis_group) %>% summarise(med = median(minor_cluster_score, na.rm = TRUE), .groups = "drop") %>% arrange(desc(med)) %>% pull(diagnosis_group)
      plot_dx_df$diagnosis_group <- factor(plot_dx_df$diagnosis_group, levels = dx_order)
      p_score_dx <- ggplot(plot_dx_df, aes(diagnosis_group, minor_cluster_score, fill = diagnosis_group)) + geom_boxplot(outlier.size = 0.2) + coord_flip() + theme_bw() + guides(fill = "none") + labs(title = "Minor-cluster myogenic signature score by diagnosis group", x = NULL, y = "Mean z-score across minor-cluster markers")
      ggsave(file.path(ADDON_PLOT_PNG_DIR, "additional_minor_cluster_myogenic_score_by_diagnosis.png"), plot = p_score_dx, width = 8, height = 6, dpi = PLOT_DPI)
      ggsave(file.path(ADDON_PLOT_PDF_DIR, "additional_minor_cluster_myogenic_score_by_diagnosis.pdf"), plot = p_score_dx, width = 8, height = 6, device = "pdf")
    }
  }
}

all_gse_counts <- sample_master2 %>% count(gse_id, name = "n_total") %>% arrange(desc(n_total), gse_id)
loo_candidates <- all_gse_counts$gse_id[all_gse_counts$n_total >= LEAVE_ONE_OUT_MIN_SAMPLES]
loo_results <- list()
loo_metrics_list <- list()

for (g in loo_candidates) {
  keep_ids <- sample_master2$sampleid[sample_master2$gse_id != g]
  keep_batch <- sub("_.*$", "", keep_ids)
  x_sub <- exprs_all[, keep_ids, drop = FALSE]
  if (length(unique(keep_batch)) >= 2) {
    x_sub2 <- try(sva::ComBat(dat = x_sub, batch = keep_batch, par.prior = TRUE, prior.plots = FALSE), silent = TRUE)
    if (inherits(x_sub2, "try-error")) x_sub2 <- NULL
    if (is.null(x_sub2)) x_sub2 <- limma::removeBatchEffect(x_sub, batch = keep_batch)
    x_sub <- x_sub2
  }
  storage.mode(x_sub) <- "double"
  rv_k <- matrixStats::rowVars(x_sub)
  ord_k <- order(rv_k, decreasing = TRUE)
  n_k_features <- min(KSUPPORT_TOP_VAR_ADDON, nrow(x_sub))
  xk <- x_sub[ord_k[seq_len(n_k_features)], , drop = FALSE]
  pc_fit_sub <- prcomp(t(xk), center = TRUE, scale. = FALSE)
  n_pc_use_sub <- min(KSUPPORT_N_PC_ADDON, ncol(pc_fit_sub$x))
  if (n_pc_use_sub < 2L) next
  x_eval <- t(pc_fit_sub$x[, seq_len(n_pc_use_sub), drop = FALSE])
  dist_mat_sub <- dist(t(x_eval))
  hc_sub <- hclust(dist_mat_sub, method = "ward.D2")
  k_metrics_sub <- list()
  for (k in K_GRID_ADDON) {
    cl_k <- cutree(hc_sub, k = k)
    sil <- try(cluster::silhouette(cl_k, dist_mat_sub), silent = TRUE)
    if (inherits(sil, "try-error")) sil <- NULL
    sil_avg <- if (!is.null(sil)) mean(sil[, "sil_width"], na.rm = TRUE) else NA_real_
    set.seed(123 + k)
    n_sub <- max(5L, floor(ncol(x_eval) * BOOT_PROP_ADDON))
    ari_vec <- rep(NA_real_, BOOT_B_ADDON)
    for (b in seq_len(BOOT_B_ADDON)) {
      idx_b <- sort(sample.int(ncol(x_eval), size = n_sub, replace = FALSE))
      sub_expr <- x_eval[, idx_b, drop = FALSE]
      sub_hc <- try(hclust(dist(t(sub_expr)), method = "ward.D2"), silent = TRUE)
      if (!inherits(sub_hc, "try-error") && !is.null(sub_hc)) {
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
          if (denom != 0) ari_vec[b] <- (sum(nij2) - expected) / denom
        }
      }
    }
    k_metrics_sub[[as.character(k)]] <- data.frame(k = k, avg_silhouette = sil_avg, ari_mean = mean(ari_vec, na.rm = TRUE), ari_median = median(ari_vec, na.rm = TRUE), ari_sd = sd(ari_vec, na.rm = TRUE), ari_n_valid = sum(!is.na(ari_vec)), stringsAsFactors = FALSE)
  }
  k_metrics_df_sub <- bind_rows(k_metrics_sub)
  k_metrics_df_sub$rank_silhouette <- rank(-k_metrics_df_sub$avg_silhouette, ties.method = "min")
  k_metrics_df_sub$rank_ari_mean <- rank(-k_metrics_df_sub$ari_mean, ties.method = "min")
  k_metrics_df_sub$rank_sum <- k_metrics_df_sub$rank_silhouette + k_metrics_df_sub$rank_ari_mean
  pick_sub <- k_metrics_df_sub %>% arrange(rank_sum, desc(ari_mean), desc(avg_silhouette), k) %>% slice(1)
  tmp <- k_metrics_df_sub
  tmp$excluded_gse <- g
  loo_metrics_list[[g]] <- tmp
  k2_row <- tmp[tmp$k == 2, , drop = FALSE]
  loo_results[[g]] <- data.frame(excluded_gse = g, n_remaining = length(keep_ids), selected_k = pick_sub$k[1], k2_silhouette = if (nrow(k2_row) == 1) k2_row$avg_silhouette else NA_real_, k2_ari_mean = if (nrow(k2_row) == 1) k2_row$ari_mean else NA_real_, stringsAsFactors = FALSE)
}

if (length(loo_results) > 0) {
  loo_summary <- bind_rows(loo_results)
  write.csv(loo_summary, file.path(ADDON_TAB_DIR, "additional_leave_one_dataset_out_summary.csv"), row.names = FALSE)
  p_loo_k <- ggplot(loo_summary, aes(reorder(excluded_gse, -n_remaining), selected_k)) + geom_point(size = 3) + theme_bw() + labs(title = "Leave-one-dataset-out selected k", x = "Excluded GEO Series", y = "Selected k")
  ggsave(file.path(ADDON_PLOT_PNG_DIR, "additional_leave_one_dataset_out_selected_k.png"), plot = p_loo_k, width = 8, height = 4, dpi = PLOT_DPI)
  ggsave(file.path(ADDON_PLOT_PDF_DIR, "additional_leave_one_dataset_out_selected_k.pdf"), plot = p_loo_k, width = 8, height = 4, device = "pdf")
  p_loo_ari <- ggplot(loo_summary, aes(reorder(excluded_gse, k2_ari_mean), k2_ari_mean)) + geom_col() + coord_flip() + theme_bw() + labs(title = "Leave-one-dataset-out stability of k = 2", x = "Excluded GEO Series", y = "Mean ARI at k = 2")
  ggsave(file.path(ADDON_PLOT_PNG_DIR, "additional_leave_one_dataset_out_k2_ari.png"), plot = p_loo_ari, width = 8, height = 5, dpi = PLOT_DPI)
  ggsave(file.path(ADDON_PLOT_PDF_DIR, "additional_leave_one_dataset_out_k2_ari.pdf"), plot = p_loo_ari, width = 8, height = 5, device = "pdf")
}
if (length(loo_metrics_list) > 0) write.csv(bind_rows(loo_metrics_list), file.path(ADDON_TAB_DIR, "additional_leave_one_dataset_out_kmetrics.csv"), row.names = FALSE)

sens_specs <- list(
  sensitivity_excluding_GSE34800 = c("GSE34800"),
  sensitivity_excluding_GSE94321 = c("GSE94321"),
  sensitivity_excluding_both_uncertain_gse = c("GSE34800", "GSE94321")
)
sens_list <- list()
for (label_stub in names(sens_specs)) {
  excluded_gse <- sens_specs[[label_stub]]
  keep_ids <- sample_master2$sampleid[!sample_master2$gse_id %in% excluded_gse]
  if (length(keep_ids) < 30) next
  keep_batch <- sub("_.*$", "", keep_ids)
  x_sub <- exprs_all[, keep_ids, drop = FALSE]
  if (length(unique(keep_batch)) >= 2) {
    x_sub2 <- try(sva::ComBat(dat = x_sub, batch = keep_batch, par.prior = TRUE, prior.plots = FALSE), silent = TRUE)
    if (inherits(x_sub2, "try-error")) x_sub2 <- NULL
    if (is.null(x_sub2)) x_sub2 <- limma::removeBatchEffect(x_sub, batch = keep_batch)
    x_sub <- x_sub2
  }
  storage.mode(x_sub) <- "double"
  rv_k <- matrixStats::rowVars(x_sub)
  ord_k <- order(rv_k, decreasing = TRUE)
  n_k_features <- min(KSUPPORT_TOP_VAR_ADDON, nrow(x_sub))
  xk <- x_sub[ord_k[seq_len(n_k_features)], , drop = FALSE]
  pc_fit_sub <- prcomp(t(xk), center = TRUE, scale. = FALSE)
  n_pc_use_sub <- min(KSUPPORT_N_PC_ADDON, ncol(pc_fit_sub$x))
  if (n_pc_use_sub < 2L) next
  x_eval <- t(pc_fit_sub$x[, seq_len(n_pc_use_sub), drop = FALSE])
  dist_mat_sub <- dist(t(x_eval))
  hc_sub <- hclust(dist_mat_sub, method = "ward.D2")
  k_metrics_sub <- list()
  for (k in K_GRID_ADDON) {
    cl_k <- cutree(hc_sub, k = k)
    sil <- try(cluster::silhouette(cl_k, dist_mat_sub), silent = TRUE)
    if (inherits(sil, "try-error")) sil <- NULL
    sil_avg <- if (!is.null(sil)) mean(sil[, "sil_width"], na.rm = TRUE) else NA_real_
    set.seed(123 + k)
    n_sub <- max(5L, floor(ncol(x_eval) * BOOT_PROP_ADDON))
    ari_vec <- rep(NA_real_, BOOT_B_ADDON)
    for (b in seq_len(BOOT_B_ADDON)) {
      idx_b <- sort(sample.int(ncol(x_eval), size = n_sub, replace = FALSE))
      sub_expr <- x_eval[, idx_b, drop = FALSE]
      sub_hc <- try(hclust(dist(t(sub_expr)), method = "ward.D2"), silent = TRUE)
      if (!inherits(sub_hc, "try-error") && !is.null(sub_hc)) {
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
          if (denom != 0) ari_vec[b] <- (sum(nij2) - expected) / denom
        }
      }
    }
    k_metrics_sub[[as.character(k)]] <- data.frame(k = k, avg_silhouette = sil_avg, ari_mean = mean(ari_vec, na.rm = TRUE), ari_median = median(ari_vec, na.rm = TRUE), ari_sd = sd(ari_vec, na.rm = TRUE), ari_n_valid = sum(!is.na(ari_vec)), stringsAsFactors = FALSE)
  }
  k_metrics_df_sub <- bind_rows(k_metrics_sub)
  k_metrics_df_sub$rank_silhouette <- rank(-k_metrics_df_sub$avg_silhouette, ties.method = "min")
  k_metrics_df_sub$rank_ari_mean <- rank(-k_metrics_df_sub$ari_mean, ties.method = "min")
  k_metrics_df_sub$rank_sum <- k_metrics_df_sub$rank_silhouette + k_metrics_df_sub$rank_ari_mean
  pick_sub <- k_metrics_df_sub %>% arrange(rank_sum, desc(ari_mean), desc(avg_silhouette), k) %>% slice(1)
  write.csv(k_metrics_df_sub, file.path(ADDON_TAB_DIR, paste0("additional_", label_stub, "_kmetrics.csv")), row.names = FALSE)
  write.csv(data.frame(excluded_gse = paste(excluded_gse, collapse = ";"), n_remaining = length(keep_ids), selected_k = pick_sub$k[1], stringsAsFactors = FALSE), file.path(ADDON_TAB_DIR, paste0("additional_", label_stub, "_summary.csv")), row.names = FALSE)
  sens_list[[label_stub]] <- data.frame(label = label_stub, excluded_gse = paste(excluded_gse, collapse = ";"), n_remaining = length(keep_ids), selected_k = pick_sub$k[1], stringsAsFactors = FALSE)
}
if (length(sens_list) > 0) write.csv(bind_rows(sens_list), file.path(ADDON_TAB_DIR, "additional_sensitivity_summary.csv"), row.names = FALSE)

if (file.exists(minor_deg_file) && file.exists(major_deg_file)) {
  minor_top <- read.csv(minor_deg_file, stringsAsFactors = FALSE, check.names = FALSE)
  minor_top <- minor_top[seq_len(min(20, nrow(minor_top))), c("GeneSymbol", "logFC", "adj.P.Val"), drop = FALSE]
  minor_top$Cluster <- minor_cluster
  major_top <- read.csv(major_deg_file, stringsAsFactors = FALSE, check.names = FALSE)
  major_top <- major_top[seq_len(min(20, nrow(major_top))), c("GeneSymbol", "logFC", "adj.P.Val"), drop = FALSE]
  major_top$Cluster <- major_cluster
  additional_marker_table <- bind_rows(minor_top, major_top) %>% select(Cluster, GeneSymbol, logFC, adj.P.Val)
  write.csv(additional_marker_table, file.path(ADDON_TAB_DIR, "additional_top20_markers_major_vs_minor_clusters.csv"), row.names = FALSE)
}

if (RUN_NETWORK) {
  sig_file <- file.path(ADDON_TAB_DIR, "additional_minor_cluster_myogenic_signature_genes.csv")
  if (file.exists(sig_file)) {
    sig_df <- read.csv(sig_file, stringsAsFactors = FALSE, check.names = FALSE)
    signature_genes <- intersect(unique(sig_df$GeneSymbol), rownames(exprs_combat))
    if (length(signature_genes) >= 5) {
      groups <- factor(cluster_df$Cluster[match(colnames(exprs_combat), cluster_df$SampleID)], levels = unique(cluster_df$Cluster))
      E_all <- exprs_combat[signature_genes, , drop = FALSE]
      rho_all <- suppressWarnings(cor(t(E_all), method = "spearman", use = "pairwise.complete.obs"))
      groups_resid <- as.factor(groups)
      X_resid <- model.matrix(~ groups_resid)
      E_res <- E_all
      for (i in seq_len(nrow(E_all))) {
        fit_resid <- stats::lm.fit(x = X_resid, y = as.numeric(E_all[i, ]))
        E_res[i, ] <- fit_resid$residuals
      }
      rho_res <- suppressWarnings(cor(t(E_res), method = "spearman", use = "pairwise.complete.obs"))
      gnames <- rownames(E_all)
      idx <- which(upper.tri(rho_all), arr.ind = TRUE)
      edge_df <- data.frame(from = gnames[idx[, 1]], to = gnames[idx[, 2]], overall_rho = rho_all[idx], residual_rho = rho_res[idx], stringsAsFactors = FALSE)
      edge_df$sign_agree_n <- ifelse(sign(edge_df$overall_rho) == sign(edge_df$residual_rho), 2L, 1L)
      edge_df$keep_robust <- with(edge_df, abs(overall_rho) >= NET_OVERALL_THR & abs(residual_rho) >= NET_RESIDUAL_THR & sign_agree_n >= 2)
      write.csv(edge_df, file.path(ADDON_TAB_DIR, "network_edge_metrics_cluster_adjusted.csv"), row.names = FALSE)
      robust_edges <- subset(edge_df, keep_robust)
      write.csv(robust_edges, file.path(ADDON_TAB_DIR, "network_edges_robust_cluster_adjusted.csv"), row.names = FALSE)
      if (nrow(robust_edges) > 0) {
        g_net <- igraph::graph_from_data_frame(robust_edges[, c("from", "to", "overall_rho", "residual_rho"), drop = FALSE], directed = FALSE)
        net_cent <- data.frame(gene = igraph::V(g_net)$name, degree = igraph::degree(g_net), betweenness = igraph::betweenness(g_net), closeness = igraph::closeness(g_net), stringsAsFactors = FALSE) %>% arrange(desc(degree), desc(betweenness))
        write.csv(net_cent, file.path(ADDON_TAB_DIR, "network_centrality_cluster_adjusted.csv"), row.names = FALSE)
        p_net <- ggraph::ggraph(g_net, layout = "fr") + ggraph::geom_edge_link(aes(width = abs(residual_rho), color = residual_rho), alpha = 0.4) + ggraph::geom_node_point(size = 4) + ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3) + ggraph::scale_edge_colour_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0) + theme_bw() + labs(title = "Robust coexpression network after cluster adjustment")
        ggsave(file.path(ADDON_PLOT_PNG_DIR, "network_cluster_adjusted_robust.png"), plot = p_net, width = 10, height = 9, dpi = PLOT_DPI)
        ggsave(file.path(ADDON_PLOT_PDF_DIR, "network_cluster_adjusted_robust.pdf"), plot = p_net, width = 10, height = 9, device = "pdf")
      }
    }
  }
}

writeLines(capture.output(sessionInfo()), con = file.path(ADDON_TAB_DIR, "discovery_additional_analyses_sessionInfo.txt"))
msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "03_discovery_additional_analyses finished")
cat(msg, "\n")
cat(msg, "\n", file = LOG_FILE, append = TRUE)
