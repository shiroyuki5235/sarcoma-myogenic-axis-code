# Sarcoma reported statistics: clinical models and paired signature comparisons
# Set the working directory to S1_Code. Run numbered modules in order.
# Analyses use the frozen S1 Data workbook; no expression preprocessing is refitted here.

#Module 00. Environment, packages and paths ----
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))
CODE_DIR <- normalizePath(Sys.getenv("SARCOMA_CODE_DIR", unset = getwd()), winslash = "/", mustWork = TRUE)
S1_DATA_FILE <- file.path(dirname(CODE_DIR), "S1_Data.xlsx")
if (!file.exists(S1_DATA_FILE)) stop("Place S1_Data.xlsx beside the extracted S1_Code folder.")
required_packages <- c("readxl", "survival", "pROC")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) install.packages(missing_packages, repos = "https://cloud.r-project.org")
stopifnot(all(vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)))
library(survival)
RESULTS_DIR <- file.path(CODE_DIR, "Results", paste0("Statistics_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (dir.exists(RESULTS_DIR)) stop("Results directory already exists; use a fresh run.")
OUTPUT_DIR <- file.path(RESULTS_DIR, "Clinical")
dir.create(OUTPUT_DIR, recursive = TRUE)
#Module 01. Read and verify the frozen inputs --------------------------------
tcga_data <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "TCGA_scores", skip = 4, n_max = 206))
existing_results <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "Cox_models", skip = 4, n_max = 48))
required_columns <- c("PATIENT_ID", "SAMPLE_ID", "histology", "age", "sex", "grade",
                      "myogenic_score", "score_z", "projected_axis", "cor_C1", "cor_C2",
                      "OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS")
if (!all(required_columns %in% names(tcga_data))) stop("Required TCGA columns are missing.")
if (nrow(tcga_data) != 206L || anyNA(tcga_data[required_columns])) stop("Expected 206 complete frozen patient records.")
if (anyDuplicated(tcga_data$PATIENT_ID) || anyDuplicated(tcga_data$SAMPLE_ID)) stop("Duplicate patient or sample IDs.")
numeric_columns <- c("age", "grade", "myogenic_score", "score_z", "cor_C1", "cor_C2", "OS_MONTHS", "DFS_MONTHS")
for (column_name in numeric_columns) {
  if (!is.numeric(tcga_data[[column_name]]) || any(!is.finite(tcga_data[[column_name]]))) stop("Invalid numeric values: ", column_name)
}
if (!all(tcga_data$sex %in% c("female", "male"))) stop("Unrecognized sex value.")
if (!all(tcga_data$grade %in% 1:3)) stop("Unrecognized grade value.")
if (!all(tcga_data$histology %in% c("STLMS", "ULMS", "DDLPS", "UPS", "MFS", "SS", "MPNST"))) stop("Unrecognized histology.")
if (!all(tcga_data$OS_STATUS %in% c("0:LIVING", "1:DECEASED"))) stop("Unrecognized OS status.")
if (!all(tcga_data$DFS_STATUS %in% c("0:DiseaseFree", "1:Recurred/Progressed"))) stop("Unrecognized DFS status.")
if (any(tcga_data$OS_MONTHS <= 0 | tcga_data$DFS_MONTHS <= 0)) stop("Nonpositive follow-up.")
if (max(abs(as.numeric(scale(tcga_data$myogenic_score)) - tcga_data$score_z)) > 1e-7) stop("Score standardization differs from the frozen 206-patient scale.")
expected_projection <- ifelse(tcga_data$cor_C2 >= tcga_data$cor_C1, "C2_like", "C1_like")
if (!identical(expected_projection, tcga_data$projected_axis)) stop("Projection labels do not match stored centroid correlations.")

#Module 02. Define the unchanged primary cohort -------------------------------
primary_data <- tcga_data[!tcga_data$histology %in% c("ULMS", "MPNST"), ]
primary_data$sex <- factor(primary_data$sex, levels = c("female", "male"))
primary_data$grade_factor <- factor(primary_data$grade, levels = c(1, 2, 3))
primary_data$projected_axis <- factor(primary_data$projected_axis, levels = c("C1_like", "C2_like"))
primary_data$histology <- factor(primary_data$histology)
primary_data$OS_event <- as.integer(primary_data$OS_STATUS == "1:DECEASED")
primary_data$DFS_event <- as.integer(primary_data$DFS_STATUS == "1:Recurred/Progressed")
if (nrow(primary_data) != 174L || sum(primary_data$OS_event) != 64L || sum(primary_data$DFS_event) != 103L) {
  stop("The primary cohort or endpoint counts differ from the manuscript.")
}
projection_composition <- as.data.frame(table(primary_data$histology, primary_data$projected_axis))
names(projection_composition) <- c("histology", "projected_axis", "n")
write.csv(projection_composition, file.path(OUTPUT_DIR, "Projection_composition.csv"), row.names = FALSE)
write.csv(primary_data[, c("PATIENT_ID", "SAMPLE_ID", "histology", "projected_axis", "grade", "OS_event", "DFS_event")],
          file.path(OUTPUT_DIR, "Analysis_population.csv"), row.names = FALSE)

#Module 03. Refit the original 27-exposure-test family ----
# Use the previously defined models and their complete-case cohorts; no model selection.
existing_exposures <- existing_results[existing_results$exposure_test, ]
stopifnot(nrow(existing_exposures) == 27L)
primary_coefficient_rows <- list(); primary_ph_rows <- list()
for (model_index in seq_len(nrow(existing_exposures))) {
  planned_model <- existing_exposures[model_index, ]
  model_data <- tcga_data
  if (planned_model$scenario == "primary_174_definition") model_data <- model_data[!model_data$histology %in% c("ULMS", "MPNST"), ]
  if (planned_model$scenario == "plus_ULMS") model_data <- model_data[model_data$histology != "MPNST", ]
  if (planned_model$scenario == "LMS_only") model_data <- model_data[model_data$histology %in% c("STLMS", "ULMS"), ]
  if (planned_model$scenario == "STLMS_only") model_data <- model_data[model_data$histology == "STLMS", ]
  model_data$sex <- factor(model_data$sex, levels = c("female", "male"))
  model_data$projected_axis <- factor(model_data$projected_axis, levels = c("C1_like", "C2_like"))
  model_data$histology <- factor(model_data$histology)
  model_data$time_months <- model_data[[paste0(planned_model$endpoint, "_MONTHS")]]
  model_data$event <- if (planned_model$endpoint == "OS") as.integer(model_data$OS_STATUS == "1:DECEASED") else as.integer(model_data$DFS_STATUS == "1:Recurred/Progressed")
  if (planned_model$model == "continuous_unadjusted") model_formula <- Surv(time_months, event) ~ score_z
  if (planned_model$model == "continuous_adjusted") model_formula <- Surv(time_months, event) ~ score_z + age + sex + grade + strata(histology)
  if (planned_model$model == "projected_axis_unadjusted") model_formula <- Surv(time_months, event) ~ projected_axis
  stopifnot(planned_model$model %in% c("continuous_unadjusted", "continuous_adjusted", "projected_axis_unadjusted"))
  primary_fit <- survival::coxph(model_formula, data = model_data, ties = "efron", na.action = na.fail, singular.ok = FALSE, x = TRUE, model = TRUE)
  primary_summary <- summary(primary_fit)
  stopifnot(primary_fit$n == planned_model$n, primary_fit$nevent == planned_model$events,
            abs(primary_summary$conf.int[planned_model$term, "exp(coef)"] - planned_model$HR) < 1e-6,
            abs(primary_summary$coefficients[planned_model$term, "Pr(>|z|)"] - planned_model$p_value) < 1e-6)
  primary_coefficient_rows[[model_index]] <- data.frame(scenario = planned_model$scenario, endpoint = planned_model$endpoint,
    model = planned_model$model, term = rownames(primary_summary$coefficients), n = primary_fit$n, events = primary_fit$nevent,
    HR = primary_summary$conf.int[, "exp(coef)"], CI_low = primary_summary$conf.int[, "lower .95"], CI_high = primary_summary$conf.int[, "upper .95"],
    p_value = primary_summary$coefficients[, "Pr(>|z|)"], exposure_test = rownames(primary_summary$coefficients) == planned_model$term, row.names = NULL)
  primary_ph <- survival::cox.zph(primary_fit, transform = "km", terms = TRUE, global = TRUE)
  primary_ph_rows[[model_index]] <- data.frame(scenario = planned_model$scenario, endpoint = planned_model$endpoint,
    model = planned_model$model, term = rownames(primary_ph$table), chisq = primary_ph$table[, "chisq"],
    df = primary_ph$table[, "df"], p = primary_ph$table[, "p"], row.names = NULL)
}
primary_coefficients <- do.call(rbind, primary_coefficient_rows)
primary_coefficients$all_exposure_family_BH <- NA_real_
primary_coefficients$all_exposure_family_BH[primary_coefficients$exposure_test] <- p.adjust(primary_coefficients$p_value[primary_coefficients$exposure_test], method = "BH", n = 27L)
stopifnot(max(abs(primary_coefficients$all_exposure_family_BH[primary_coefficients$exposure_test] - existing_exposures$all_exposure_family_BH)) < 1e-6)
write.csv(primary_coefficients, file.path(OUTPUT_DIR, "Primary_Cox_models.csv"), row.names = FALSE)
write.csv(do.call(rbind, primary_ph_rows), file.path(OUTPUT_DIR, "Primary_PH_diagnostics.csv"), row.names = FALSE)

#Module 04. Fit four planned sensitivity comparisons -------------------------
# Grade is categorical, allowing separate grade 2 and grade 3 coefficients.
# Each model has five coefficients: exposure, age, sex, two grade indicators.
# Histology gets separate baseline hazards. No new cutoffs or gene selection.
analysis_plan <- data.frame(
  endpoint = c("OS", "OS", "DFS", "DFS"),
  exposure = c("score_z", "projected_axis", "score_z", "projected_axis"),
  exposure_term = c("score_z", "projected_axisC2_like", "score_z", "projected_axisC2_like")
)
results <- data.frame(analysis_plan, n = 174L, events = c(64L, 64L, 103L, 103L),
                      parameters = 5L, events_per_parameter = c(64, 64, 103, 103) / 5,
                      HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_,
                      BH_q_four_planned_tests = NA_real_, status = "not_fitted", detail = "")
model_list <- vector("list", nrow(analysis_plan))
ph_list <- vector("list", nrow(analysis_plan))
coefficient_tables <- list()
ph_tables <- list()
for (model_index in seq_len(nrow(analysis_plan))) {
  model_data <- primary_data
  endpoint_name <- analysis_plan$endpoint[model_index]
  exposure_name <- analysis_plan$exposure[model_index]
  exposure_term <- analysis_plan$exposure_term[model_index]
  model_data$time_months <- model_data[[paste0(endpoint_name, "_MONTHS")]]
  model_data$event <- model_data[[paste0(endpoint_name, "_event")]]
  model_formula <- as.formula(paste("Surv(time_months, event) ~", exposure_name,
                                    "+ age + sex + grade_factor + strata(histology)"))
  fit_warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    survival::coxph(model_formula, data = model_data, ties = "efron", na.action = na.fail,
                   singular.ok = FALSE, x = TRUE, y = TRUE, model = TRUE),
    warning = function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(fit, "error") || length(fit_warnings) > 0L) {
    results$status[model_index] <- "unavailable_review_needed"
    results$detail[model_index] <- if (inherits(fit, "error")) conditionMessage(fit) else paste(fit_warnings, collapse = "; ")
    next
  }
  fit_summary <- summary(fit)
  if (length(coef(fit)) != 5L || fit$n != 174L || !exposure_term %in% names(coef(fit)) ||
      any(!is.finite(fit_summary$coefficients)) || any(!is.finite(fit_summary$conf.int))) {
    results$status[model_index] <- "unstable_or_unexpected_model"
    next
  }
  results$HR[model_index] <- fit_summary$conf.int[exposure_term, "exp(coef)"]
  results$CI_low[model_index] <- fit_summary$conf.int[exposure_term, "lower .95"]
  results$CI_high[model_index] <- fit_summary$conf.int[exposure_term, "upper .95"]
  results$p_value[model_index] <- fit_summary$coefficients[exposure_term, "Pr(>|z|)"]
  results$status[model_index] <- "fitted"
  model_list[[model_index]] <- fit
  coefficient_tables[[length(coefficient_tables) + 1L]] <- data.frame(
    endpoint = endpoint_name, exposure = exposure_name, term = rownames(fit_summary$coefficients),
    HR = fit_summary$conf.int[, "exp(coef)"], CI_low = fit_summary$conf.int[, "lower .95"],
    CI_high = fit_summary$conf.int[, "upper .95"], p_value = fit_summary$coefficients[, "Pr(>|z|)"], row.names = NULL)
  ph_result <- tryCatch(survival::cox.zph(fit, transform = "km", terms = TRUE, global = TRUE), error = function(e) e)
  if (inherits(ph_result, "error")) {
    results$detail[model_index] <- paste("PH diagnostic unavailable:", conditionMessage(ph_result))
  } else {
    ph_list[[model_index]] <- ph_result
    ph_tables[[length(ph_tables) + 1L]] <- data.frame(
      endpoint = endpoint_name, exposure = exposure_name, term = rownames(ph_result$table),
      chisq = ph_result$table[, "chisq"], df = ph_result$table[, "df"], p_value = ph_result$table[, "p"], row.names = NULL)
  }
}

#Module 05. Multiplicity and numerical outputs -------------------------------
# Retain the planned denominator of four even if one model is unavailable.
# These exploratory sensitivities do not replace the existing 27-test family.
available <- is.finite(results$p_value)
results$BH_q_four_planned_tests[available] <- p.adjust(results$p_value[available], method = "BH", n = 4L)
write.csv(results, file.path(OUTPUT_DIR, "Clinical_sensitivity_results.csv"), row.names = FALSE, na = "")
if (length(coefficient_tables) > 0L) write.csv(do.call(rbind, coefficient_tables), file.path(OUTPUT_DIR, "Model_coefficients.csv"), row.names = FALSE)
if (length(ph_tables) > 0L) write.csv(do.call(rbind, ph_tables), file.path(OUTPUT_DIR, "PH_diagnostics.csv"), row.names = FALSE)
saveRDS(list(plan = analysis_plan, results = results, models = model_list, PH = ph_list), file.path(OUTPUT_DIR, "Fitted_models.rds"))


stopifnot(all(results$status == "fitted"))

#Module 06. Submitted versus revised signature comparison settings ----
TABLE_DIR <- file.path(RESULTS_DIR, "Paired_signature"); dir.create(TABLE_DIR)
SEED <- 20261102L; set.seed(SEED)
BOOTSTRAP_REPLICATES <- 2000L
#Module 07. Exact paired TCGA matching and fixed signature inputs --------------
revised <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "TCGA_scores", skip = 4, n_max = 206))
original <- revised[, c("PATIENT_ID", "SAMPLE_ID", "submitted_myogenic_score")]
names(original)[3] <- "minor_cluster_score"
signature_membership <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "Signature_genes", skip = 4, n_max = 34))
original_genes <- signature_membership$GeneSymbol[signature_membership$in_original]
revised_genes <- signature_membership$GeneSymbol[signature_membership$in_revised]
stopifnot(nrow(original) == 206L, nrow(revised) == 206L,
  !anyNA(original$SAMPLE_ID), !anyNA(revised$SAMPLE_ID),
  !anyDuplicated(original$SAMPLE_ID), !anyDuplicated(revised$SAMPLE_ID),
  !anyDuplicated(original$PATIENT_ID), !anyDuplicated(revised$PATIENT_ID),
  setequal(original$SAMPLE_ID, revised$SAMPLE_ID),
  !anyNA(original_genes), !anyNA(revised_genes),
  !anyDuplicated(original_genes), !anyDuplicated(revised_genes),
  all(nzchar(original_genes)), all(nzchar(revised_genes)), length(revised_genes) == 30L)
original <- original[match(revised$SAMPLE_ID, original$SAMPLE_ID), , drop = FALSE]
stopifnot(identical(original$SAMPLE_ID, revised$SAMPLE_ID), identical(original$PATIENT_ID, revised$PATIENT_ID),
  all(is.finite(original$minor_cluster_score)), all(is.finite(revised$myogenic_score)),
  all(revised$histology %in% c("STLMS", "ULMS", "DDLPS", "UPS", "MFS", "SS", "MPNST")))
paired <- revised
paired$original_score <- original$minor_cluster_score
paired$revised_score <- revised$myogenic_score
# Both SD scales are fixed across the same full 206 patients before subsetting.
paired$original_score_z <- as.numeric(scale(paired$original_score))
paired$revised_score_z <- as.numeric(scale(paired$revised_score))
stopifnot(all(is.finite(paired$original_score_z)), all(is.finite(paired$revised_score_z)),
  max(abs(paired$revised_score_z - paired$score_z)) < 1e-7)
cohort_masks <- list(primary_174_definition = !paired$histology %in% c("ULMS", "MPNST"),
  plus_ULMS = paired$histology != "MPNST", all_histologies = rep(TRUE, nrow(paired)))
stopifnot(identical(as.integer(vapply(cohort_masks, sum, integer(1))), c(174L, 201L, 206L)))
write.csv(paired, file.path(TABLE_DIR, "02_paired_TCGA_inputs.csv"), row.names = FALSE)

#Module 08. Gene overlap and paired score agreement ---------------------------
gene_union <- union(original_genes, revised_genes)
gene_membership <- data.frame(GeneSymbol = gene_union,
  in_original = gene_union %in% original_genes, in_revised = gene_union %in% revised_genes)
write.csv(gene_membership, file.path(TABLE_DIR, "03_signature_gene_membership.csv"), row.names = FALSE)
write.csv(data.frame(original_n = length(original_genes), revised_n = length(revised_genes),
  overlap_n = length(intersect(original_genes, revised_genes)), union_n = length(gene_union),
  Jaccard = length(intersect(original_genes, revised_genes)) / length(gene_union)),
  file.path(TABLE_DIR, "03_signature_overlap_summary.csv"), row.names = FALSE)
agreement_rows <- list()
for (scenario in names(cohort_masks)) {
  d <- paired[cohort_masks[[scenario]], , drop = FALSE]
  agreement_rows[[scenario]] <- data.frame(scenario = scenario, n = nrow(d),
    Pearson = cor(d$original_score, d$revised_score, method = "pearson"),
    Spearman = cor(d$original_score, d$revised_score, method = "spearman"),
    mean_difference_z = mean(d$revised_score_z - d$original_score_z),
    median_absolute_difference_z = median(abs(d$revised_score_z - d$original_score_z)))
}
write.csv(do.call(rbind, agreement_rows), file.path(TABLE_DIR, "03_score_agreement.csv"), row.names = FALSE)

#Module 09. Paired AUC comparison without outcome-based signature selection ---
# Higher scores must indicate LMS for BOTH signatures: direction is fixed '<'.
# Bootstrap samples the SAME patients for both scores within each histology class.
# AUC and its difference concern histology, not survival prediction.
# pROC documentation: https://search.r-project.org/CRAN/refmans/pROC/html/roc.test.html
auc_rows <- list(); bootstrap_rows <- list()
for (scenario in names(cohort_masks)) {
  d <- paired[cohort_masks[[scenario]], , drop = FALSE]
  is_lms <- as.integer(d$histology %in% c("STLMS", "ULMS"))
  stopifnot(all(table(is_lms) > 1L))
  old_roc <- pROC::roc(is_lms, d$original_score, levels = c(0, 1), direction = "<", quiet = TRUE)
  new_roc <- pROC::roc(is_lms, d$revised_score, levels = c(0, 1), direction = "<", quiet = TRUE)
  comparison <- pROC::roc.test(old_roc, new_roc, method = "delong", paired = TRUE, alternative = "two.sided")
  draws <- matrix(NA_real_, nrow = BOOTSTRAP_REPLICATES, ncol = 3,
    dimnames = list(NULL, c("original_AUC", "revised_AUC", "difference_revised_minus_original")))
  control_indices <- which(is_lms == 0L); lms_indices <- which(is_lms == 1L)
  for (b in seq_len(BOOTSTRAP_REPLICATES)) {
    j <- c(sample(control_indices, length(control_indices), replace = TRUE), sample(lms_indices, length(lms_indices), replace = TRUE))
    draws[b, 1] <- as.numeric(pROC::auc(pROC::roc(is_lms[j], d$original_score[j], levels = c(0, 1), direction = "<", quiet = TRUE)))
    draws[b, 2] <- as.numeric(pROC::auc(pROC::roc(is_lms[j], d$revised_score[j], levels = c(0, 1), direction = "<", quiet = TRUE)))
    draws[b, 3] <- draws[b, 2] - draws[b, 1]
  }
  ci <- apply(draws, 2, quantile, probs = c(.025, .975), names = FALSE)
  auc_rows[[scenario]] <- data.frame(scenario = scenario, n = nrow(d), LMS_n = sum(is_lms),
    original_AUC = as.numeric(pROC::auc(old_roc)), original_CI_low = ci[1, 1], original_CI_high = ci[2, 1],
    revised_AUC = as.numeric(pROC::auc(new_roc)), revised_CI_low = ci[1, 2], revised_CI_high = ci[2, 2],
    difference_revised_minus_original = as.numeric(pROC::auc(new_roc)) - as.numeric(pROC::auc(old_roc)),
    difference_CI_low = ci[1, 3], difference_CI_high = ci[2, 3], paired_DeLong_P = comparison$p.value,
    bootstrap_replicates = BOOTSTRAP_REPLICATES)
  bootstrap_rows[[scenario]] <- data.frame(scenario = scenario, draw = seq_len(BOOTSTRAP_REPLICATES), draws)
}
auc_comparison <- do.call(rbind, auc_rows)
auc_comparison$BH_three_scenario_tests <- p.adjust(auc_comparison$paired_DeLong_P, method = "BH")
write.csv(auc_comparison, file.path(TABLE_DIR, "04_paired_AUC_comparison.csv"), row.names = FALSE)
write.csv(do.call(rbind, bootstrap_rows), file.path(TABLE_DIR, "04_paired_bootstrap_draws.csv"), row.names = FALSE)


#Module 10. Primary versus online comparison settings ----
OUTPUT_DIR <- file.path(RESULTS_DIR, "Online_comparison"); dir.create(OUTPUT_DIR)
#Module 11. Load and align IDs ----
tcga_data <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "TCGA_scores", skip = 4, n_max = 206))
geo_data <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "GEO_records", skip = 4, n_max = 1176))
signature_data <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "Signature_genes", skip = 4, n_max = 34))
reference_results <- as.data.frame(readxl::read_excel(S1_DATA_FILE, sheet = "AUC_comparison", skip = 18, n_max = 3))
stopifnot(nrow(tcga_data) == 206L, !anyDuplicated(tcga_data$PATIENT_ID), !anyDuplicated(tcga_data$SAMPLE_ID),
          !anyNA(tcga_data[, c("myogenic_score", "online_score", "projected_axis", "online_projection")]))
geo_data <- geo_data[geo_data$Sample_Decision == "KEEP", ]
stopifnot(nrow(geo_data) == 1115L, !anyDuplicated(geo_data$SampleID), !anyNA(geo_data$online_all_record_label))

#Module 12. GEO label agreement and signature overlap ----
geo_cross_tabulation <- as.data.frame(table(primary = geo_data$legacy_cluster, online = geo_data$online_all_record_label))
geo_changed <- geo_data[geo_data$legacy_cluster != geo_data$online_all_record_label, ]
signature_overlap <- data.frame(primary_n = sum(signature_data$in_revised), online_n = sum(signature_data$in_online),
                                shared_n = sum(signature_data$in_revised & signature_data$in_online))
stopifnot(nrow(geo_changed) == 2L, signature_overlap$primary_n == 30L, signature_overlap$online_n == 30L,
          signature_overlap$shared_n == 26L)
write.csv(geo_cross_tabulation, file.path(OUTPUT_DIR, "GEO_label_comparison.csv"), row.names = FALSE)
write.csv(geo_changed, file.path(OUTPUT_DIR, "Changed_GEO_labels.csv"), row.names = FALSE)
write.csv(signature_overlap, file.path(OUTPUT_DIR, "Signature_overlap.csv"), row.names = FALSE)

#Module 13. Paired score, AUC and centroid comparisons ----
scenario_masks <- list(primary_174 = !tcga_data$histology %in% c("ULMS", "MPNST"),
                       plus_ULMS = tcga_data$histology != "MPNST", all_histologies = rep(TRUE, nrow(tcga_data)))
scenario_results <- list()
for (scenario_index in seq_along(scenario_masks)) {
  d <- tcga_data[scenario_masks[[scenario_index]], ]
  is_lms <- d$histology %in% c("STLMS", "ULMS")
  primary_score_ranks <- rank(d$myogenic_score, ties.method = "average")
  online_score_ranks <- rank(d$online_score, ties.method = "average")
  # Empirical AUC = Mann-Whitney rank statistic; higher score indicates LMS.
  n_lms <- sum(is_lms); n_other <- sum(!is_lms)
  primary_auc <- (sum(primary_score_ranks[is_lms]) - n_lms * (n_lms + 1) / 2) / (n_lms * n_other)
  online_auc <- (sum(online_score_ranks[is_lms]) - n_lms * (n_lms + 1) / 2) / (n_lms * n_other)
  projection_changed <- (d$projected_axis == "C2_like") != (d$online_projection == "minor_like")
  scenario_results[[scenario_index]] <- data.frame(scenario = names(scenario_masks)[scenario_index], n = nrow(d), n_LMS = n_lms,
    primary_fixed_AUC = primary_auc, online_AUC = online_auc, Pearson_r = cor(d$myogenic_score, d$online_score),
    Spearman_rho = cor(d$myogenic_score, d$online_score, method = "spearman"), projection_changed = sum(projection_changed),
    online_minor_n = sum(d$online_projection == "minor_like"))
  write.csv(as.data.frame(table(primary = d$projected_axis, online = d$online_projection)),
            file.path(OUTPUT_DIR, paste0("Centroid_assignments_", names(scenario_masks)[scenario_index], ".csv")), row.names = FALSE)
}
comparison_results <- do.call(rbind, scenario_results)
comparison_columns <- setdiff(names(comparison_results), "scenario")
stopifnot(all(abs(as.matrix(comparison_results[comparison_columns]) - as.matrix(reference_results[comparison_columns])) < 1e-10))
write.csv(comparison_results, file.path(OUTPUT_DIR, "TCGA_scenario_comparison.csv"), row.names = FALSE)


#Module 14. Recompute reported online enrichment from fixed inputs ----
# These ranks come from the online C2-versus-C1 limma contrast; they are not primary GEO ranks.
GSEA_INPUT_DIR <- file.path(CODE_DIR, "Inputs", "Primary_project", "Input")
GSEA_OUTPUT_DIR <- file.path(RESULTS_DIR, "GSEA")
dir.create(GSEA_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
gsea_packages <- c("fgsea", "BiocParallel")
if (any(!vapply(gsea_packages, requireNamespace, logical(1), quietly = TRUE))) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(gsea_packages[!vapply(gsea_packages, requireNamespace, logical(1), quietly = TRUE)], ask = FALSE, update = FALSE)
}
online_de <- read.csv(file.path(GSEA_INPUT_DIR, "Online_GSEA", "DEG_C2_vsRest.csv"), check.names = FALSE)
online_ranks <- sort(setNames(online_de$t[is.finite(online_de$t)], online_de$GeneSymbol[is.finite(online_de$t)]), decreasing = TRUE)
stopifnot(length(online_ranks) > 1000L, !anyDuplicated(names(online_ranks)), !anyNA(online_ranks))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
for (collection in c("Hallmark", "Reactome")) {
  collection_pathways <- readRDS(file.path(GSEA_INPUT_DIR, "GeneSets", paste0(collection, ".rds")))
  collection_source <- readRDS(file.path(GSEA_INPUT_DIR, "GeneSets", paste0(collection, "_source_table.rds")))
  stopifnot(identical(unique(as.character(collection_source$db_version)), "2026.1.Hs"))
  set.seed(123L + match(collection, c("Hallmark", "Reactome")))
  collection_result <- fgsea::fgseaMultilevel(pathways = collection_pathways, stats = online_ranks,
    minSize = 10L, maxSize = 500L, BPPARAM = BiocParallel::SerialParam(progressbar = FALSE))
  collection_result <- as.data.frame(collection_result)
  stopifnot(nrow(collection_result) > 0L, all(is.finite(collection_result$NES)), all(is.finite(collection_result$padj)))
  collection_result$leadingEdge <- vapply(collection_result$leadingEdge, paste, character(1), collapse = ";")
  collection_result <- collection_result[order(collection_result$padj, -abs(collection_result$NES)), , drop = FALSE]
  write.csv(collection_result, file.path(GSEA_OUTPUT_DIR, paste0("GSEA_", collection, "_C2_vsRest.csv")), row.names = FALSE)
  print(data.frame(collection, release = unique(collection_source$db_version), pathways = nrow(collection_result),
    seed = 123L + match(collection, c("Hallmark", "Reactome")), fgsea = as.character(packageVersion("fgsea"))))
  if (collection == "Hallmark") print(collection_result[collection_result$pathway == "HALLMARK_MYOGENESIS", c("pathway", "NES", "ES", "pval", "padj", "size")])
}

#Module 15. Save the R session and display the principal results ----
writeLines(capture.output(sessionInfo()), file.path(RESULTS_DIR, "R_sessionInfo.txt"))
print(results); print(auc_comparison); print(comparison_results)
message("Reported statistics completed: ", RESULTS_DIR)
