run_step <- function(context) {
  source(file.path(context$revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)
  config <- load_revision_config(context)
  step <- "06_cutoff_and_dca"
  log_message(context, step, "Preparing IHS cutoff, label audit, and DCA analyses")

  train_expr <- read_expression_matrix(resolve_revision_path(context, config$paths$train_expr))
  train_pheno <- read_table_file(resolve_revision_path(context, config$paths$train_pheno))
  train_aligned <- align_expression_with_pheno(
    train_expr, train_pheno,
    config$cohorts$train$sample_col,
    config$cohorts$train$group_col,
    config$cohorts$train$positive_label
  )
  ext_expr <- read_expression_matrix(resolve_revision_path(context, config$paths$external_expr))
  ext_pheno <- read_table_file(resolve_revision_path(context, config$paths$external_pheno))
  ext_aligned <- align_expression_with_pheno(
    ext_expr, ext_pheno,
    config$cohorts$external$sample_col,
    config$cohorts$external$group_col,
    config$cohorts$external$positive_label
  )

  train_ihs_raw <- compute_ihs(train_aligned$expr, config$ihs$genes)
  ext_ihs_raw <- compute_ihs(ext_aligned$expr, config$ihs$genes)
  train_slc2a1 <- extract_feature_matrix(train_aligned$expr, c("SLC2A1"))[, 1]
  ext_slc2a1 <- extract_feature_matrix(ext_aligned$expr, c("SLC2A1"))[, 1]

  folds <- make_folds(train_aligned$y, 5, config$revision$seed + 200)
  oof_ihs_raw <- rep(NA_real_, length(train_aligned$y))
  for (fold_id in seq_along(folds)) {
    test_idx <- folds[[fold_id]]
    oof_ihs_raw[test_idx] <- train_ihs_raw[test_idx]
  }
  train_oof_auc_raw <- roc_auc_ci(train_aligned$y, oof_ihs_raw)["auc"]

  flip_sign <- is.finite(train_oof_auc_raw) && (train_oof_auc_raw < 0.5)
  locked_multiplier <- if (flip_sign) -1 else 1
  direction_reason <- if (flip_sign) {
    "Training OOF raw-score AUC < 0.5; locked IHS direction flipped (IHS := -IHS) using training cohort only."
  } else {
    "Training OOF raw-score AUC >= 0.5; kept original IHS direction."
  }

  direction_df <- data.frame(
    decision_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    training_cohort_id = config$cohorts$train$id,
    training_group_col = config$cohorts$train$group_col,
    training_positive_label = config$cohorts$train$positive_label,
    training_control_label = if (!is.null(config$cohorts$train$control_label)) config$cohorts$train$control_label else "non-positive_label",
    training_oof_auc_rawscore = train_oof_auc_raw,
    flip_sign = flip_sign,
    locked_multiplier = locked_multiplier,
    decision_rule = "Flip sign only if training OOF raw-score AUC < 0.5; external cohorts never used for direction decision.",
    rationale = direction_reason,
    external_used_for_decision = FALSE,
    stringsAsFactors = FALSE
  )
  direction_path <- file.path(context$tables_dir, "IHS_direction_decision.csv")
  write.csv(direction_df, direction_path, row.names = FALSE)

  train_ihs_locked <- train_ihs_raw * locked_multiplier
  ext_ihs_locked <- ext_ihs_raw * locked_multiplier

  build_audit_row <- function(cohort_cfg, aligned_obj, score_raw, score_locked, group_values) {
    raw_auc <- roc_auc_ci(aligned_obj$y, score_raw)["auc"]
    locked_auc <- roc_auc_ci(aligned_obj$y, score_locked)["auc"]
    group_tab <- table(group_values)
    group_detail <- paste(sprintf("%s=%d", names(group_tab), as.integer(group_tab)), collapse = "; ")
    data.frame(
      cohort_id = cohort_cfg$id,
      cohort_label = cohort_cfg$label,
      sample_col_source = cohort_cfg$sample_col,
      group_col_source = cohort_cfg$group_col,
      case_definition = cohort_cfg$positive_label,
      control_definition = if (!is.null(cohort_cfg$control_label)) cohort_cfg$control_label else paste0("not_", cohort_cfg$positive_label),
      label_regex_or_rule = "exact label match in configured group column",
      n_case = sum(aligned_obj$y == 1),
      n_control = sum(aligned_obj$y == 0),
      group_value_counts = group_detail,
      ihs_auc_raw = raw_auc,
      ihs_auc_locked = locked_auc,
      stringsAsFactors = FALSE
    )
  }

  audit_df <- rbind(
    build_audit_row(config$cohorts$train, train_aligned, train_ihs_raw, train_ihs_locked, train_aligned$pheno[[config$cohorts$train$group_col]]),
    build_audit_row(config$cohorts$external, ext_aligned, ext_ihs_raw, ext_ihs_locked, ext_aligned$pheno[[config$cohorts$external$group_col]])
  )
  audit_path <- file.path(context$tables_dir, "IHS_label_mapping_audit.csv")
  write.csv(audit_df, audit_path, row.names = FALSE)

  for (i in seq_len(nrow(audit_df))) {
    log_message(
      context,
      step,
      sprintf(
        "Label audit | %s | group_col=%s | case=%s | control=%s | n_case=%d n_control=%d | AUC_raw=%.3f AUC_locked=%.3f",
        audit_df$cohort_id[i],
        audit_df$group_col_source[i],
        audit_df$case_definition[i],
        audit_df$control_definition[i],
        audit_df$n_case[i],
        audit_df$n_control[i],
        audit_df$ihs_auc_raw[i],
        audit_df$ihs_auc_locked[i]
      )
    )
  }

  oof_prob_ihs <- rep(NA_real_, length(train_aligned$y))
  oof_prob_slc2a1 <- rep(NA_real_, length(train_aligned$y))
  for (fold_id in seq_along(folds)) {
    test_idx <- folds[[fold_id]]
    train_idx <- setdiff(seq_along(train_aligned$y), test_idx)
    ihs_fit <- glm(y ~ ihs, family = binomial(), data = data.frame(y = train_aligned$y[train_idx], ihs = train_ihs_locked[train_idx]))
    slc_fit <- glm(y ~ slc2a1, family = binomial(), data = data.frame(y = train_aligned$y[train_idx], slc2a1 = train_slc2a1[train_idx]))
    oof_prob_ihs[test_idx] <- predict(ihs_fit, newdata = data.frame(ihs = train_ihs_locked[test_idx]), type = "response")
    oof_prob_slc2a1[test_idx] <- predict(slc_fit, newdata = data.frame(slc2a1 = train_slc2a1[test_idx]), type = "response")
  }

  final_ihs_fit <- glm(y ~ ihs, family = binomial(), data = data.frame(y = train_aligned$y, ihs = train_ihs_locked))
  final_slc_fit <- glm(y ~ slc2a1, family = binomial(), data = data.frame(y = train_aligned$y, slc2a1 = train_slc2a1))
  ext_prob_ihs <- predict(final_ihs_fit, newdata = data.frame(ihs = ext_ihs_locked), type = "response")
  ext_prob_slc2a1 <- predict(final_slc_fit, newdata = data.frame(slc2a1 = ext_slc2a1), type = "response")

  youden_prob <- youden_threshold(train_aligned$y, oof_prob_ihs)
  sens90_prob <- find_threshold_for_sensitivity(train_aligned$y, oof_prob_ihs, target = 0.90)
  coef_ihs <- coef(final_ihs_fit)
  youden_ihs_cutoff <- (qlogis(clip_probs(youden_prob)) - coef_ihs[1]) / coef_ihs[2]
  sens90_ihs_cutoff <- (qlogis(clip_probs(sens90_prob)) - coef_ihs[1]) / coef_ihs[2]

  cutoff_df <- rbind(
    cbind(data.frame(dataset = "OOF_train", cutoff_type = "Youden"), threshold_with_ci(train_aligned$y, oof_prob_ihs, youden_prob, "IHS_probability")),
    cbind(data.frame(dataset = "External", cutoff_type = "Youden"), threshold_with_ci(ext_aligned$y, ext_prob_ihs, youden_prob, "IHS_probability")),
    cbind(data.frame(dataset = "OOF_train", cutoff_type = "Sensitivity_0.90"), threshold_with_ci(train_aligned$y, oof_prob_ihs, sens90_prob, "IHS_probability")),
    cbind(data.frame(dataset = "External", cutoff_type = "Sensitivity_0.90"), threshold_with_ci(ext_aligned$y, ext_prob_ihs, sens90_prob, "IHS_probability"))
  )
  cutoff_df$locked_probability_threshold <- ifelse(cutoff_df$cutoff_type == "Youden", youden_prob, sens90_prob)
  cutoff_df$locked_ihs_score_cutoff <- ifelse(cutoff_df$cutoff_type == "Youden", youden_ihs_cutoff, sens90_ihs_cutoff)
  cutoff_df$ihs_direction_multiplier <- locked_multiplier
  cutoff_path <- file.path(context$tables_dir, "Table_Cutoff_Performance_Metrics.csv")
  write.csv(cutoff_df, cutoff_path, row.names = FALSE)

  dca_thresholds <- seq(0.05, 0.90, by = 0.05)
  ihs_model_name <- "Locked 10-gene IHS score"
  dca_df <- rbind(
    decision_curve_df(ext_aligned$y, ext_prob_ihs, model_name = ihs_model_name, thresholds = dca_thresholds),
    decision_curve_df(ext_aligned$y, ext_prob_slc2a1, model_name = "SLC2A1-only", thresholds = dca_thresholds),
    baseline_decision_curves(ext_aligned$y, thresholds = dca_thresholds)
  )
  dca_path <- file.path(context$tables_dir, "Table_DCA_NetBenefit_Values.csv")
  write.csv(dca_df, dca_path, row.names = FALSE)

  # Incremental DCA summary at matched thresholds to quantify descriptive gains.
  get_model_nb <- function(model_name) {
    sub <- dca_df[dca_df$model == model_name, c("threshold", "net_benefit"), drop = FALSE]
    colnames(sub)[2] <- paste0("net_benefit_", gsub("[^A-Za-z0-9]+", "_", tolower(model_name)))
    sub
  }
  dca_cmp <- Reduce(
    function(x, y) merge(x, y, by = "threshold", all = TRUE),
    list(
      get_model_nb(ihs_model_name),
      get_model_nb("SLC2A1-only"),
      get_model_nb("Treat all"),
      get_model_nb("Treat none")
    )
  )
  ihs_nb_col <- paste0("net_benefit_", gsub("[^A-Za-z0-9]+", "_", tolower(ihs_model_name)))
  dca_cmp$delta_ihs_vs_slc2a1 <- dca_cmp[[ihs_nb_col]] - dca_cmp$net_benefit_slc2a1_only
  dca_cmp$delta_ihs_vs_treat_all <- dca_cmp[[ihs_nb_col]] - dca_cmp$net_benefit_treat_all
  dca_cmp$delta_ihs_vs_treat_none <- dca_cmp[[ihs_nb_col]] - dca_cmp$net_benefit_treat_none
  dca_cmp <- dca_cmp[order(dca_cmp$threshold), , drop = FALSE]
  dca_cmp_path <- file.path(context$tables_dir, "Table_DCA_Incremental_vs_SLC2A1.csv")
  write.csv(dca_cmp, dca_cmp_path, row.names = FALSE)

  key_thresholds <- c(0.10, 0.20, 0.30, 0.50, 0.70)
  key_mask <- vapply(dca_cmp$threshold, function(th) any(abs(th - key_thresholds) < 1e-10), logical(1))
  dca_key <- dca_cmp[key_mask, c(
    "threshold",
    ihs_nb_col,
    "net_benefit_slc2a1_only",
    "delta_ihs_vs_slc2a1",
    "net_benefit_treat_all",
    "net_benefit_treat_none"
  ), drop = FALSE]
  dca_key_path <- file.path(context$tables_dir, "Table_DCA_Incremental_vs_SLC2A1_KeyThresholds.csv")
  write.csv(dca_key, dca_key_path, row.names = FALSE)

  dca_key_latex <- dca_key
  key_num <- vapply(dca_key_latex, is.numeric, logical(1))
  dca_key_latex[key_num] <- lapply(dca_key_latex[key_num], function(x) sprintf("%.3f", x))
  colnames(dca_key_latex) <- c("threshold", "NBIHS", "NBSLC2A1", "DeltaIHSvsSLC2A1", "NBTreatAll", "NBTreatNone")
  dca_key_latex_path <- file.path(context$tables_dir, "Table_DCA_Incremental_vs_SLC2A1_KeyThresholds.tex")
  simple_latex_table(
    dca_key_latex,
    dca_key_latex_path,
    caption = "External-cohort DCA at representative thresholds using an auxiliary reviewer-boundary score versus SLC2A1-only and treat-all/none baselines (not the primary locked 10-gene model claim).",
    label = "tab:rev1_dca_keythresholds"
  )

  better_vs_slc2a1 <- sum(dca_cmp$delta_ihs_vs_slc2a1 > 0, na.rm = TRUE)
  better_vs_none <- sum(dca_cmp$delta_ihs_vs_treat_none > 0, na.rm = TRUE)
  dca_summary_path <- file.path(context$tables_dir, "DCA_Incremental_Summary.txt")
  write_text_file(dca_summary_path, c(
    sprintf("Threshold points evaluated: %d", nrow(dca_cmp)),
    sprintf("Threshold points where locked 10-gene IHS net benefit > SLC2A1-only: %d", better_vs_slc2a1),
    sprintf("Threshold points where locked 10-gene IHS net benefit > Treat-none: %d", better_vs_none),
    "Interpretation: in this case-enriched external cohort, DCA is descriptive and not intended for population-level net-benefit inference."
  ))

  p <- ggplot2::ggplot(dca_df, ggplot2::aes(x = threshold, y = net_benefit, color = model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = "Decision-curve analysis on external MASLD cohort",
      subtitle = "FIB-4/NFS unavailable in public transcriptomic metadata; comparator reduced to SLC2A1-only baseline",
      x = "Risk threshold",
      y = "Net benefit"
    ) +
    ggplot2::theme_bw(base_size = 11)
  dca_fig_path <- file.path(context$figures_supp_dir, "DCA_curve_external_realprob.pdf")
  ggplot2::ggsave(dca_fig_path, plot = p, width = 7.5, height = 5.2)

  ihs_pred_df <- data.frame(
    sample_id = c(train_aligned$sample_ids, ext_aligned$sample_ids),
    cohort = c(rep(config$cohorts$train$id, length(train_aligned$sample_ids)), rep(config$cohorts$external$id, length(ext_aligned$sample_ids))),
    split = c(rep("train", length(train_aligned$sample_ids)), rep("external", length(ext_aligned$sample_ids))),
    observed = c(train_aligned$y, ext_aligned$y),
    ihs_score_raw = c(train_ihs_raw, ext_ihs_raw),
    ihs_score_locked = c(train_ihs_locked, ext_ihs_locked),
    ihs_probability = c(oof_prob_ihs, ext_prob_ihs),
    slc2a1_probability = c(oof_prob_slc2a1, ext_prob_slc2a1),
    stringsAsFactors = FALSE
  )
  write.csv(ihs_pred_df, file.path(context$tables_dir, "IHS_predictions_train_external.csv"), row.names = FALSE)

  guidance_lines <- c(
    sprintf("IHS direction lock multiplier: %d", locked_multiplier),
    sprintf("Direction decision basis (training OOF raw-score AUC): %.3f", train_oof_auc_raw),
    sprintf("Primary locked IHS cutoff (Youden-derived): IHS score >= %.3f", youden_ihs_cutoff),
    sprintf("Secondary high-sensitivity cutoff: IHS score >= %.3f", sens90_ihs_cutoff),
    "Clinical use note:",
    "1. Compute the raw auxiliary reviewer-boundary score using the configured marker set in config.yml.",
    "2. Apply the locked direction multiplier from training OOF to obtain direction-locked IHS.",
    "3. Apply the locked threshold from training OOF analysis to external cohorts without re-optimizing.",
    "4. Because FIB-4 and NFS covariates were unavailable in these public transcriptomic cohorts, current DCA compares this auxiliary score only against a single-gene baseline and treat-all/treat-none strategies.",
    "5. Primary manuscript claims remain tied to the locked 10-gene model output; this auxiliary score is retained only for reviewer-boundary stress testing."
  )
  guidance_path <- file.path(context$tables_dir, "Clinical_Guidance_IHS_Use.txt")
  write_text_file(guidance_path, guidance_lines)

  log_message(context, step, sprintf("Direction decision written: %s", direction_path))
  log_message(context, step, sprintf("Label mapping audit written: %s", audit_path))
  log_message(context, step, sprintf("Cutoff metrics written: %s", cutoff_path))
  log_message(context, step, sprintf("DCA incremental table written: %s", dca_key_path))
  log_message(context, step, sprintf("DCA figure written: %s", dca_fig_path))
}
