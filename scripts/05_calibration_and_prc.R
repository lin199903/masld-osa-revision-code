run_step <- function(context) {
  source(file.path(context$revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)
  config <- load_revision_config(context)
  step <- "05_calibration_and_prc"
  log_message(context, step, "Loading OOF and holdout predictions")

  # Seed locked manuscript Figure 4 assets without overwriting existing files.
  locked_seed <- ensure_locked_main_fig4_assets(context, config, overwrite = FALSE, strict = FALSE)
  if (!isTRUE(locked_seed$skipped) && length(locked_seed$copied) > 0L) {
    log_message(context, step, sprintf("Seeded locked Figure 4 assets: %d file(s)", length(locked_seed$copied)))
  }
  if (!isTRUE(locked_seed$skipped) && length(locked_seed$missing) > 0L) {
    log_message(
      context,
      step,
      sprintf("Locked Figure 4 source assets not found yet (non-fatal at Step 05): %s", paste(locked_seed$missing, collapse = "; ")),
      level = "WARN"
    )
  }

  oof_df <- read.csv(file.path(context$tables_dir, "OOF_predictions_with_probabilities.csv"), check.names = FALSE)
  holdout_df <- read.csv(file.path(context$tables_dir, "External_test_predictions.csv"), check.names = FALSE)
  external_df <- holdout_df[holdout_df$cohort == "GSE135251", , drop = FALSE]

  metrics_df <- do.call(rbind, lapply(split(rbind(oof_df, external_df), c(rep("OOF", nrow(oof_df)), rep("External", nrow(external_df)))), function(df) {
    auc_ci <- roc_auc_ci(df$observed, df$predicted_prob)
    cal <- calibration_metrics(df$observed, df$predicted_prob)
    data.frame(
      dataset = unique(ifelse(df$cohort == "GSE135251", "External", "OOF")),
      auc = auc_ci["auc"],
      auc_ci_low = auc_ci["ci_low"],
      auc_ci_high = auc_ci["ci_high"],
      pr_auc_case = pr_auc_case(df$observed, df$predicted_prob),
      pr_auc_control = pr_auc_case(1L - df$observed, 1 - df$predicted_prob),
      calibration_intercept = cal["intercept"],
      calibration_slope = cal["slope"],
      brier = cal["brier"],
      stringsAsFactors = FALSE
    )
  }))
  metrics_path <- file.path(context$tables_dir, "Table_Calibration_Metrics.csv")
  write.csv(metrics_df, metrics_path, row.names = FALSE)

  oof_roc <- pROC::roc(oof_df$observed, oof_df$predicted_prob, quiet = TRUE, direction = "<")
  ext_roc <- pROC::roc(external_df$observed, external_df$predicted_prob, quiet = TRUE, direction = "<")
  oof_cal <- calibration_curve_df(oof_df$observed, oof_df$predicted_prob)
  ext_cal <- calibration_curve_df(external_df$observed, external_df$predicted_prob)

  combined_path <- file.path(context$figures_main_dir, "Combined_ROC_Calibration_SVM_RFE_LDA_AUX_canonical.pdf")
  combined_main_path <- file.path(context$figures_main_dir, "Combined_ROC_Calibration_SVM_RFE_LDA.pdf")
  grDevices::pdf(combined_path, width = 10, height = 4.8)
  old_par <- graphics::par(no.readonly = TRUE)
  graphics::par(
    mfrow = c(1, 2),
    mar = c(4, 4, 3, 1),
    cex.main = 1.15,
    cex.lab = 1.1,
    cex.axis = 1.05
  )
  graphics::plot(oof_roc, col = "#C44E52", lwd = 2, main = "ROC", legacy.axes = TRUE)
  graphics::plot(ext_roc, col = "#2F6C8F", lwd = 2, add = TRUE, legacy.axes = TRUE)
  graphics::abline(0, 1, lty = 2, col = "grey70")
  graphics::legend(
    "bottomright",
    legend = c(
      sprintf("OOF AUC = %.3f", as.numeric(pROC::auc(oof_roc))),
      sprintf("External AUC = %.3f", as.numeric(pROC::auc(ext_roc)))
    ),
    col = c("#C44E52", "#2F6C8F"),
    lwd = 2,
    cex = 1.24,
    bty = "n"
  )
  graphics::plot(c(0, 1), c(0, 1), type = "n", xlab = "Mean predicted risk", ylab = "Observed event rate", main = "Calibration")
  graphics::abline(0, 1, lty = 2, col = "grey70")
  graphics::lines(oof_cal$prob, oof_cal$obs, type = "b", pch = 16, col = "#C44E52", lwd = 2)
  graphics::lines(ext_cal$prob, ext_cal$obs, type = "b", pch = 17, col = "#2F6C8F", lwd = 2)
  graphics::legend("topleft", legend = c("OOF", "External"), col = c("#C44E52", "#2F6C8F"), pch = c(16, 17), lwd = 2, cex = 1.16, bty = "n")
  graphics::par(old_par)
  grDevices::dev.off()
  file.copy(combined_path, combined_main_path, overwrite = TRUE)

  prc_case_path <- file.path(context$figures_main_dir, "PRC_Combined_SVM_RFE_LDA_AUX_canonical.pdf")
  prc_case_main_path <- file.path(context$figures_main_dir, "PRC_Combined_SVM_RFE_LDA.pdf")
  grDevices::pdf(prc_case_path, width = 6.5, height = 5.2)
  old_par <- graphics::par(no.readonly = TRUE)
  graphics::par(cex.main = 1.15, cex.lab = 1.1, cex.axis = 1.05)
  oof_pr <- PRROC::pr.curve(scores.class0 = oof_df$predicted_prob[oof_df$observed == 1], scores.class1 = oof_df$predicted_prob[oof_df$observed == 0], curve = TRUE)
  ext_pr <- PRROC::pr.curve(scores.class0 = external_df$predicted_prob[external_df$observed == 1], scores.class1 = external_df$predicted_prob[external_df$observed == 0], curve = TRUE)
  graphics::plot(oof_pr, color = "#C44E52", lwd = 2, main = "Case-positive precision-recall")
  graphics::lines(ext_pr$curve[, 1], ext_pr$curve[, 2], col = "#2F6C8F", lwd = 2)
  graphics::legend(
    "bottomleft",
    legend = c(
      sprintf("OOF PRAUC = %.3f", oof_pr$auc.integral),
      sprintf("External PRAUC = %.3f", ext_pr$auc.integral)
    ),
    col = c("#C44E52", "#2F6C8F"),
    lwd = 2,
    cex = 1.24,
    bty = "n"
  )
  graphics::par(old_par)
  grDevices::dev.off()
  file.copy(prc_case_path, prc_case_main_path, overwrite = TRUE)

  prc_control_path <- file.path(context$figures_main_dir, "PRC_Inverted_ControlPositive_SVM_RFE_LDA_AUX_canonical.pdf")
  prc_control_main_path <- file.path(context$figures_main_dir, "PRC_Inverted_ControlPositive_SVM_RFE_LDA.pdf")
  grDevices::pdf(prc_control_path, width = 6.5, height = 5.2)
  old_par <- graphics::par(no.readonly = TRUE)
  graphics::par(cex.main = 1.15, cex.lab = 1.1, cex.axis = 1.05)
  oof_pr_ctrl <- PRROC::pr.curve(scores.class0 = 1 - oof_df$predicted_prob[oof_df$observed == 0], scores.class1 = 1 - oof_df$predicted_prob[oof_df$observed == 1], curve = TRUE)
  ext_pr_ctrl <- PRROC::pr.curve(scores.class0 = 1 - external_df$predicted_prob[external_df$observed == 0], scores.class1 = 1 - external_df$predicted_prob[external_df$observed == 1], curve = TRUE)
  graphics::plot(oof_pr_ctrl, color = "#4F6D7A", lwd = 2, main = "Control-positive precision-recall")
  graphics::lines(ext_pr_ctrl$curve[, 1], ext_pr_ctrl$curve[, 2], col = "#A0A0A0", lwd = 2)
  graphics::legend(
    "bottomleft",
    inset = 0.02,
    legend = c(
      sprintf("OOF control-PRAUC = %.3f", oof_pr_ctrl$auc.integral),
      sprintf("External control-PRAUC = %.3f", ext_pr_ctrl$auc.integral)
    ),
    col = c("#4F6D7A", "#A0A0A0"),
    lwd = 2,
    cex = 1.24,
    bty = "o",
    bg = grDevices::adjustcolor("white", alpha.f = 0.9),
    box.col = "grey70"
  )
  graphics::par(old_par)
  grDevices::dev.off()
  file.copy(prc_control_path, prc_control_main_path, overwrite = TRUE)

  log_message(context, step, sprintf("Calibration table written: %s", metrics_path))
  log_message(context, step, sprintf("Auxiliary ROC/calibration figure written (non-manuscript Figure 4): %s", combined_path))
  log_message(context, step, sprintf("Auxiliary case-positive PRC figure written (non-manuscript Figure 4): %s", prc_case_path))
  log_message(context, step, sprintf("Auxiliary control-positive PRC figure written (non-manuscript Figure 4): %s", prc_control_path))
  log_message(context, step, sprintf("Figure 4B updated from new plotting code: %s", combined_main_path))
  log_message(context, step, sprintf("Figure 4C updated from new plotting code: %s", prc_case_main_path))
  log_message(context, step, sprintf("Figure 4D updated from new plotting code: %s", prc_control_main_path))
}
