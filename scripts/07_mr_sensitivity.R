run_step <- function(context) {
  source(file.path(context$revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)
  config <- load_revision_config(context)
  step <- "07_mr_sensitivity"
  mr_root <- resolve_revision_path(context, config$paths$mr_root)
  log_message(context, step, sprintf("Using MR result root: %s", mr_root))

  mr_complete <- read.csv(file.path(mr_root, "Publication_Tables", "MR_Complete_Results.csv"), check.names = FALSE)
  sensitivity_summary <- read.csv(file.path(mr_root, "Sensitivity", "MR_Sensitivity_Summary.csv"), check.names = FALSE)
  egger_osa <- read.csv(file.path(mr_root, "Egger_Intercept_OSA_to_NAFLD.csv"), check.names = FALSE)
  egger_nafld <- read.csv(file.path(mr_root, "Egger_Intercept_NAFLD_to_OSA.csv"), check.names = FALSE)
  het_osa <- read.csv(file.path(mr_root, "Heterogeneity_OSA_to_NAFLD.csv"), check.names = FALSE)
  het_nafld <- read.csv(file.path(mr_root, "Heterogeneity_NAFLD_to_OSA.csv"), check.names = FALSE)
  strength <- read.csv(file.path(mr_root, "Publication_Tables", "Instrument_Strength_Summary.csv"), check.names = FALSE)
  loo_osa <- read.csv(file.path(mr_root, "LeaveOneOut_OSA_to_NAFLD.csv"), check.names = FALSE)
  loo_nafld <- read.csv(file.path(mr_root, "LeaveOneOut_NAFLD_to_OSA.csv"), check.names = FALSE)

  normalize_direction <- function(x) {
    x <- as.character(x)
    if (grepl("OSA", x, ignore.case = TRUE) && grepl("MASLD|NAFLD", x, ignore.case = TRUE) && !grepl("^MASLD|^NAFLD", x, ignore.case = TRUE)) {
      return("OSA_to_MASLD")
    }
    if (grepl("MASLD|NAFLD", x, ignore.case = TRUE) && grepl("OSA", x, ignore.case = TRUE) && (grepl("^MASLD|^NAFLD", x, ignore.case = TRUE) || grepl("to_OSA", x, ignore.case = TRUE))) {
      return("MASLD_to_OSA")
    }
    x
  }

  mr_complete$direction_std <- vapply(mr_complete$Direction, normalize_direction, character(1))
  sensitivity_summary$direction_std <- ifelse(grepl("^OSA", sensitivity_summary$Direction), "OSA_to_MASLD", "MASLD_to_OSA")
  strength$direction_std <- ifelse(grepl("^OSA", strength$Direction), "OSA_to_MASLD", "MASLD_to_OSA")

  add_direction_summary <- function(direction_std, egger, heterogeneity, loo, presso_row, steiger_row, strength_row) {
    subset_complete <- mr_complete[mr_complete$direction_std == direction_std, , drop = FALSE]
    subset_complete$egger_intercept <- egger$Intercept[1]
    subset_complete$egger_intercept_p <- egger$P[1]
    subset_complete$heterogeneity_q <- heterogeneity$Q[1]
    subset_complete$heterogeneity_p <- heterogeneity$Q_pval[1]
    subset_complete$I2 <- heterogeneity$I2[1]
    subset_complete$steiger_prop_correct <- steiger_row$Steiger_Prop_Correct[1]
    subset_complete$steiger_p <- steiger_row$Steiger_P[1]
    subset_complete$steiger_interpretation <- steiger_row$Steiger_Interpretation[1]
    subset_complete$presso_global_p <- presso_row$PRESSO_Global_P[1]
    subset_complete$presso_n_outliers <- presso_row$PRESSO_N_Outliers[1]
    subset_complete$presso_distortion_p <- presso_row$PRESSO_Distortion_P[1]
    subset_complete$f_min <- strength_row$F_min[1]
    subset_complete$f_median <- strength_row$F_median[1]
    subset_complete$weak_ivs <- strength_row$Weak_IVs[1]
    subset_complete$loo_or_min <- min(loo$OR, na.rm = TRUE)
    subset_complete$loo_or_max <- max(loo$OR, na.rm = TRUE)
    subset_complete
  }

  osa_row <- sensitivity_summary[sensitivity_summary$direction_std == "OSA_to_MASLD", , drop = FALSE]
  nafld_row <- sensitivity_summary[sensitivity_summary$direction_std == "MASLD_to_OSA", , drop = FALSE]
  strength_osa <- strength[strength$direction_std == "OSA_to_MASLD", , drop = FALSE]
  strength_nafld <- strength[strength$direction_std == "MASLD_to_OSA", , drop = FALSE]

  combined_df <- rbind(
    add_direction_summary("OSA_to_MASLD", egger_osa, het_osa, loo_osa, osa_row, osa_row, strength_osa),
    add_direction_summary("MASLD_to_OSA", egger_nafld, het_nafld, loo_nafld, nafld_row, nafld_row, strength_nafld)
  )
  combined_df$direction_label <- ifelse(combined_df$direction_std == "OSA_to_MASLD", "OSA to MASLD", "MASLD to OSA")
  combined_df$Direction <- combined_df$direction_label

  out_path <- file.path(context$tables_dir, "Table_MR_Sensitivity_Summary.csv")
  csv_df <- combined_df[, c(
    "Direction", "Method", "N_SNPs", "Beta", "SE", "OR", "OR_lci", "OR_uci", "P",
    "egger_intercept", "egger_intercept_p", "heterogeneity_q", "heterogeneity_p", "I2",
    "steiger_prop_correct", "steiger_p", "steiger_interpretation",
    "presso_global_p", "presso_n_outliers", "presso_distortion_p",
    "f_min", "f_median", "weak_ivs", "loo_or_min", "loo_or_max", "direction_std"
  )]
  write.csv(csv_df, out_path, row.names = FALSE)
  latex_df <- combined_df[, c("Direction", "Method", "N_SNPs", "OR", "OR_lci", "OR_uci", "P", "egger_intercept_p", "presso_global_p", "steiger_p", "I2", "loo_or_min", "loo_or_max")]
  numeric_cols <- vapply(latex_df, is.numeric, logical(1))
  latex_df[numeric_cols] <- lapply(latex_df[numeric_cols], function(x) sprintf("%.3g", x))
  latex_path <- file.path(context$tables_dir, "Table_MR_Sensitivity_Summary.tex")
  simple_latex_table(latex_df, latex_path, caption = "MR sensitivity summary including directionality, pleiotropy, heterogeneity, and leave-one-out ranges.", label = "tab:rev1_mr")

  plot_df <- combined_df[, c("Direction", "Method", "OR", "OR_lci", "OR_uci")]
  plot_df$Method <- factor(plot_df$Method, levels = rev(unique(plot_df$Method)))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Method, y = OR, ymin = OR_lci, ymax = OR_uci, color = Direction)) +
    ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "grey60") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "MR sensitivity summary", x = NULL, y = "Odds ratio (log scale)", color = "Direction") +
    ggplot2::theme_bw(base_size = 11)
  figure_path <- file.path(context$figures_supp_dir, "rev1_mr_sensitivity_summary.pdf")
  ggplot2::ggsave(figure_path, plot = p, width = 8.5, height = 5.5)

  log_message(context, step, sprintf("MR sensitivity table written: %s", out_path))
  log_message(context, step, sprintf("MR sensitivity figure written: %s", figure_path))
}
