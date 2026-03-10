run_step <- function(context) {
  source(file.path(context$revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)
  config <- load_revision_config(context)
  step <- "02_common_genes_report"
  log_message(context, step, "Computing IHS specificity summaries")

  ihs_genes <- config$ihs$genes

  infer_binary_groups <- function(group_values, candidate, text_values = NULL) {
    v <- tolower(trimws(as.character(group_values)))
    case_hit <- grepl(candidate$case_regex, v, ignore.case = TRUE, perl = TRUE) | v %in% c("disease", "case")
    control_hit <- grepl(candidate$control_regex, v, ignore.case = TRUE, perl = TRUE) | v %in% c("control", "healthy", "normal")

    if (!is.null(text_values)) {
      txt <- tolower(trimws(as.character(text_values)))
      case_hit <- case_hit | grepl(candidate$case_regex, txt, ignore.case = TRUE, perl = TRUE)
      control_hit <- control_hit | grepl(candidate$control_regex, txt, ignore.case = TRUE, perl = TRUE)
    }
    ifelse(case_hit & !control_hit, "Disease", ifelse(control_hit & !case_hit, "Control", NA_character_))
  }

  infer_geo_groups_from_pheno <- function(pheno_df, candidate) {
    col_priority <- c(
      "title", "source_name_ch1", "description", "characteristics_ch1", "characteristics_ch1.1",
      "disease state", "disease_state", "disease", "status", "group", "condition"
    )
    col_priority <- unique(col_priority[col_priority %in% colnames(pheno_df)])
    if (length(col_priority) == 0L) {
      col_priority <- colnames(pheno_df)
    }
    text_mat <- lapply(col_priority, function(col) as.character(pheno_df[[col]]))
    joined <- apply(do.call(cbind, text_mat), 1, function(x) paste(x, collapse = " | "))
    group <- infer_binary_groups(joined, candidate)
    data.frame(group = group, label_text = joined, source_columns = paste(col_priority, collapse = ";"), stringsAsFactors = FALSE)
  }

  load_negative_control_cohort <- function(candidate, manual_root) {
    cohort_dir <- file.path(manual_root, candidate$manual_subdir)
    expr_path <- file.path(cohort_dir, "expression_matrix.csv")
    anno_path <- file.path(cohort_dir, "sample_annotation.csv")
    series_path <- file.path(cohort_dir, paste0(candidate$accession, "_series_matrix.txt.gz"))

    if (file.exists(expr_path) && file.exists(anno_path)) {
      expr <- read_expression_matrix(expr_path)
      pheno <- read_table_file(anno_path)
      sample_col_candidates <- c("sample_id", "Sample_ID", "Sample", "sample", "gsm", "GSM", "geo_accession")
      group_col_candidates <- c("group", "Group", "condition", "Condition", "status", "Status")
      sample_col <- sample_col_candidates[sample_col_candidates %in% colnames(pheno)][1]
      group_col <- group_col_candidates[group_col_candidates %in% colnames(pheno)][1]
      if (is.na(sample_col) || is.na(group_col)) {
        return(list(ok = FALSE, reason = "sample_annotation_missing_required_columns"))
      }

      pheno$sample_id_resolved <- as.character(pheno[[sample_col]])
      pheno$group_binary <- infer_binary_groups(pheno[[group_col]], candidate, text_values = pheno[["label_text"]])
      pheno <- pheno[!is.na(pheno$group_binary), , drop = FALSE]
      if (nrow(pheno) == 0L) {
        return(list(ok = FALSE, reason = "group_mapping_failed"))
      }
      pheno$GroupBinary <- pheno$group_binary
      pheno$Sample_ID <- pheno$sample_id_resolved
      return(list(
        ok = TRUE,
        source = "manual_matrix_annotation",
        expr = expr,
        pheno = pheno,
        sample_col = "Sample_ID",
        group_col = "GroupBinary",
        positive_label = "Disease"
      ))
    }

    if (file.exists(series_path)) {
      # Default to strict offline behavior: a probe-level series matrix alone is
      # not considered analyzable unless users explicitly opt in to GEOquery parsing.
      if (isTRUE(getOption("rev1.allow_geoquery_series_parse", FALSE)) &&
          requireNamespace("GEOquery", quietly = TRUE)) {
        series_obj <- tryCatch(GEOquery::getGEO(filename = series_path), error = function(e) e)
        if (inherits(series_obj, "error")) {
          return(list(ok = FALSE, reason = "series_matrix_parse_failed"))
        }
        expr <- Biobase::exprs(series_obj)
        pheno <- Biobase::pData(series_obj)
        sample_id <- if ("geo_accession" %in% colnames(pheno)) as.character(pheno$geo_accession) else rownames(pheno)
        pheno$Sample_ID <- sample_id
        group_df <- infer_geo_groups_from_pheno(pheno, candidate)
        pheno$GroupBinary <- group_df$group
        pheno$label_text <- group_df$label_text
        pheno$label_source_columns <- group_df$source_columns
        pheno <- pheno[!is.na(pheno$GroupBinary), , drop = FALSE]
        expr <- expr[, pheno$Sample_ID, drop = FALSE]

        expr_df <- data.frame(gene = rownames(expr), expr, check.names = FALSE, stringsAsFactors = FALSE)
        write.csv(expr_df, file.path(cohort_dir, "expression_matrix.csv"), row.names = FALSE)
        write.csv(
          pheno[, c("Sample_ID", "GroupBinary", "label_text", "label_source_columns"), drop = FALSE],
          file.path(cohort_dir, "sample_annotation.csv"),
          row.names = FALSE
        )

        return(list(
          ok = TRUE,
          source = "series_matrix_parsed",
          expr = expr,
          pheno = pheno,
          sample_col = "Sample_ID",
          group_col = "GroupBinary",
          positive_label = "Disease"
        ))
      }
      return(list(ok = FALSE, reason = "series_matrix_present_probe_level_manual_conversion_required"))
    }

    list(ok = FALSE, reason = "manual_files_not_available")
  }

  summarise_ihs <- function(dataset_id, dataset_label, cohort_type, expr, pheno, sample_col, group_col, positive_label) {
    aligned <- align_expression_with_pheno(expr, pheno, sample_col, group_col, positive_label)
    score <- compute_ihs(aligned$expr, ihs_genes)
    score_sd <- stats::sd(score)
    score_z <- if (is.finite(score_sd) && score_sd > 0) (score - mean(score)) / score_sd else rep(0, length(score))
    auc_ci <- roc_auc_ci(aligned$y, score)
    glm_fit <- tryCatch(glm(aligned$y ~ score_z, family = binomial()), error = function(e) NULL)
    cohen_d <- tryCatch({
      case <- score[aligned$y == 1]
      control <- score[aligned$y == 0]
      pooled_sd <- sqrt(((length(case) - 1) * stats::var(case) + (length(control) - 1) * stats::var(control)) / (length(case) + length(control) - 2))
      (mean(case) - mean(control)) / pooled_sd
    }, error = function(e) NA_real_)
    sample_df <- data.frame(
      dataset_id = dataset_id,
      dataset_label = dataset_label,
      cohort_type = cohort_type,
      sample_id = aligned$sample_ids,
      group = as.character(aligned$pheno[[group_col]]),
      y = aligned$y,
      ihs_score = score,
      stringsAsFactors = FALSE
    )
    summary_df <- data.frame(
      dataset_id = dataset_id,
      dataset_label = dataset_label,
      cohort_type = cohort_type,
      n_case = sum(aligned$y == 1),
      n_control = sum(aligned$y == 0),
      ihs_case_mean = safe_mean(score[aligned$y == 1]),
      ihs_control_mean = safe_mean(score[aligned$y == 0]),
      cohen_d = cohen_d,
      auc = auc_ci["auc"],
      auc_ci_low = auc_ci["ci_low"],
      auc_ci_high = auc_ci["ci_high"],
      or_per_sd = if (is.null(glm_fit)) NA_real_ else unname(exp(coef(glm_fit)[2])),
      p_value = if (is.null(glm_fit)) NA_real_ else summary(glm_fit)$coefficients[2, 4],
      stringsAsFactors = FALSE
    )
    list(samples = sample_df, summary = summary_df)
  }

  train_expr <- read_expression_matrix(resolve_revision_path(context, config$paths$train_expr))
  train_pheno <- read_table_file(resolve_revision_path(context, config$paths$train_pheno))
  ext_expr <- read_expression_matrix(resolve_revision_path(context, config$paths$external_expr))
  ext_pheno <- read_table_file(resolve_revision_path(context, config$paths$external_pheno))

  target_train <- summarise_ihs(
    dataset_id = config$cohorts$train$id,
    dataset_label = config$cohorts$train$label,
    cohort_type = "target",
    expr = train_expr,
    pheno = train_pheno,
    sample_col = config$cohorts$train$sample_col,
    group_col = config$cohorts$train$group_col,
    positive_label = config$cohorts$train$positive_label
  )
  target_external <- summarise_ihs(
    dataset_id = config$cohorts$external$id,
    dataset_label = config$cohorts$external$label,
    cohort_type = "target",
    expr = ext_expr,
    pheno = ext_pheno,
    sample_col = config$cohorts$external$sample_col,
    group_col = config$cohorts$external$group_col,
    positive_label = config$cohorts$external$positive_label
  )

  manual_root <- resolve_revision_path(context, config$paths$manual_negative_control_dir)
  negative_summaries <- list()
  descriptive_only_summaries <- list()
  negative_samples <- list()
  status_rows <- list()

  for (candidate in config$negative_controls$candidate_gse) {
    loaded <- load_negative_control_cohort(candidate, manual_root)
    if (isTRUE(loaded$ok)) {
      res <- summarise_ihs(
        dataset_id = candidate$accession,
        dataset_label = candidate$title,
        cohort_type = "negative_control",
        expr = loaded$expr,
        pheno = loaded$pheno,
        sample_col = loaded$sample_col,
        group_col = loaded$group_col,
        positive_label = loaded$positive_label
      )
      negative_samples[[length(negative_samples) + 1L]] <- res$samples
      n_case <- as.integer(res$summary$n_case[1])
      n_control <- as.integer(res$summary$n_control[1])

      if (n_case > 0L && n_control > 0L) {
        negative_summaries[[length(negative_summaries) + 1L]] <- res$summary
        status_rows[[length(status_rows) + 1L]] <- data.frame(
          accession = candidate$accession,
          status = "analyzable",
          source = loaded$source,
          n_case = n_case,
          n_control = n_control,
          stringsAsFactors = FALSE
        )
      } else if (n_case > 0L && n_control == 0L) {
        descriptive_only_summaries[[length(descriptive_only_summaries) + 1L]] <- res$summary
        status_rows[[length(status_rows) + 1L]] <- data.frame(
          accession = candidate$accession,
          status = "descriptive_only_no_control",
          source = loaded$source,
          n_case = n_case,
          n_control = n_control,
          stringsAsFactors = FALSE
        )
      } else if (n_case == 0L && n_control > 0L) {
        descriptive_only_summaries[[length(descriptive_only_summaries) + 1L]] <- res$summary
        status_rows[[length(status_rows) + 1L]] <- data.frame(
          accession = candidate$accession,
          status = "descriptive_only_no_disease",
          source = loaded$source,
          n_case = n_case,
          n_control = n_control,
          stringsAsFactors = FALSE
        )
      } else {
        status_rows[[length(status_rows) + 1L]] <- data.frame(
          accession = candidate$accession,
          status = "group_mapping_failed",
          source = loaded$source,
          n_case = n_case,
          n_control = n_control,
          stringsAsFactors = FALSE
        )
      }
    } else {
      status_rows[[length(status_rows) + 1L]] <- data.frame(
        accession = candidate$accession,
        status = loaded$reason,
        source = NA_character_,
        n_case = NA_integer_,
        n_control = NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- rbind(
    target_train$summary,
    target_external$summary,
    if (length(negative_summaries) > 0L) do.call(rbind, negative_summaries) else NULL
  )
  score_df <- rbind(
    target_train$samples,
    target_external$samples,
    if (length(negative_samples) > 0L) do.call(rbind, negative_samples) else NULL
  )

  summary_path <- file.path(context$tables_dir, "Table_NegControl_IHS_Performance_Comparison.csv")
  write.csv(summary_df, summary_path, row.names = FALSE)
  descriptive_path <- file.path(context$tables_dir, "Table_NegControl_DescriptiveOnly.csv")
  if (length(descriptive_only_summaries) > 0L) {
    write.csv(do.call(rbind, descriptive_only_summaries), descriptive_path, row.names = FALSE)
  } else {
    write.csv(data.frame(), descriptive_path, row.names = FALSE)
  }
  score_path <- file.path(context$tables_dir, "IHS_target_negative_control_scores.csv")
  write.csv(score_df, score_path, row.names = FALSE)
  status_df <- do.call(rbind, status_rows)
  status_path <- file.path(context$tables_dir, "negative_control_analysis_status.csv")
  write.csv(status_df, status_path, row.names = FALSE)

  latex_df <- summary_df
  numeric_cols <- vapply(latex_df, is.numeric, logical(1))
  latex_df[numeric_cols] <- lapply(latex_df[numeric_cols], function(x) sprintf("%.3f", x))
  colnames(latex_df) <- c("dataset", "label", "cohort", "nCase", "nControl", "meanCase", "meanControl", "cohenD", "AUC", "AUCloCI", "AUChiCI", "ORperSD", "pValue")
  latex_path <- file.path(context$tables_dir, "Table_NegControl_IHS_Performance_Comparison.tex")
  simple_latex_table(
    latex_df,
    latex_path,
    caption = "Auxiliary reviewer-boundary score performance in target MASLD cohorts and available ALD/AH negative-control cohorts. These values are boundary-setting contextual checks from the negative-control module, not the primary locked-model performance estimates reported in Fig.~4, and should not be interpreted as replacing the main 10-gene claim.",
    label = "tab:rev1_negcontrol",
    resize_to_textwidth = TRUE
  )

  # Additional transparency outputs for descriptive-only cohorts and status audit.
  descriptive_ids <- status_df$accession[status_df$status %in% c("descriptive_only_no_control", "descriptive_only_no_disease")]
  descriptive_dist_path <- file.path(context$tables_dir, "Table_NegControl_DescriptiveOnly_Distribution.csv")
  if (length(descriptive_ids) > 0L) {
    dist_rows <- list()
    for (dataset_id in unique(descriptive_ids)) {
      sub <- score_df[score_df$dataset_id == dataset_id, , drop = FALSE]
      if (nrow(sub) == 0L) next
      status_val <- status_df$status[match(dataset_id, status_df$accession)]
      qv <- stats::quantile(sub$ihs_score, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
      dist_rows[[length(dist_rows) + 1L]] <- data.frame(
        dataset_id = dataset_id,
        status = status_val,
        n_samples = nrow(sub),
        available_group_levels = paste(sort(unique(as.character(sub$group))), collapse = ";"),
        ihs_mean = safe_mean(sub$ihs_score),
        ihs_sd = stats::sd(sub$ihs_score, na.rm = TRUE),
        ihs_min = suppressWarnings(min(sub$ihs_score, na.rm = TRUE)),
        ihs_p25 = qv[1],
        ihs_median = qv[2],
        ihs_p75 = qv[3],
        ihs_max = suppressWarnings(max(sub$ihs_score, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }
    descriptive_dist_df <- do.call(rbind, dist_rows)
    write.csv(descriptive_dist_df, descriptive_dist_path, row.names = FALSE)
    descriptive_latex <- descriptive_dist_df
    desc_num <- vapply(descriptive_latex, is.numeric, logical(1))
    descriptive_latex[desc_num] <- lapply(descriptive_latex[desc_num], function(x) sprintf("%.3f", x))
    colnames(descriptive_latex) <- c("dataset", "status", "n", "groups", "mean", "sd", "min", "q25", "median", "q75", "max")
    descriptive_latex_path <- file.path(context$tables_dir, "Table_NegControl_DescriptiveOnly_Distribution.tex")
    simple_latex_table(
      descriptive_latex,
      descriptive_latex_path,
      caption = "Descriptive-only negative-control cohorts (single-class) with auxiliary reviewer-boundary score distribution summaries.",
      label = "tab:rev1_negcontrol_descriptive",
      resize_to_textwidth = TRUE
    )
  } else {
    write.csv(data.frame(), descriptive_dist_path, row.names = FALSE)
  }

  status_latex <- status_df
  colnames(status_latex) <- c("accession", "status", "source", "nCase", "nControl")
  status_latex_path <- file.path(context$tables_dir, "Table_NegControl_Analysis_Status.tex")
  simple_latex_table(
    status_latex,
    status_latex_path,
    caption = "Negative-control analysis status by cohort (analyzable vs descriptive-only vs pending).",
    label = "tab:rev1_negcontrol_status"
  )

  pending <- status_df
  no_control_ids <- pending$accession[pending$status == "descriptive_only_no_control"]
  no_disease_ids <- pending$accession[pending$status == "descriptive_only_no_disease"]
  pending_ids <- pending$accession[!(pending$status %in% c("analyzable", "descriptive_only_no_control", "descriptive_only_no_disease"))]
  subtitle_parts <- character(0)
  if (length(no_control_ids) > 0L) {
    subtitle_parts <- c(subtitle_parts, paste("descriptive-only cohort with no controls:", paste(no_control_ids, collapse = ", ")))
  }
  if (length(no_disease_ids) > 0L) {
    subtitle_parts <- c(subtitle_parts, paste("descriptive-only cohort with no disease cases:", paste(no_disease_ids, collapse = ", ")))
  }
  if (length(pending_ids) > 0L) {
    subtitle_parts <- c(subtitle_parts, paste("Pending manual negative-control files:", paste(pending_ids, collapse = ", ")))
  }
  subtitle_text <- if (length(subtitle_parts) > 0L) {
    paste(subtitle_parts, collapse = " | ")
  } else {
    "Available target and ALD/AH negative-control cohorts"
  }

  boxplot_data <- score_df
  boxplot_data$group <- factor(boxplot_data$group)
  facet_levels <- unique(c(as.character(summary_df$dataset_label), as.character(score_df$dataset_label)))
  facet_levels <- facet_levels[!is.na(facet_levels) & nzchar(facet_levels)]
  boxplot_data$dataset_label <- factor(boxplot_data$dataset_label, levels = facet_levels)
  plot_seed <- suppressWarnings(as.integer(config$revision$seed))
  if (!is.finite(plot_seed)) {
    plot_seed <- 20260306L
  }
  jitter_pos <- ggplot2::position_jitter(width = 0.1, height = 0, seed = plot_seed)

  p <- ggplot2::ggplot(boxplot_data, ggplot2::aes(x = group, y = ihs_score, fill = group)) +
    ggplot2::geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.8) +
    ggplot2::geom_jitter(position = jitter_pos, alpha = 0.6, size = 1.5) +
    ggplot2::facet_wrap(~ dataset_label, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("Control" = "#4F6D7A", "NAFLD" = "#C44E52", "Disease" = "#C44E52")) +
    ggplot2::labs(
      title = "Auxiliary reviewer-boundary score across target MASLD and ALD/AH cohorts",
      subtitle = subtitle_text,
      x = NULL,
      y = "Auxiliary score"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(face = "bold"))

  figure_path <- file.path(context$figures_supp_dir, "rev1_negative_control_specificity.pdf")
  ggplot2::ggsave(figure_path, plot = p, width = 10, height = 6)

  interp_lines <- c(
    sprintf("Target-train AUC: %.3f", target_train$summary$auc),
    sprintf("Target-external AUC: %.3f", target_external$summary$auc)
  )
  if (length(descriptive_only_summaries) > 0L) {
    desc_df <- do.call(rbind, descriptive_only_summaries)
    for (i in seq_len(nrow(desc_df))) {
      interp_lines <- c(
        interp_lines,
        sprintf(
          "Descriptive-only cohort (no discrimination metrics reported): %s | n_case=%d | n_control=%d",
          desc_df$dataset_id[i],
          as.integer(desc_df$n_case[i]),
          as.integer(desc_df$n_control[i])
        )
      )
    }
  }
  if (length(negative_summaries) == 0L) {
    interp_lines <- c(
      interp_lines,
      "No eligible ALD/AH negative-control cohort was analyzable from local or auto-downloaded files.",
      "Specificity remains unresolved until manual ALD/AH cohort files are placed under rev1_major_revision/manual_data/negative_controls/.",
      "Suggested manuscript downgrade: current evidence supports within-MASLD transportability but does not establish disease specificity versus non-MASLD inflammatory liver disease."
    )
  }
  interp_path <- file.path(context$tables_dir, "negative_control_interpretation.txt")
  write_text_file(interp_path, interp_lines)

  log_message(context, step, sprintf("Summary table written: %s", summary_path))
  log_message(context, step, sprintf("Specificity figure written: %s", figure_path))
}
