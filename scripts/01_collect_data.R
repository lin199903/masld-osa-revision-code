run_step <- function(context) {
  source(file.path(context$revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)
  config <- load_revision_config(context)
  step <- "01_collect_data"
  log_message(context, step, "Loading negative-control configuration (ALD/AH priority cohorts)")

  manual_root <- ensure_dir(resolve_revision_path(context, config$paths$manual_negative_control_dir))
  cache_root <- ensure_dir(resolve_revision_path(context, config$paths$geo_cache_dir))

  audit_rows <- lapply(config$negative_controls$local_misannotated_series, function(item) {
    local_file <- resolve_project_path(context, paste0(item$accession, "_series_matrix.txt.gz"))
    data.frame(
      accession = item$accession,
      local_file_present = file.exists(local_file),
      actual_annotation = item$actual_annotation,
      eligible_as_negative_control = FALSE,
      stringsAsFactors = FALSE
    )
  })
  audit_df <- do.call(rbind, audit_rows)
  audit_path <- file.path(context$tables_dir, "negative_control_local_audit.csv")
  write.csv(audit_df, audit_path, row.names = FALSE)
  log_message(context, step, sprintf("Wrote local audit: %s", audit_path))

  build_geo_url <- function(accession) {
    prefix <- sub("[0-9]{3}$", "nnn", accession)
    sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz", prefix, accession, accession)
  }

  infer_group_vector <- function(pheno_df, candidate) {
    column_priority <- c(
      "title", "source_name_ch1", "description", "characteristics_ch1", "characteristics_ch1.1",
      "disease state", "disease_state", "disease", "status", "group", "condition"
    )
    column_priority <- unique(column_priority[column_priority %in% colnames(pheno_df)])
    if (length(column_priority) == 0L) {
      column_priority <- colnames(pheno_df)
    }

    text_mat <- lapply(column_priority, function(col) {
      as.character(pheno_df[[col]])
    })
    names(text_mat) <- column_priority
    joined_text <- apply(do.call(cbind, text_mat), 1, function(x) paste(x, collapse = " | "))

    case_regex <- candidate$case_regex
    control_regex <- candidate$control_regex
    is_case <- grepl(case_regex, joined_text, ignore.case = TRUE, perl = TRUE)
    is_control <- grepl(control_regex, joined_text, ignore.case = TRUE, perl = TRUE)

    group <- ifelse(is_case & !is_control, "Disease", ifelse(is_control & !is_case, "Control", NA_character_))
    source_col <- rep(paste(column_priority, collapse = ";"), length(group))

    list(group = group, joined_text = joined_text, source_col = source_col)
  }

  try_auto_geoquery <- function(accession, candidate, cache_dir, manual_dir) {
    if (!requireNamespace("GEOquery", quietly = TRUE)) {
      return(list(
        status = "geoquery_not_installed",
        message = "GEOquery is not installed; skipped automatic download.",
        n_case = NA_integer_,
        n_control = NA_integer_,
        files_ready = FALSE
      ))
    }

    gse_obj <- tryCatch(
      GEOquery::getGEO(accession, GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE, destdir = cache_dir),
      error = function(e) e
    )
    if (inherits(gse_obj, "error")) {
      return(list(
        status = "geoquery_failed",
        message = conditionMessage(gse_obj),
        n_case = NA_integer_,
        n_control = NA_integer_,
        files_ready = FALSE
      ))
    }

    if (is.list(gse_obj)) {
      n_samples <- vapply(gse_obj, function(x) ncol(Biobase::exprs(x)), numeric(1))
      keep_idx <- which.max(n_samples)
      eset <- gse_obj[[keep_idx]]
      platform_msg <- names(gse_obj)[keep_idx]
    } else {
      eset <- gse_obj
      platform_msg <- "single_platform"
    }

    expr <- Biobase::exprs(eset)
    pheno <- Biobase::pData(eset)
    pheno$sample_id <- if ("geo_accession" %in% colnames(pheno)) as.character(pheno$geo_accession) else rownames(pheno)

    expr_sample_ids <- colnames(expr)
    match_idx <- match(expr_sample_ids, pheno$sample_id)
    keep <- !is.na(match_idx)
    expr <- expr[, keep, drop = FALSE]
    pheno <- pheno[match_idx[keep], , drop = FALSE]
    rownames(pheno) <- pheno$sample_id

    group_info <- infer_group_vector(pheno, candidate)
    anno_df <- data.frame(
      sample_id = colnames(expr),
      group = group_info$group,
      label_text = group_info$joined_text,
      label_source_columns = group_info$source_col,
      stringsAsFactors = FALSE
    )
    anno_df <- anno_df[!is.na(anno_df$group), , drop = FALSE]
    expr <- expr[, anno_df$sample_id, drop = FALSE]

    n_case <- sum(anno_df$group == "Disease")
    n_control <- sum(anno_df$group == "Control")
    if (n_case == 0L || n_control == 0L) {
      return(list(
        status = "geoquery_downloaded_group_unresolved",
        message = sprintf("Downloaded via GEOquery (%s) but failed to resolve both Disease and Control from metadata.", platform_msg),
        n_case = n_case,
        n_control = n_control,
        files_ready = FALSE
      ))
    }

    expr_df <- data.frame(gene = rownames(expr), expr, check.names = FALSE, stringsAsFactors = FALSE)
    write.csv(expr_df, file.path(manual_dir, "expression_matrix.csv"), row.names = FALSE)
    write.csv(anno_df, file.path(manual_dir, "sample_annotation.csv"), row.names = FALSE)
    saveRDS(eset, file.path(cache_dir, paste0(accession, "_eset.rds")))

    list(
      status = "geoquery_downloaded_and_parsed",
      message = sprintf("Downloaded and parsed via GEOquery (%s).", platform_msg),
      n_case = n_case,
      n_control = n_control,
      files_ready = TRUE
    )
  }

  candidate_rows <- list()
  manual_lines <- c(
    "Candidate GSE list for reviewer-response negative-control liver cohorts:",
    "",
    sprintf("Manual placement root: %s", normalizePath(manual_root, winslash = "/", mustWork = FALSE)),
    "Primary negative controls for Revision 1.1: GSE28619 and GSE103580 (ALD/AH cohorts).",
    "Place each cohort under rev1_major_revision/manual_data/negative_controls/<GSE>/.",
    "Accepted file names per cohort directory:",
    "  1) <GSE>_series_matrix.txt.gz",
    "  2) expression_matrix.csv + sample_annotation.csv",
    "Sample annotation must contain a binary disease label (Disease vs Control) or equivalent text metadata.",
    ""
  )

  for (candidate in config$negative_controls$candidate_gse) {
    accession <- candidate$accession
    manual_dir <- ensure_dir(file.path(manual_root, candidate$manual_subdir))
    cache_dir <- ensure_dir(file.path(cache_root, accession))

    local_series <- file.path(manual_dir, paste0(accession, "_series_matrix.txt.gz"))
    local_expr <- file.path(manual_dir, "expression_matrix.csv")
    local_anno <- file.path(manual_dir, "sample_annotation.csv")
    local_ready <- file.exists(local_series) || (file.exists(local_expr) && file.exists(local_anno))

    auto_result <- list(
      status = if (local_ready) "manual_files_already_present" else "not_attempted",
      message = if (local_ready) "Manual files already present; skipped GEOquery." else "Will attempt GEOquery download.",
      n_case = NA_integer_,
      n_control = NA_integer_,
      files_ready = local_ready
    )

    if (!local_ready) {
      auto_result <- try_auto_geoquery(accession, candidate, cache_dir, manual_dir)
      log_message(
        context,
        step,
        sprintf("[%s] Auto status: %s | %s", accession, auto_result$status, auto_result$message),
        level = if (grepl("failed|unresolved|not_installed", auto_result$status)) "WARNING" else "INFO"
      )
    } else {
      log_message(context, step, sprintf("[%s] Manual files already present.", accession))
    }

    final_ready <- file.exists(local_series) || (file.exists(local_expr) && file.exists(local_anno))
    candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
      accession = accession,
      disease = candidate$disease,
      title = candidate$title,
      priority = candidate$priority,
      geo_series_url = build_geo_url(accession),
      cache_dir = normalizePath(cache_dir, winslash = "/", mustWork = FALSE),
      local_manual_dir = normalizePath(manual_dir, winslash = "/", mustWork = FALSE),
      download_status = auto_result$status,
      download_message = auto_result$message,
      auto_n_case = auto_result$n_case,
      auto_n_control = auto_result$n_control,
      files_ready = final_ready,
      stringsAsFactors = FALSE
    )

    manual_lines <- c(
      manual_lines,
      sprintf("- %s | %s | %s", accession, candidate$disease, candidate$title),
      sprintf("  Priority: %s", as.character(candidate$priority)),
      sprintf("  Manual target: rev1_major_revision/manual_data/negative_controls/%s/", candidate$manual_subdir),
      sprintf("  GEO matrix URL attempted: %s", build_geo_url(accession)),
      sprintf("  Auto status: %s", auto_result$status),
      sprintf("  Auto message: %s", auto_result$message),
      "  Required for manual fallback: <GSE>_series_matrix.txt.gz OR expression_matrix.csv + sample_annotation.csv",
      ""
    )
  }

  candidate_df <- do.call(rbind, candidate_rows)
  candidate_df <- candidate_df[order(candidate_df$priority), , drop = FALSE]
  candidate_path <- file.path(context$tables_dir, "negative_control_candidates.csv")
  write.csv(candidate_df, candidate_path, row.names = FALSE)

  manual_path <- file.path(context$tables_dir, "negative_control_manual_download_instructions.txt")
  write_text_file(manual_path, manual_lines)

  missing_count <- sum(!candidate_df$files_ready)
  if (missing_count > 0L) {
    message("Automatic negative-control acquisition incomplete; manual placement instructions generated.")
    message(sprintf("Manual placement root: %s", normalizePath(manual_root, winslash = "/", mustWork = FALSE)))
  }

  log_message(context, step, sprintf("Candidate table written: %s", candidate_path))
  log_message(context, step, sprintf("Manual instructions written: %s", manual_path))
}
