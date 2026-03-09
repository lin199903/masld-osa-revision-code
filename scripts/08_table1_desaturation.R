run_step <- function(context) {
  source(file.path(context$revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)
  config <- load_revision_config(context)
  step <- "08_table1_desaturation"
  log_message(context, step, "Extracting hypoxia and desaturation metadata")

  read_gz_lines <- function(path) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    readLines(con, warn = FALSE)
  }

  extract_geo_accessions <- function(lines) {
    line <- grep("^!Sample_geo_accession", lines, value = TRUE)[1]
    if (is.na(line)) {
      return(character())
    }
    parts <- strsplit(line, "\t")[[1]]
    gsub('"', "", parts[-1])
  }

  extract_numeric_characteristic <- function(lines, pattern) {
    line <- grep(pattern, lines, ignore.case = TRUE, value = TRUE)[1]
    if (is.na(line)) {
      return(NULL)
    }
    parts <- strsplit(line, "\t")[[1]][-1]
    as.numeric(sub(".*: *", "", gsub('"', "", parts)))
  }

  rows <- list()

  gse75097_lines <- read_gz_lines(resolve_revision_path(context, config$paths$gse75097_series))
  gse75097_accessions <- extract_geo_accessions(gse75097_lines)
  gse75097_ahi <- extract_numeric_characteristic(gse75097_lines, "apnea hyponea index")
  gse75097_pdata <- read.csv(resolve_revision_path(context, config$paths$gse75097_pdata), check.names = FALSE)
  if (!is.null(gse75097_ahi) && length(gse75097_accessions) == length(gse75097_ahi)) {
    ahi_df <- data.frame(Sample = gse75097_accessions, AHI_events_per_h = gse75097_ahi, stringsAsFactors = FALSE)
    merged <- merge(gse75097_pdata, ahi_df, by = "Sample", all.x = TRUE)
    for (grp in unique(merged$Group)) {
      sub <- merged[merged$Group == grp, , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = "GSE75097",
        tissue = "PBMC",
        group = grp,
        n = nrow(sub),
        AHI_events_per_h = safe_mean(sub$AHI_events_per_h),
        RDI_events_per_h = NA_real_,
        ODI_events_per_h = NA_real_,
        T90 = NA_real_,
        mean_SpO2_pct = NA_real_,
        nadir_SpO2_pct = NA_real_,
        source = "Sample-level GEO characteristic",
        limitation_note = "Only AHI was available in GEO series_matrix metadata.",
        stringsAsFactors = FALSE
      )
    }
  }

  gse38792_lines <- read_gz_lines(resolve_revision_path(context, config$paths$gse38792_series))
  series_summary <- grep("^!Series_summary", gse38792_lines, value = TRUE)
  summary_text <- paste(series_summary, collapse = " ")
  rdi_vals <- regmatches(summary_text, regexec("respiratory disturbance index \\(([-0-9.]+) vs\\. ([-0-9.]+)", summary_text, ignore.case = TRUE))[[1]]
  nadir_vals <- regmatches(summary_text, regexec("minimum oxygen saturation ([-0-9.]+)% vs\\. ([-0-9.]+)%", summary_text, ignore.case = TRUE))[[1]]
  if (length(rdi_vals) == 3L) {
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = "GSE38792",
      tissue = "Visceral adipose",
      group = "OSA",
      n = sum(read.csv(resolve_revision_path(context, config$paths$gse38792_pdata), check.names = FALSE)$Group == "OSA"),
      AHI_events_per_h = NA_real_,
      RDI_events_per_h = as.numeric(rdi_vals[2]),
      ODI_events_per_h = NA_real_,
      T90 = NA_real_,
      mean_SpO2_pct = NA_real_,
      nadir_SpO2_pct = if (length(nadir_vals) == 3L) as.numeric(nadir_vals[2]) else NA_real_,
      source = "Series summary group-level text",
      limitation_note = "Series summary reported RDI and minimum oxygen saturation, not sample-level AHI/ODI/T90.",
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = "GSE38792",
      tissue = "Visceral adipose",
      group = "Control",
      n = sum(read.csv(resolve_revision_path(context, config$paths$gse38792_pdata), check.names = FALSE)$Group == "Control"),
      AHI_events_per_h = NA_real_,
      RDI_events_per_h = as.numeric(rdi_vals[3]),
      ODI_events_per_h = NA_real_,
      T90 = NA_real_,
      mean_SpO2_pct = NA_real_,
      nadir_SpO2_pct = if (length(nadir_vals) == 3L) as.numeric(nadir_vals[3]) else NA_real_,
      source = "Series summary group-level text",
      limitation_note = "Series summary reported RDI and minimum oxygen saturation, not sample-level AHI/ODI/T90.",
      stringsAsFactors = FALSE
    )
  }

  for (cohort_id in c("GSE135917", "GSE89632", "GSE135251")) {
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = cohort_id,
      tissue = ifelse(cohort_id == "GSE135917", "Subcutaneous adipose", "Liver"),
      group = "NA",
      n = NA_integer_,
      AHI_events_per_h = NA_real_,
      RDI_events_per_h = NA_real_,
      ODI_events_per_h = NA_real_,
      T90 = NA_real_,
      mean_SpO2_pct = NA_real_,
      nadir_SpO2_pct = NA_real_,
      source = "No usable desaturation variables in local series_matrix metadata",
      limitation_note = "Public transcriptomic metadata did not contain extractable ODI/T90/mean SpO2/nadir SpO2/AHI fields.",
      stringsAsFactors = FALSE
    )
  }

  hypoxia_df <- do.call(rbind, rows)
  out_csv <- file.path(context$tables_dir, "Table1_Supplement_Hypoxia_Parameters_wide.csv")
  write.csv(hypoxia_df, out_csv, row.names = FALSE)
  out_alias <- file.path(context$tables_dir, "Table1_Supplement_Hypoxia_Parameters.csv")
  write.csv(hypoxia_df, out_alias, row.names = FALSE)

  latex_df <- hypoxia_df
  numeric_cols <- vapply(latex_df, is.numeric, logical(1))
  latex_df[numeric_cols] <- lapply(latex_df[numeric_cols], function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x)))
  latex_path <- file.path(context$tables_dir, "Table1_Supplement_Hypoxia_Parameters.tex")
  simple_latex_table(latex_df, latex_path, caption = "Available hypoxia/desaturation parameters extracted from local public metadata. NA indicates metadata not available in the local repository.", label = "tab:rev1_hypoxia")

  metric_cols <- c("AHI_events_per_h", "RDI_events_per_h", "ODI_events_per_h", "T90", "mean_SpO2_pct", "nadir_SpO2_pct")
  completeness_rows <- lapply(metric_cols, function(metric) {
    values <- hypoxia_df[[metric]]
    available <- !is.na(values)
    available_cohorts <- unique(hypoxia_df$cohort[available])
    data.frame(
      metric = metric,
      non_missing_rows = sum(available),
      total_rows = length(values),
      coverage_pct = 100 * sum(available) / length(values),
      available_cohorts = if (length(available_cohorts) == 0L) "none" else paste(available_cohorts, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  completeness_df <- do.call(rbind, completeness_rows)
  completeness_csv <- file.path(context$tables_dir, "Table1_Hypoxia_Metadata_Completeness.csv")
  write.csv(completeness_df, completeness_csv, row.names = FALSE)
  completeness_latex <- completeness_df
  completeness_latex$coverage_pct <- sprintf("%.1f", completeness_latex$coverage_pct)
  colnames(completeness_latex) <- c("metric", "nNonMissing", "nTotal", "coveragePct", "availableCohorts")
  completeness_tex <- file.path(context$tables_dir, "Table1_Hypoxia_Metadata_Completeness.tex")
  simple_latex_table(
    completeness_latex,
    completeness_tex,
    caption = "Metadata completeness audit for nocturnal hypoxia/desaturation variables across cohorts included in Supplementary Table S12.",
    label = "tab:rev1_hypoxia_completeness"
  )

  limitation_path <- file.path(context$tables_dir, "Table1_Hypoxia_Limitation.txt")
  write_text_file(limitation_path, c(
    "ODI/T90/mean SpO2/nadir SpO2/AHI extraction was limited by public GEO metadata availability.",
    "GSE75097 contributed sample-level AHI only.",
    "GSE38792 contributed group-level RDI and minimum oxygen saturation from the series summary text.",
    "GSE135917, GSE89632, and GSE135251 did not expose usable desaturation fields in the local metadata files and were left as NA."
  ))

  log_message(context, step, sprintf("Hypoxia table written: %s", out_csv))
  log_message(context, step, sprintf("Hypoxia completeness audit written: %s", completeness_csv))
}
