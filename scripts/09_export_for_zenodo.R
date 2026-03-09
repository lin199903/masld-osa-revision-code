run_step <- function(context) {
  source(file.path(context$revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)
  config <- load_revision_config(context)
  step <- "09_export_for_zenodo"
  log_message(context, step, "Creating Zenodo release package and LaTeX mappings")

  locked_status <- ensure_locked_main_fig4_assets(context, config, overwrite = TRUE, strict = TRUE)
  log_message(
    context,
    step,
    sprintf(
      "Locked manuscript Figure 4 assets refreshed in outputs/figures_main: copied=%d missing=%d",
      length(locked_status$copied),
      length(locked_status$missing)
    )
  )

  zenodo_dir <- ensure_dir(file.path(context$project_root, "zenodo_release"))
  latex_table_dir <- ensure_dir(file.path(context$project_root, "Frontiers_LaTeX_Templates", "rev1_tables"))
  mapping_rows <- list()

  for (alias_name in names(config$latex$figure_aliases)) {
    source_rel <- config$latex$figure_aliases[[alias_name]]
    source_path <- resolve_revision_path(context, source_rel)
    if (!file.exists(source_path)) {
      mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
        asset_name = alias_name,
        source_relative = source_rel,
        latex_target = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }
    target_path <- sync_latex_alias(context, source_rel, alias_name)
    mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
      asset_name = alias_name,
      source_relative = source_rel,
      latex_target = normalizePath(target_path, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  }

  for (asset_name in names(config$latex$supplementary_assets)) {
    source_rel <- config$latex$supplementary_assets[[asset_name]]
    source_path <- resolve_revision_path(context, source_rel)
    if (grepl("\\.pdf$", asset_name, ignore.case = TRUE)) {
      target_path <- file.path(context$latex_figures_dir, asset_name)
    } else {
      target_path <- file.path(latex_table_dir, asset_name)
    }
    if (file.exists(source_path)) {
      file.copy(source_path, target_path, overwrite = TRUE)
      mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
        asset_name = asset_name,
        source_relative = source_rel,
        latex_target = normalizePath(target_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    } else {
      mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
        asset_name = asset_name,
        source_relative = source_rel,
        latex_target = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  mapping_df <- do.call(rbind, mapping_rows)
  mapping_path <- file.path(context$tables_dir, "LaTeX_Revision1_Mapping.csv")
  write.csv(mapping_df, mapping_path, row.names = FALSE)

  submission_pack_dir <- file.path(context$revision_root, "submission_pack")
  if (dir.exists(submission_pack_dir)) {
    submission_pack_fig_dir <- ensure_dir(file.path(submission_pack_dir, "figures"))
    synced_pack <- 0L
    for (alias_name in names(config$latex$figure_aliases)) {
      source_rel <- config$latex$figure_aliases[[alias_name]]
      source_path <- resolve_revision_path(context, source_rel)
      if (!file.exists(source_path)) {
        next
      }
      ok_root <- file.copy(source_path, file.path(submission_pack_dir, alias_name), overwrite = TRUE)
      ok_fig <- file.copy(source_path, file.path(submission_pack_fig_dir, alias_name), overwrite = TRUE)
      if (ok_root || ok_fig) {
        synced_pack <- synced_pack + 1L
      }
    }
    log_message(context, step, sprintf("Synced figure aliases to submission_pack: %d file(s)", synced_pack))
  }

  readme_lines <- c(
    "# MASLD+OSA Reviewer 1 Major Revision Evidence Pack",
    "",
    "This Zenodo-ready bundle contains the reproducible code and derived outputs generated under rev1_major_revision/.",
    "",
    "## Run command",
    "",
    "```bash",
    "Rscript rev1_major_revision/run_revision_pack.R",
    "```",
    "",
    "## Key inputs",
    "- GSE89632 liver expression matrix and phenotype table",
    "- GSE135251 liver external validation cohort",
    "- Discovery-stage cross-disease 29-gene candidate list",
    "- Existing MR sensitivity outputs under step5/MR_Results_2026-02-02_11-45",
    "",
    "## Negative-control fallback",
    "If automatic GEO download fails, place manual files under rev1_major_revision/manual_data/negative_controls/<GSE>/ and rerun the pipeline.",
    "See rev1_major_revision/outputs/tables/negative_control_manual_download_instructions.txt for the expected layout.",
    "",
    "## Derived outputs included here",
    "- code_clean.zip",
    "- tables_machine_readable.zip",
    "- figures_main.zip",
    "- figures_supp.zip",
    "- environment.yml",
    "- LICENSE",
    "",
    "## Notes",
    "- All generated paths are relative to the project root.",
    "- Session info for every major step is stored under rev1_major_revision/outputs/sessionInfo/.",
    "- Discovery-era GSE126848/GSE130970 holdout analyses are sensitivity analyses because the local matrix was pre-harmonized in the legacy workflow.",
    "- FIB-4/NFS comparisons were not possible because the public transcriptomic metadata lacked the required clinical covariates."
  )
  readme_path <- file.path(zenodo_dir, "README.md")
  write_text_file(readme_path, readme_lines)

  license_lines <- c(
    "Creative Commons Attribution 4.0 International (CC BY 4.0)",
    "",
    "This Zenodo update covers the derived figures/tables and the revision scripts in this repository.",
    "Original public datasets remain subject to their original source licenses and terms of use."
  )
  license_path <- file.path(zenodo_dir, "LICENSE")
  write_text_file(license_path, license_lines)

  pkg_list <- c("yaml", "ggplot2", "dplyr", "pROC", "PRROC", "glmnet", "caret", "MASS", "fastshap", "zip")
  env_path <- file.path(zenodo_dir, "environment.yml")
  build_environment_yml(env_path, pkg_list)
  file.copy(env_path, file.path(context$session_dir, "environment.yml"), overwrite = TRUE)

  code_files <- c(
    file.path(context$revision_root, "config.yml"),
    file.path(context$revision_root, "run_revision_pack.R"),
    list.files(file.path(context$revision_root, "scripts"), full.names = TRUE),
    list.files(file.path(context$revision_root, "manual_data", "locked_main_figure4"), full.names = TRUE),
    file.path(context$revision_root, "CHANGELOG_revision1.md")
  )
  code_zip <- file.path(zenodo_dir, "code_clean.zip")
  if (file.exists(code_zip)) unlink(code_zip)
  zip::zipr(code_zip, files = code_files, root = context$project_root)

  table_files <- list.files(context$tables_dir, full.names = TRUE)
  table_zip <- file.path(zenodo_dir, "tables_machine_readable.zip")
  if (file.exists(table_zip)) unlink(table_zip)
  zip::zipr(table_zip, files = table_files, root = context$project_root)

  main_fig_files <- list.files(context$figures_main_dir, full.names = TRUE)
  main_zip <- file.path(zenodo_dir, "figures_main.zip")
  if (file.exists(main_zip)) unlink(main_zip)
  zip::zipr(main_zip, files = main_fig_files, root = context$project_root)

  supp_fig_files <- list.files(context$figures_supp_dir, full.names = TRUE)
  supp_zip <- file.path(zenodo_dir, "figures_supp.zip")
  if (file.exists(supp_zip)) unlink(supp_zip)
  zip::zipr(supp_zip, files = supp_fig_files, root = context$project_root)

  log_message(context, step, sprintf("LaTeX mapping written: %s", mapping_path))
  log_message(context, step, sprintf("Zenodo release prepared under: %s", zenodo_dir))
}
