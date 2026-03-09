#!/usr/bin/env Rscript

script_path <- tryCatch({
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(arg) > 0L) {
    normalizePath(sub("^--file=", "", arg[1]), winslash = "/", mustWork = TRUE)
  } else {
    normalizePath("rev1_major_revision/run_revision_pack.R", winslash = "/", mustWork = FALSE)
  }
}, error = function(e) normalizePath("rev1_major_revision/run_revision_pack.R", winslash = "/", mustWork = FALSE))

revision_root <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(dirname(revision_root), winslash = "/", mustWork = TRUE)
setwd(revision_root)

outputs_dir <- file.path(revision_root, "outputs")
figures_main_dir <- file.path(outputs_dir, "figures_main")
figures_supp_dir <- file.path(outputs_dir, "figures_supp")
tables_dir <- file.path(outputs_dir, "tables")
logs_dir <- file.path(outputs_dir, "logs")
session_dir <- file.path(outputs_dir, "sessionInfo")

for (path in c(outputs_dir, figures_main_dir, figures_supp_dir, tables_dir, logs_dir, session_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

REVISION_CONTEXT <- list(
  project_root = project_root,
  revision_root = revision_root,
  outputs_dir = outputs_dir,
  figures_main_dir = figures_main_dir,
  figures_supp_dir = figures_supp_dir,
  tables_dir = tables_dir,
  logs_dir = logs_dir,
  session_dir = session_dir,
  config_path = file.path(revision_root, "config.yml"),
  latex_figures_dir = file.path(project_root, "Frontiers_LaTeX_Templates", "figures"),
  master_log = file.path(logs_dir, sprintf("revision_pack_run_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
)
assign("REVISION_CONTEXT", REVISION_CONTEXT, envir = .GlobalEnv)

source(file.path(revision_root, "scripts", "00_revision_utils.R"), local = .GlobalEnv)

required_pkgs <- c("yaml", "ggplot2", "dplyr", "pROC", "PRROC", "glmnet", "caret", "MASS", "fastshap", "zip")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop(sprintf("Missing required packages: %s", paste(missing_pkgs, collapse = ", ")))
}

log_message(REVISION_CONTEXT, "run_revision_pack", sprintf("Project root: %s", project_root))
log_message(REVISION_CONTEXT, "run_revision_pack", sprintf("Revision root: %s", revision_root))

steps <- c(
  "01_collect_data.R",
  "02_common_genes_report.R",
  "03_ml_nested_pipeline.R",
  "04_feature_importance_and_shap.R",
  "05_calibration_and_prc.R",
  "06_cutoff_and_dca.R",
  "07_mr_sensitivity.R",
  "08_table1_desaturation.R",
  "09_export_for_zenodo.R"
)

results <- data.frame(step = steps, status = "pending", message = NA_character_, stringsAsFactors = FALSE)

for (i in seq_along(steps)) {
  step_file <- steps[i]
  step_id <- sub("\\.R$", "", step_file)
  step_path <- file.path(revision_root, "scripts", step_file)
  log_message(REVISION_CONTEXT, step_id, "Starting step")
  step_env <- new.env(parent = .GlobalEnv)
  status <- tryCatch({
    sys.source(step_path, envir = step_env)
    if (!exists("run_step", envir = step_env, inherits = FALSE)) {
      stop("run_step() is not defined")
    }
    step_env$run_step(REVISION_CONTEXT)
    write_session_info(REVISION_CONTEXT, step_id)
    list(ok = TRUE, message = "completed")
  }, error = function(e) {
    write_session_info(REVISION_CONTEXT, step_id)
    list(ok = FALSE, message = conditionMessage(e))
  })
  results$status[i] <- if (status$ok) "completed" else "failed"
  results$message[i] <- status$message
  log_message(
    REVISION_CONTEXT,
    step_id,
    if (status$ok) "Completed step" else paste("Failed step:", status$message),
    level = if (status$ok) "INFO" else "ERROR"
  )
}

results_path <- file.path(tables_dir, "revision_pack_step_status.csv")
write.csv(results, results_path, row.names = FALSE)

completed <- results$step[results$status == "completed"]
failed <- results$step[results$status == "failed"]

cat("\nRevision pack completed.\n")
cat("Completed steps:\n")
for (step_name in completed) {
  cat(sprintf("  - %s\n", step_name))
}
if (length(failed) > 0L) {
  cat("Failed steps:\n")
  for (step_name in failed) {
    cat(sprintf("  - %s\n", step_name))
  }
}
cat(sprintf("Status table: %s\n", results_path))
