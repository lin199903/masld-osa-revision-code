options(stringsAsFactors = FALSE)

ensure_revision_context <- function(context = NULL) {
  if (!is.null(context)) {
    return(context)
  }
  if (exists("REVISION_CONTEXT", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("REVISION_CONTEXT", envir = .GlobalEnv, inherits = FALSE))
  }
  stop("REVISION_CONTEXT is not available.")
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

resolve_revision_path <- function(context, relative_path) {
  normalizePath(file.path(context$revision_root, relative_path), winslash = "/", mustWork = FALSE)
}

resolve_project_path <- function(context, relative_path) {
  normalizePath(file.path(context$project_root, relative_path), winslash = "/", mustWork = FALSE)
}

load_revision_config <- function(context = NULL) {
  context <- ensure_revision_context(context)
  yaml::read_yaml(context$config_path)
}

step_log_path <- function(context, step) {
  file.path(context$logs_dir, paste0(step, ".log"))
}

log_message <- function(context, step, message, level = "INFO") {
  context <- ensure_revision_context(context)
  line <- sprintf("[%s] [%s] [%s] %s",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  level,
                  step,
                  message)
  cat(line, "\n")
  cat(line, "\n", file = context$master_log, append = TRUE)
  cat(line, "\n", file = step_log_path(context, step), append = TRUE)
}

write_text_file <- function(path, lines) {
  ensure_dir(dirname(path))
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

write_session_info <- function(context, step) {
  context <- ensure_revision_context(context)
  path <- file.path(context$session_dir, paste0(step, "_sessionInfo.txt"))
  capture.output(sessionInfo(), file = path)
  invisible(path)
}

clip_probs <- function(prob, eps = 1e-6) {
  pmin(pmax(as.numeric(prob), eps), 1 - eps)
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L || all(!is.finite(x))) {
    return(NA_real_)
  }
  mean(x[is.finite(x)])
}

read_expression_matrix <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  gene_col <- if (nzchar(colnames(df)[1])) colnames(df)[1] else 1
  genes <- toupper(trimws(as.character(df[[gene_col]])))
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  if (anyDuplicated(rownames(mat))) {
    counts <- table(rownames(mat))
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
    mat <- mat / as.numeric(counts[rownames(mat)])
  }
  mat
}

read_table_file <- function(path) {
  read.csv(path, check.names = FALSE)
}

align_expression_with_pheno <- function(expr, pheno, sample_col, group_col, positive_label) {
  sample_ids <- as.character(pheno[[sample_col]])
  keep <- sample_ids %in% colnames(expr)
  pheno <- pheno[keep, , drop = FALSE]
  sample_ids <- as.character(pheno[[sample_col]])
  expr <- expr[, sample_ids, drop = FALSE]
  y <- as.integer(as.character(pheno[[group_col]]) == positive_label)
  rownames(pheno) <- sample_ids
  list(expr = expr, pheno = pheno, y = y, sample_ids = sample_ids)
}

extract_feature_matrix <- function(expr, genes, sample_ids = NULL) {
  genes <- toupper(genes)
  if (!is.null(sample_ids)) {
    expr <- expr[, sample_ids, drop = FALSE]
  }
  missing <- setdiff(genes, rownames(expr))
  if (length(missing) > 0L) {
    stop(sprintf("Missing genes: %s", paste(missing, collapse = ", ")))
  }
  x <- t(expr[genes, , drop = FALSE])
  colnames(x) <- genes
  x
}

compute_ihs <- function(expr, genes, sample_ids = NULL) {
  x <- extract_feature_matrix(expr, genes, sample_ids)
  rowMeans(x)
}

roc_auc_ci <- function(y, prob) {
  if (length(unique(y)) < 2L) {
    return(c(auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  roc_obj <- pROC::roc(response = y, predictor = prob, quiet = TRUE, direction = "<")
  ci <- as.numeric(pROC::ci.auc(roc_obj))
  c(auc = as.numeric(pROC::auc(roc_obj)), ci_low = ci[1], ci_high = ci[3])
}

pr_auc_case <- function(y, prob) {
  pos <- prob[y == 1]
  neg <- prob[y == 0]
  if (length(pos) == 0L || length(neg) == 0L) {
    return(NA_real_)
  }
  PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg, curve = TRUE)$auc.integral
}

pr_curve_df <- function(y, prob, control_positive = FALSE) {
  if (control_positive) {
    y_use <- 1L - y
    prob_use <- 1 - prob
  } else {
    y_use <- y
    prob_use <- prob
  }
  pos <- prob_use[y_use == 1]
  neg <- prob_use[y_use == 0]
  if (length(pos) == 0L || length(neg) == 0L) {
    return(data.frame())
  }
  curve <- PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg, curve = TRUE)
  data.frame(recall = curve$curve[, 1], precision = curve$curve[, 2])
}

calibration_metrics <- function(y, prob) {
  prob <- clip_probs(prob)
  lp <- qlogis(prob)
  intercept <- tryCatch({
    suppressWarnings(coef(glm(y ~ 1, family = binomial(), offset = lp))[1])
  }, error = function(e) NA_real_)
  slope <- tryCatch({
    suppressWarnings(coef(glm(y ~ lp, family = binomial()))[2])
  }, error = function(e) NA_real_)
  c(
    intercept = as.numeric(intercept),
    slope = as.numeric(slope),
    brier = mean((prob - y) ^ 2)
  )
}

threshold_metrics <- function(y, prob, threshold) {
  if (!is.finite(threshold)) {
    return(data.frame(
      threshold = threshold,
      tp = NA_real_,
      fp = NA_real_,
      tn = NA_real_,
      fn = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_,
      ppv = NA_real_,
      npv = NA_real_,
      f1 = NA_real_,
      accuracy = NA_real_
    ))
  }
  pred <- as.integer(prob >= threshold)
  tp <- sum(pred == 1 & y == 1, na.rm = TRUE)
  fp <- sum(pred == 1 & y == 0, na.rm = TRUE)
  tn <- sum(pred == 0 & y == 0, na.rm = TRUE)
  fn <- sum(pred == 0 & y == 1, na.rm = TRUE)
  sensitivity <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  ppv <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  npv <- if ((tn + fn) == 0) NA_real_ else tn / (tn + fn)
  f1 <- if ((2 * tp + fp + fn) == 0) NA_real_ else (2 * tp) / (2 * tp + fp + fn)
  accuracy <- (tp + tn) / length(y)
  data.frame(
    threshold = threshold,
    tp = tp,
    fp = fp,
    tn = tn,
    fn = fn,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    f1 = f1,
    accuracy = accuracy
  )
}

threshold_with_ci <- function(y, prob, threshold, label) {
  metrics <- threshold_metrics(y, prob, threshold)
  sens_ci <- if ((metrics$tp + metrics$fn) > 0) binom.test(metrics$tp, metrics$tp + metrics$fn)$conf.int else c(NA_real_, NA_real_)
  spec_ci <- if ((metrics$tn + metrics$fp) > 0) binom.test(metrics$tn, metrics$tn + metrics$fp)$conf.int else c(NA_real_, NA_real_)
  cbind(
    data.frame(label = label),
    metrics,
    sensitivity_ci_low = sens_ci[1],
    sensitivity_ci_high = sens_ci[2],
    specificity_ci_low = spec_ci[1],
    specificity_ci_high = spec_ci[2]
  )
}

youden_threshold <- function(y, prob) {
  roc_obj <- pROC::roc(response = y, predictor = prob, quiet = TRUE, direction = "<")
  coords <- pROC::coords(roc_obj, x = "best", best.method = "youden", transpose = FALSE)
  as.numeric(coords["threshold"])
}

find_threshold_for_sensitivity <- function(y, prob, target = 0.90) {
  thresholds <- sort(unique(prob), decreasing = TRUE)
  metrics <- do.call(rbind, lapply(thresholds, function(th) threshold_metrics(y, prob, th)))
  eligible <- metrics[metrics$sensitivity >= target, , drop = FALSE]
  if (nrow(eligible) == 0L) {
    return(NA_real_)
  }
  eligible <- eligible[order(-eligible$specificity, eligible$threshold), , drop = FALSE]
  eligible$threshold[1]
}

calibration_curve_df <- function(y, prob, n_bins = 10L) {
  prob <- clip_probs(prob)
  n_use <- max(2L, min(n_bins, length(prob)))
  bins <- cut(rank(prob, ties.method = "first"), breaks = n_use, include.lowest = TRUE, labels = FALSE)
  out <- aggregate(data.frame(prob = prob, obs = y), by = list(bin = bins), FUN = mean)
  counts <- aggregate(y, by = list(bin = bins), FUN = length)
  out$n <- counts$x[match(out$bin, counts$bin)]
  out
}

roc_curve_df <- function(y, prob) {
  roc_obj <- pROC::roc(response = y, predictor = prob, quiet = TRUE, direction = "<")
  data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities,
    false_positive_rate = 1 - roc_obj$specificities
  )
}

make_folds <- function(y, k, seed) {
  set.seed(seed)
  caret::createFolds(factor(y), k = k, returnTrain = FALSE)
}

fit_inner_glmnet <- function(x, y, alpha_grid, inner_folds, seed) {
  best <- NULL
  k_use <- max(2L, min(inner_folds, min(sum(y == 1), sum(y == 0))))
  for (alpha in alpha_grid) {
    set.seed(seed + as.integer(alpha * 100))
    cv_fit <- glmnet::cv.glmnet(
      x = x,
      y = y,
      family = "binomial",
      alpha = alpha,
      type.measure = "auc",
      nfolds = k_use,
      standardize = TRUE
    )
    score <- max(cv_fit$cvm, na.rm = TRUE)
    if (is.null(best) || score > best$inner_auc) {
      best <- list(alpha = alpha, cv_fit = cv_fit, inner_auc = score)
    }
  }
  lambda_name <- if (identical(best$cv_fit$lambda.1se, NA_real_)) "lambda.min" else "lambda.1se"
  lambda_value <- if (lambda_name == "lambda.1se") best$cv_fit$lambda.1se else best$cv_fit$lambda.min
  final_fit <- glmnet::glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = best$alpha,
    lambda = lambda_value,
    standardize = TRUE
  )
  beta <- as.matrix(stats::coef(final_fit))
  nonzero <- rownames(beta)[abs(beta[, 1]) > 0]
  nonzero <- setdiff(nonzero, "(Intercept)")
  list(
    alpha = best$alpha,
    lambda = lambda_value,
    inner_auc = best$inner_auc,
    fit = final_fit,
    coefficients = beta,
    selected_features = nonzero
  )
}

predict_glmnet_prob <- function(model, x) {
  as.numeric(glmnet::predict.glmnet(model, newx = x, type = "response"))
}

decision_curve_df <- function(y, prob, model_name, thresholds = seq(0.05, 0.95, by = 0.05)) {
  n <- length(y)
  prevalence <- mean(y)
  model_nb <- vapply(thresholds, function(th) {
    pred <- prob >= th
    tp <- sum(pred & y == 1) / n
    fp <- sum(pred & y == 0) / n
    tp - fp * (th / (1 - th))
  }, numeric(1))
  data.frame(
    threshold = thresholds,
    model = model_name,
    net_benefit = model_nb,
    prevalence = prevalence
  )
}

baseline_decision_curves <- function(y, thresholds = seq(0.05, 0.95, by = 0.05)) {
  prevalence <- mean(y)
  rbind(
    data.frame(threshold = thresholds, model = "Treat none", net_benefit = 0, prevalence = prevalence),
    data.frame(
      threshold = thresholds,
      model = "Treat all",
      net_benefit = prevalence - (1 - prevalence) * (thresholds / (1 - thresholds)),
      prevalence = prevalence
    )
  )
}

simple_latex_table <- function(df, path, caption = NULL, label = NULL) {
  latex_escape <- function(x) {
    x <- as.character(x)
    x <- gsub("\\\\", "\\\\textbackslash{}", x, perl = TRUE)
    x <- gsub("([%&_#\\$\\{\\}])", "\\\\\\1", x, perl = TRUE)
    x
  }
  escaped_df <- as.data.frame(lapply(df, latex_escape), stringsAsFactors = FALSE, check.names = FALSE)
  header <- paste(latex_escape(colnames(escaped_df)), collapse = " & ")
  rows <- apply(escaped_df, 1, function(x) paste(x, collapse = " & "))
  lines <- c("\\begin{table}[ht]", "\\centering")
  if (!is.null(caption)) {
    lines <- c(lines, paste0("\\caption{", caption, "}"))
  }
  if (!is.null(label)) {
    lines <- c(lines, paste0("\\label{", label, "}"))
  }
  lines <- c(
    lines,
    paste0("\\begin{tabular}{", paste(rep("l", ncol(df)), collapse = ""), "}"),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline",
    paste0(rows, " \\\\"),
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )
  write_text_file(path, lines)
  invisible(path)
}

sync_latex_alias <- function(context, source_relative, alias_name) {
  source_path <- resolve_revision_path(context, source_relative)
  target_path <- file.path(context$latex_figures_dir, alias_name)
  if (!file.exists(source_path)) {
    stop(sprintf("Cannot sync missing asset: %s", source_path))
  }
  ensure_dir(dirname(target_path))
  ok <- file.copy(source_path, target_path, overwrite = TRUE)
  if (!ok) {
    stop(sprintf("Failed to copy %s to %s", source_path, target_path))
  }
  invisible(target_path)
}

ensure_locked_main_fig4_assets <- function(context, config, overwrite = TRUE, strict = TRUE) {
  locked_cfg <- config$manuscript_locked_fig4
  if (is.null(locked_cfg) || isFALSE(locked_cfg$enabled)) {
    return(list(copied = character(0), missing = character(0), skipped = TRUE))
  }
  if (is.null(locked_cfg$source_dir) || !nzchar(locked_cfg$source_dir)) {
    stop("config.yml: manuscript_locked_fig4.source_dir is missing.")
  }
  if (is.null(locked_cfg$files) || length(locked_cfg$files) == 0L) {
    stop("config.yml: manuscript_locked_fig4.files is empty.")
  }

  source_dir <- resolve_revision_path(context, locked_cfg$source_dir)
  ensure_dir(context$figures_main_dir)

  copied <- character(0)
  missing <- character(0)
  for (asset_name in as.character(unlist(locked_cfg$files))) {
    source_path <- normalizePath(file.path(source_dir, asset_name), winslash = "/", mustWork = FALSE)
    target_path <- normalizePath(file.path(context$figures_main_dir, asset_name), winslash = "/", mustWork = FALSE)
    if (!file.exists(source_path)) {
      missing <- c(missing, source_path)
      next
    }
    if (file.exists(target_path) && !isTRUE(overwrite)) {
      copied <- c(copied, target_path)
      next
    }
    ok <- file.copy(source_path, target_path, overwrite = overwrite)
    if (!ok) {
      stop(sprintf("Failed to copy locked Figure 4 asset: %s -> %s", source_path, target_path))
    }
    copied <- c(copied, target_path)
  }

  if (length(missing) > 0L && isTRUE(strict)) {
    stop(sprintf(
      "Locked Figure 4 assets missing (%d): %s",
      length(missing),
      paste(missing, collapse = "; ")
    ))
  }

  list(copied = copied, missing = missing, skipped = FALSE)
}

build_environment_yml <- function(path, packages) {
  lines <- c(
    "name: masld-osa-revision1",
    "channels:",
    "  - conda-forge",
    "dependencies:",
    sprintf("  - r-base=%s", paste(R.version$major, sub(" .*", "", R.version$minor), sep = ".")),
    "  - r-essentials",
    sprintf("  - r-%s", sort(unique(packages)))
  )
  write_text_file(path, lines)
}
