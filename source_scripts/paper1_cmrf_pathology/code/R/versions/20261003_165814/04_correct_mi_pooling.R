#!/usr/bin/env Rscript
# Author: mjwang. Reuse the frozen 50 imputations and formulas; correct k only.
suppressPackageStartupMessages({
  library(mice)
  library(digest)
})
version_id <- "20261003_165814"
source_version <- "20260809_152651"
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(arg) != 1L) stop("Run using Rscript.")
script_path <- normalizePath(sub("^--file=", "", arg))
paper_dir <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))
source_dir <- file.path(paper_dir, "results", "versions", source_version)
out_dir <- file.path(paper_dir, "results", "versions", version_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
target_names <- c("main_models_multiple_imputation.csv", "component_models_multiple_imputation.csv",
  "glycaemia_decomposition_multiple_imputation.csv", "multiple_imputation_precision_diagnostics.csv",
  "Table2_main_models.csv", "Table3_component_models.csv", "robustness_models.csv",
  "mi_pooling_old_new_differences.csv", "mi_pooling_parameter_checks.csv",
  "mi_pooling_per_imputation_estimates.csv", "mi_pooling_input_manifest.csv",
  "mi_pooling_method_note.md", "sessionInfo_mi_pooling.txt")
if (any(file.exists(file.path(out_dir, target_names)))) stop("Refusing to overwrite correction outputs.")
old_script <- file.path(paper_dir, "code", "R", "versions", source_version, "01_reanalysis.R")
source_names <- c("multiple_imputation_object.rds", "main_models_multiple_imputation.csv",
  "component_models_multiple_imputation.csv", "glycaemia_decomposition_multiple_imputation.csv",
  "multiple_imputation_precision_diagnostics.csv", "main_models_complete_component.csv",
  "component_models_complete_component.csv", "robustness_models.csv")
source_paths <- c(file.path(source_dir, source_names), old_script)
input_manifest <- data.frame(
  relative_path = c(file.path("results", "versions", source_version, source_names),
    file.path("code", "R", "versions", source_version, "01_reanalysis.R")),
  bytes = file.info(source_paths)$size,
  sha256 = vapply(source_paths, function(p) digest(p, algo = "sha256", file = TRUE), character(1)),
  stringsAsFactors = FALSE
)
# Evaluate only the two function definitions, not the generating source pipeline.
source_expressions <- parse(old_script)
for (fn in c("derive_cmrf", "pool_scalar_rubin")) {
  matches <- vapply(source_expressions, function(e) is.call(e) &&
    identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(fn)), logical(1))
  if (sum(matches) != 1L) stop("Expected unique source function: ", fn)
  eval(source_expressions[[which(matches)]], envir = .GlobalEnv)
}
imp <- readRDS(file.path(source_dir, "multiple_imputation_object.rds"))
stopifnot(inherits(imp, "mids"), nrow(imp$data) == 1085L, imp$m == 50L, imp$iteration == 20L,
  identical(sort(levels(imp$data$center_code)), c("HZ", "JL", "RJ")))

fits <- list()
parameter_rows <- list()
per_imputation <- list()
for (i in seq_len(imp$m)) {
  d <- mice::complete(imp, i)
  for (v in c("sex", "dm", "htn", "alcohol", "hcv")) d[[v]] <- as.integer(as.character(d[[v]]))
  d$sex_label <- factor(d$sex, levels = c(0, 1), labels = c("Male", "Female"))
  d <- cbind(d, derive_cmrf(d))
  stopifnot(!anyNA(d[paste0("cmrf", 1:5)]))
  d$f_ge2 <- as.integer(d$fibrosis >= 2)
  d$f_ge3 <- as.integer(d$fibrosis >= 3)
  d$glycaemia_category <- factor(ifelse(d$dm == 1, "Recorded T2DM",
    ifelse(d$cmrf2 == 1, "Dysglycaemia without recorded T2DM", "Normoglycaemia")),
    levels = c("Normoglycaemia", "Dysglycaemia without recorded T2DM", "Recorded T2DM"))
  for (outcome in c("f_ge2", "f_ge3")) {
    rhs <- list(
      main = c("cmrf_count", "age", "sex_label", "center_code"),
      component = c(paste0("cmrf", 1:5), "age", "sex_label", "center_code"),
      glycaemia = c("glycaemia_category", "cmrf1", "cmrf3", "cmrf4", "cmrf5", "age", "sex_label", "center_code")
    )
    for (kind in names(rhs)) {
      fit <- glm(reformulate(rhs[[kind]], outcome), data = d, family = binomial())
      stopifnot(fit$converged, nobs(fit) == 1085L,
        fit$rank == nobs(fit) - fit$df.residual,
        fit$rank == c(main = 6L, component = 10L, glycaemia = 11L)[[kind]])
      parameter_rows[[length(parameter_rows)+1L]] <- data.frame(imputation = i, model = kind,
        outcome = outcome, N = nobs(fit), rank = fit$rank, df_residual = fit$df.residual,
        old_k = c(main = 7L, component = 11L, glycaemia = 12L)[[kind]],
        new_k = fit$rank, formula = paste(deparse(formula(fit)), collapse = " "),
        author = "mjwang", stringsAsFactors = FALSE)
      terms <- switch(kind, main = "cmrf_count", component = paste0("cmrf", 1:5),
        glycaemia = grep("^glycaemia_category", names(coef(fit)), value = TRUE))
      for (term in terms) {
        key <- paste(kind, outcome, term, sep = "::")
        z <- data.frame(imputation = i, model = kind, outcome = outcome, term = term,
          N = nobs(fit), k = fit$rank, df_residual = fit$df.residual,
          Q = unname(coef(fit)[term]), U = unname(vcov(fit)[term, term]),
          author = "mjwang", stringsAsFactors = FALSE)
        fits[[key]][[i]] <- z
        per_imputation[[length(per_imputation)+1L]] <- z
      }
    }
  }
}

old_files <- c(main = "main_models_multiple_imputation.csv",
  component = "component_models_multiple_imputation.csv",
  glycaemia = "glycaemia_decomposition_multiple_imputation.csv")
new_tables <- list()
diff_rows <- list()
max_pool_difference <- 0
for (kind in names(old_files)) {
  old <- read.csv(file.path(source_dir, old_files[[kind]]), check.names = FALSE,
    stringsAsFactors = FALSE)
  updated <- old
  updated$pooling_k <- NA_integer_
  updated$pooling_df_complete <- NA_integer_
  for (j in seq_len(nrow(old))) {
    key <- paste(kind, old$outcome[j], old$term[j], sep = "::")
    estimates <- do.call(rbind, fits[[key]])
    stopifnot(nrow(estimates) == 50L, length(unique(estimates$N)) == 1L,
      length(unique(estimates$k)) == 1L)
    n <- estimates$N[1]
    k <- estimates$k[1]
    old_k <- c(main = 7L, component = 11L, glycaemia = 12L)[[kind]]
    # Confirm the old exported result is reproduced using its legacy k.
    previous <- pool_scalar_rubin(estimates$Q, estimates$U, n, old_k)
    for (field in names(previous)) stopifnot(abs(unname(previous[field]) - old[[field]][j]) < 1e-10)
    correct <- pool_scalar_rubin(estimates$Q, estimates$U, n, k)
    reference <- mice::pool.scalar(estimates$Q, estimates$U, n = n, k = k)
    differences <- c(correct["estimate"] - reference$qbar,
      correct["std_error"] - sqrt(reference$t), correct["df"] - reference$df,
      correct["fmi"] - reference$fmi)
    max_pool_difference <- max(max_pool_difference, abs(differences))
    stopifnot(max(abs(differences)) < 1e-10,
      abs(correct["estimate"] - old$estimate[j]) < 1e-12,
      abs(correct["std_error"] - old$std_error[j]) < 1e-12)
    for (field in names(correct)) updated[[field]][j] <- unname(correct[field])
    updated$pooling_k[j] <- k
    updated$pooling_df_complete[j] <- n-k
    for (field in names(correct)) diff_rows[[length(diff_rows)+1L]] <- data.frame(
      model = kind, outcome = old$outcome[j], term = old$term[j], field = field,
      old_value = old[[field]][j], corrected_value = unname(correct[field]),
      corrected_minus_old = unname(correct[field]) - old[[field]][j],
      old_k = old_k, corrected_k = k, author = "mjwang", stringsAsFactors = FALSE)
  }
  if (kind == "component") updated$p_fdr_all_10 <- p.adjust(updated$p_value, "BH")
  if (kind == "glycaemia") updated$p_holm_all_4 <- p.adjust(updated$p_value, "holm")
  adjusted_field <- if (kind == "component") "p_fdr_all_10" else if (kind == "glycaemia") "p_holm_all_4" else NULL
  if (!is.null(adjusted_field)) for (j in seq_len(nrow(old))) diff_rows[[length(diff_rows)+1L]] <- data.frame(
    model = kind, outcome = old$outcome[j], term = old$term[j], field = adjusted_field,
    old_value = old[[adjusted_field]][j], corrected_value = updated[[adjusted_field]][j],
    corrected_minus_old = updated[[adjusted_field]][j] - old[[adjusted_field]][j],
    old_k = c(main = 7L, component = 11L, glycaemia = 12L)[[kind]],
    corrected_k = updated$pooling_k[j], author = "mjwang", stringsAsFactors = FALSE)
  new_tables[[kind]] <- updated
}

precision_rows <- list()
for (outcome in c("f_ge2", "f_ge3")) {
  z <- do.call(rbind, fits[[paste("main", outcome, "cmrf_count", sep = "::")]])
  pool <- mice::pool.scalar(z$Q, z$U, n = z$N[1], k = z$k[1])
  lambda <- (1+1/pool$m)*pool$b/pool$t
  precision_rows[[length(precision_rows)+1L]] <- data.frame(outcome = outcome, m = pool$m,
    pooled_log_or = pool$qbar, pooled_se = sqrt(pool$t), within_variance = pool$ubar,
    between_variance = pool$b, fraction_missing_information_approx = lambda,
    monte_carlo_error_log_or = sqrt(pool$b/pool$m),
    monte_carlo_error_pct_pooled_se = 100*sqrt(pool$b/pool$m)/sqrt(pool$t),
    pooling_k = z$k[1], pooling_df_complete = z$df_residual[1],
    pooled_df = pool$df, pooled_fmi = pool$fmi, stringsAsFactors = FALSE)
}
precision <- do.call(rbind, precision_rows)
old_precision <- read.csv(file.path(source_dir, "multiple_imputation_precision_diagnostics.csv"))
for (field in names(old_precision)[-1]) stopifnot(max(abs(old_precision[[field]] - precision[[field]])) < 1e-12)

cc <- read.csv(file.path(source_dir, "main_models_complete_component.csv"), stringsAsFactors = FALSE)
cc <- transform(cc, m = NA_integer_, df = NA_real_, fmi = NA_real_)
table2 <- rbind(cc, new_tables$main[names(cc)])
cc_component <- read.csv(file.path(source_dir, "component_models_complete_component.csv"), stringsAsFactors = FALSE)
cc_component$dataset <- "Complete-component"
mi_component <- new_tables$component
mi_component$dataset <- "Multiple imputation"
table3 <- rbind(cc_component, mi_component[names(cc_component)])
robustness <- read.csv(file.path(source_dir, "robustness_models.csv"), stringsAsFactors = FALSE)
for (outcome in c("f_ge2", "f_ge3")) {
  idx <- which(robustness$analysis == "Adult raw-variable MI" & robustness$outcome == outcome)
  z <- new_tables$main[new_tables$main$outcome == outcome, , drop = FALSE]
  for (field in intersect(names(robustness), names(z))) {
    if (field != "analysis") robustness[[field]][idx] <- z[[field]][1]
  }
}
differences <- do.call(rbind, diff_rows)
stopifnot(all(vapply(source_paths, function(p) digest(p, algo = "sha256", file = TRUE), character(1)) == input_manifest$sha256))
for (kind in names(old_files)) write.csv(new_tables[[kind]], file.path(out_dir, old_files[[kind]]), row.names = FALSE)
write.csv(precision, file.path(out_dir, "multiple_imputation_precision_diagnostics.csv"), row.names = FALSE)
write.csv(table2, file.path(out_dir, "Table2_main_models.csv"), row.names = FALSE)
write.csv(table3, file.path(out_dir, "Table3_component_models.csv"), row.names = FALSE)
write.csv(robustness, file.path(out_dir, "robustness_models.csv"), row.names = FALSE)
write.csv(differences, file.path(out_dir, "mi_pooling_old_new_differences.csv"), row.names = FALSE)
write.csv(do.call(rbind, parameter_rows), file.path(out_dir, "mi_pooling_parameter_checks.csv"), row.names = FALSE)
write.csv(do.call(rbind, per_imputation), file.path(out_dir, "mi_pooling_per_imputation_estimates.csv"), row.names = FALSE)
write.csv(input_manifest, file.path(out_dir, "mi_pooling_input_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_mi_pooling.txt"))
summary_delta <- aggregate(abs(corrected_minus_old) ~ field, data = differences, FUN = max)
note <- c("---", "author: mjwang", "revision_author: mjwang", paste0("version: ", version_id),
  paste0("source_analysis_version: ", source_version), "---", "",
  "# Restricted three-centre MI pooling correction", "",
  "The original 50 completed datasets (20 iterations) were reused from the frozen mids object. No new imputation or stochastic draw was performed. The original pipeline set the RNG seed to 2026073012 before imputation; the stored mids seed slot is NA because the original call used set.seed externally.", "",
  "The three model formulas, 1,085 participants, component re-derivation, and outcomes were kept identical. The original derive_cmrf and Rubin pooling function definitions were loaded in isolation from the frozen script. The old results were reproduced with k=7/11/12. Each new fit recorded rank and residual degrees of freedom, establishing k=6/10/11 including the intercept for the burden/component/glycaemia-decomposition models.", "",
  "Corrected pooling used n-k from those actual model ranks. Pooled log estimates, SEs, df and FMI matched mice::pool.scalar(n=1085,k=actual rank) within 1e-10. All pooled point estimates and SEs matched the old frozen files within 1e-12. The t-based CI, P values and adjusted P values were updated. Original files were hash-checked before and after and were unchanged.", "",
  paste0("Maximum discrepancy from mice::pool.scalar: ", format(max_pool_difference, digits = 17)), "",
  "Maximum absolute corrected-minus-old difference by field:", "",
  capture.output(print(summary_delta, row.names = FALSE, digits = 12)), "",
  "Monte Carlo variance, SE and simulation error were unchanged; finite-df pooled FMI/df are additionally provided in the corrected precision file. Source hashes, per-imputation coefficients/variances and model ranks are retained. Main-table/source-data combined outputs preserve complete-component results while replacing only MI estimates.", "",
  "This correction resolves a pooling implementation residue; it does not establish MAR, clinical validity, prospective prediction, or external validation. No changes were made to preprocessing, model selection, thresholds, sample eligibility, or centre mapping.")
writeLines(note, file.path(out_dir, "mi_pooling_method_note.md"))
print(summary_delta, row.names = FALSE, digits = 12)
message("Completed restricted MI correction: ", out_dir)
