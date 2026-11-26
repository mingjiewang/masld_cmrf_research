#!/usr/bin/env Rscript

# Author: mjwang
# Derive two steatosis probability panels from the frozen complete-component
# cohort. Preserve the original three-threshold model estimates and BH3 values.
# Usage: Rscript --vanilla 01_steatosis_probabilities.R --project-root PATH
# Optional --output-version YYYYMMDD_HHMMSS writes a fresh reproducibility run.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(key, default = NULL) {
  at <- which(args == key)
  if (!length(at)) return(default)
  if (length(at) != 1L || at == length(args)) stop("Invalid argument: ", key)
  args[at + 1L]
}
project_arg <- arg_value("--project-root")
if (is.null(project_arg)) stop("Supply --project-root PATH.")
project_dir <- normalizePath(project_arg, mustWork = TRUE)
version_id <- arg_value("--output-version", "20261027_230711")
if (!grepl("^[0-9]{8}_[0-9]{6}$", version_id)) stop("Invalid output version.")
if (!requireNamespace("digest", quietly = TRUE) ||
    !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("The existing digest and jsonlite packages are required.")
}
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
paper_dir <- file.path(project_dir, "paper1_cmrf_pathology")
frozen_dir <- file.path(paper_dir, "results", "versions", "20260809_152651")
input_csv <- file.path(frozen_dir, "canonical_adult_complete_component.csv")
original_csv <- file.path(frozen_dir, "steatosis_exploratory.csv")
reference_code <- file.path(paper_dir, "code", "R", "versions",
                            "20260809_152651", "01_reanalysis.R")
out_dir <- file.path(paper_dir, "results", "versions", version_id)
source_dir <- file.path(out_dir, "source_data")
qa_dir <- file.path(out_dir, "qa")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
outputs <- c(
  file.path(source_dir, c("standardized_steatosis_probabilities.csv",
    "steatosis_distribution_by_cmrf.csv", "observed_steatosis_wilson.csv",
    "steatosis_secondary_models.csv", "original_model_reproduction.csv")),
  file.path(qa_dir, c("steatosis_bootstrap_diagnostics.csv",
    "steatosis_bootstrap_probabilities.csv", "steatosis_analysis_checks.csv",
    "steatosis_analysis_manifest.json", "steatosis_sessionInfo.txt"))
)
if (any(file.exists(outputs))) {
  stop("Analysis outputs already exist; use a new --output-version.")
}
hash_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
input_paths <- c(canonical = input_csv, original_models = original_csv,
                 reference_method = reference_code)
input_hashes <- vapply(input_paths, hash_file, character(1))
started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
checks <- list()
check <- function(label, pass, detail = "") {
  ok <- isTRUE(pass)
  checks[[length(checks) + 1L]] <<- data.frame(
    check = label, passed = ok, detail = as.character(detail),
    stringsAsFactors = FALSE
  )
  if (!ok) stop("Check failed: ", label, ": ", detail)
}
canonical <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
original <- read.csv(original_csv, stringsAsFactors = FALSE, check.names = FALSE)
check("canonical_n_935", nrow(canonical) == 935L)
check("unique_participants", !anyDuplicated(canonical$id))
check("adult_complete_component", all(canonical$age >= 18) &&
  all(canonical$complete_component) &&
  all(canonical$cmrf_count == rowSums(canonical[paste0("cmrf", 1:5)])))
check("source_prefix_matches_id", all(canonical$source_prefix == substr(canonical$id, 1, 2)))
derived_center <- ifelse(canonical$source_prefix == "ZS", "RJ", canonical$source_prefix)
check("recorded_centres_match_ZS_to_RJ", all(canonical$center_code == derived_center))
canonical$center_code <- factor(derived_center, levels = c("HZ", "JL", "RJ"))
check("three_centres_only", !anyNA(canonical$center_code) &&
  identical(as.integer(table(canonical$center_code)), c(440L, 139L, 356L)))
check("sex_label_matches_recorded_code", all(canonical$sex_label ==
  ifelse(canonical$sex == 0, "Male", "Female")))
canonical$sex_label <- factor(canonical$sex_label, levels = c("Male", "Female"))
check("model_variables_complete", !anyNA(canonical[c(
  "cmrf_count", "age", "sex_label", "center_code", "steatosis", "s_ge2", "s_ge3")]))
check("count_groups_1_to_5", identical(sort(unique(canonical$cmrf_count)), 1:5))
check("steatosis_grades_0_to_3", all(canonical$steatosis %in% 0:3))
check("steatosis_grade_counts", identical(as.integer(table(factor(
  canonical$steatosis, levels = 0:3))), c(19L, 271L, 344L, 301L)))
for (k in 1:3) check(paste0("recorded_s_ge", k, "_matches_grade"),
  all(canonical[[paste0("s_ge", k)]] == as.integer(canonical$steatosis >= k)))
check("all_original_thresholds_retained", identical(original$outcome,
  c("s_ge1", "s_ge2", "s_ge3")) && all(original$N == 935))
check("original_BH3_matches_three_tests", max(abs(original$p_fdr_3 -
  p.adjust(original$p_value, method = "BH"))) < 1e-14)

capture_call <- function(expr) {
  messages <- character()
  error <- ""
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    error <<- conditionMessage(e)
    NULL
  })
  list(value = value, warnings = messages, error = error)
}
formula_for <- function(outcome) reformulate(
  c("cmrf_count", "age", "sex_label", "center_code"), outcome)
gcomp <- function(fit, d, values = 1:5) {
  vapply(values, function(v) {
    nd <- d
    nd$cmrf_count <- v
    mean(predict(fit, newdata = nd, type = "response"))
  }, numeric(1))
}

# Reproduce the saved associations before deriving any new probabilities.
fits <- list()
reproduction <- list()
point_diagnostics <- list()
for (outcome in c("s_ge2", "s_ge3")) {
  captured <- capture_call(glm(formula_for(outcome), data = canonical, family = binomial()))
  check(paste0(outcome, "_fit_no_error"), identical(captured$error, ""), captured$error)
  fit <- captured$value
  check(paste0(outcome, "_fit_converged"), fit$converged)
  check(paste0(outcome, "_fit_no_warning"), length(captured$warnings) == 0L,
        paste(captured$warnings, collapse = " | "))
  fits[[outcome]] <- fit
  z <- coef(summary(fit))["cmrf_count", ]
  actual <- c(N = nobs(fit), events = sum(model.response(model.frame(fit))),
    estimate = unname(z["Estimate"]), std_error = unname(z["Std. Error"]),
    OR = exp(unname(z["Estimate"])),
    ci_low = exp(unname(z["Estimate"]) - 1.96 * unname(z["Std. Error"])),
    ci_high = exp(unname(z["Estimate"]) + 1.96 * unname(z["Std. Error"])),
    p_value = unname(z["Pr(>|z|)"]))
  expected <- unlist(original[original$outcome == outcome, names(actual)], use.names = TRUE)
  tolerance <- pmax(1e-13, abs(expected) * 1e-10)
  reproduced <- data.frame(outcome = outcome, metric = names(actual),
    original = as.numeric(expected), reproduced = as.numeric(actual),
    absolute_difference = abs(actual - expected), tolerance = tolerance,
    passed = abs(actual - expected) <= tolerance, row.names = NULL)
  reproduction[[outcome]] <- reproduced
  check(paste0(outcome, "_original_model_reproduced"), all(reproduced$passed),
        paste("maximum absolute difference", format(max(reproduced$absolute_difference), digits = 16)))
  point_diagnostics[[outcome]] <- list(converged = fit$converged,
    warnings = captured$warnings, error = captured$error, rank = fit$rank,
    formula = paste(deparse(formula(fit)), collapse = " "))
}

# Match Figure 2: one seed before the two-outcome loop; fixed centre sizes;
# refit each bootstrap sample, then standardize over that bootstrap sample.
set.seed(2026073011)
B <- 1000L
prob_rows <- list()
diagnostic_rows <- list()
bootstrap_rows <- list()
for (outcome in c("s_ge2", "s_ge3")) {
  point <- gcomp(fits[[outcome]], canonical)
  boots <- matrix(NA_real_, nrow = B, ncol = 5)
  for (b in seq_len(B)) {
    idx <- unlist(lapply(split(seq_len(nrow(canonical)), canonical$center_code), function(ii) {
      sample(ii, length(ii), replace = TRUE)
    }), use.names = FALSE)
    bd <- droplevels(canonical[idx, , drop = FALSE])
    check_centres <- identical(as.integer(table(bd$center_code)), c(440L, 139L, 356L))
    captured <- capture_call({
      bf <- glm(formula_for(outcome), data = bd, family = binomial())
      list(fit = bf, probabilities = gcomp(bf, bd))
    })
    success <- identical(captured$error, "") &&
      length(captured$value$probabilities) == 5L &&
      all(is.finite(captured$value$probabilities))
    if (success) boots[b, ] <- captured$value$probabilities
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
      outcome = outcome, replicate = b, success = success,
      converged = if (is.null(captured$value)) FALSE else captured$value$fit$converged,
      centre_sizes_preserved = check_centres,
      warning_count = length(captured$warnings),
      warnings = paste(captured$warnings, collapse = " | "), error = captured$error,
      stringsAsFactors = FALSE
    )
  }
  ci <- apply(boots, 2, quantile, probs = c(.025, .975), na.rm = TRUE)
  prob_rows[[outcome]] <- data.frame(
    outcome = outcome, cmrf_count = 1:5, probability = point,
    ci_low = ci[1, ], ci_high = ci[2, ],
    bootstrap_success = colSums(!is.na(boots)), bootstrap_replicates = B,
    N = nrow(canonical), events = sum(canonical[[outcome]]),
    adjustment = "Age, sex, centre; marginal standardization",
    stringsAsFactors = FALSE
  )
  bootstrap_rows[[outcome]] <- data.frame(outcome = outcome,
    replicate = rep(seq_len(B), each = 5L), cmrf_count = rep(1:5, B),
    probability = as.vector(t(boots)))
  message(outcome, ": ", min(colSums(!is.na(boots))), "/", B, " bootstrap estimates per count.")
}
probabilities <- do.call(rbind, prob_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
check("all_2000_bootstrap_runs_successful", nrow(diagnostics) == 2000L && all(diagnostics$success))
check("all_2000_bootstrap_models_converged", all(diagnostics$converged))
check("bootstrap_warning_count_zero", sum(diagnostics$warning_count) == 0L)
check("bootstrap_centre_sizes_preserved", all(diagnostics$centre_sizes_preserved))
check("probability_output_10_rows", nrow(probabilities) == 10L)
check("probabilities_and_CI_in_unit_interval", all(as.matrix(probabilities[
  c("probability", "ci_low", "ci_high")]) >= 0) &&
  all(as.matrix(probabilities[c("probability", "ci_low", "ci_high")]) <= 1))
check("confidence_interval_order", all(probabilities$ci_low < probabilities$ci_high))

distribution <- as.data.frame(table(cmrf_count = factor(canonical$cmrf_count, levels = 1:5),
  steatosis_stage = factor(canonical$steatosis, levels = 0:3)))
names(distribution)[3] <- "n"
distribution$cmrf_count <- as.integer(as.character(distribution$cmrf_count))
distribution$steatosis_stage <- as.integer(as.character(distribution$steatosis_stage))
distribution$pct <- ave(distribution$n, distribution$cmrf_count, FUN = function(z) 100 * z / sum(z))
check("distribution_all_grades_counts_including_zero", nrow(distribution) == 20L && sum(distribution$n) == 935L)
check("distribution_percentages_sum_100", all(abs(tapply(distribution$pct,
  distribution$cmrf_count, sum) - 100) < 1e-12))

observed_rows <- list()
for (outcome in c("s_ge2", "s_ge3")) for (count in 1:5) {
  d <- canonical[canonical$cmrf_count == count, , drop = FALSE]
  n <- nrow(d)
  events <- sum(d[[outcome]])
  p <- events / n
  z <- qnorm(.975)
  den <- 1 + z^2 / n
  midpoint <- (p + z^2 / (2 * n)) / den
  halfwidth <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  observed_rows[[length(observed_rows) + 1L]] <- data.frame(
    count = count, outcome = if (outcome == "s_ge2") "S≥2" else "S3",
    N = n, events = events, probability = p, ci_low = midpoint - halfwidth,
    ci_high = midpoint + halfwidth, stringsAsFactors = FALSE)
}
observed <- do.call(rbind, observed_rows)
check("observed_event_totals", identical(as.numeric(tapply(observed$events,
  observed$outcome, sum)[c("S≥2", "S3")]), c(645, 301)))
check("observed_denominators", all(observed$N == rep(c(40, 188, 313, 292, 102), 2)))

write.csv(probabilities, file.path(source_dir, "standardized_steatosis_probabilities.csv"), row.names = FALSE)
write.csv(distribution, file.path(source_dir, "steatosis_distribution_by_cmrf.csv"), row.names = FALSE)
write.csv(observed, file.path(source_dir, "observed_steatosis_wilson.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(do.call(rbind, reproduction), file.path(source_dir, "original_model_reproduction.csv"), row.names = FALSE)
check("original_three_threshold_csv_copied", file.copy(original_csv,
  file.path(source_dir, "steatosis_secondary_models.csv"), overwrite = FALSE))
check("original_three_threshold_csv_byte_identical", identical(hash_file(original_csv),
  hash_file(file.path(source_dir, "steatosis_secondary_models.csv"))))
write.csv(diagnostics, file.path(qa_dir, "steatosis_bootstrap_diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, bootstrap_rows), file.path(qa_dir, "steatosis_bootstrap_probabilities.csv"), row.names = FALSE)
check("protected_inputs_unchanged", identical(input_hashes, vapply(input_paths, hash_file, character(1))))
write.csv(do.call(rbind, checks), file.path(qa_dir, "steatosis_analysis_checks.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(qa_dir, "steatosis_sessionInfo.txt"))
relative_path <- function(path) sub(paste0(project_dir, "/"), "", path, fixed = TRUE)
manifest <- list(
  author = "mjwang", version = version_id, started_at = started_at,
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script = list(path = relative_path(script_path), sha256 = hash_file(script_path)),
  inputs = lapply(names(input_paths), function(name) list(role = name,
    path = relative_path(input_paths[[name]]), sha256 = input_hashes[[name]])),
  analysis_input = "Only frozen canonical_adult_complete_component.csv; original CSV is a numerical comparison source.",
  N = 935L, centres = as.list(table(canonical$center_code)),
  centre_mapping = "ZS to RJ in a derived centre variable; original IDs and source-prefix fields preserved.",
  outcome_events = list(s_ge2 = 645L, s_ge3 = 301L),
  formulas = lapply(point_diagnostics, function(x) x$formula),
  point_model_diagnostics = point_diagnostics,
  standardization = "Assign counts 1:5 in turn to every participant, retain age/sex/centre, average response predictions.",
  uncertainty = "Centre-stratified nonparametric bootstrap; refit and standardize over the resampled cohort; pointwise percentile 95% CIs; quantile type 7.",
  seed = 2026073011L, bootstrap_replicates_per_outcome = B,
  RNG_kind = RNGkind(), bootstrap_outcome_order = c("s_ge2", "s_ge3"),
  bootstrap_diagnostics = list(successful = sum(diagnostics$success),
    converged = sum(diagnostics$converged), warnings = sum(diagnostics$warning_count)),
  original_BH_family = "All three original steatosis thresholds retained unchanged; no additional hypothesis tests.",
  limits = "Cross-sectional prevalent-outcome probabilities, not future risk or causal effects. Complete-component cohort only; no new imputation; no BMI adjustment.",
  software = list(R = R.version.string, digest = as.character(utils::packageVersion("digest")),
    jsonlite = as.character(utils::packageVersion("jsonlite"))),
  checks_passed = length(checks),
  outputs = lapply(outputs[file.exists(outputs)], function(path) list(
    path = relative_path(path), sha256 = hash_file(path)))
)
jsonlite::write_json(manifest, file.path(qa_dir, "steatosis_analysis_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = 16, null = "null")
print(probabilities, row.names = FALSE, digits = 8)
message(length(checks), " analysis checks passed; protected inputs unchanged.")
