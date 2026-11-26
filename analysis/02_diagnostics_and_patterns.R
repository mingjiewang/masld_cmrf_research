#!/usr/bin/env Rscript
# Author: mjwang

# Paper 1 model QC and explicitly exploratory analyses.
# Reads the frozen adult canonical data produced by 01_reanalysis.R.

suppressPackageStartupMessages({
  library(emmeans)
})

version_id <- "20260809_152651"
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
paper_dir <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))
out_dir <- file.path(paper_dir, "results", "versions", version_id)
canonical_path <- file.path(out_dir, "canonical_adult_complete_component.csv")
if (!file.exists(canonical_path)) stop("Run 01_reanalysis.R first.")

d <- read.csv(canonical_path, check.names = FALSE, stringsAsFactors = FALSE)
d$center_code <- factor(d$center_code)
d$sex_label <- factor(d$sex_label, levels = c("Male", "Female"))
d$cmrf_pattern <- apply(
  d[paste0("cmrf", 1:5)], 1,
  function(z) paste0(c("A", "G", "H", "T", "L")[z == 1], collapse = "+")
)

derive_cmrf_for_mi <- function(x) {
  cmrf1 <- ifelse(
    x$bmi >= 23 |
      (!is.na(x$wc) & x$sex == 0 & x$wc >= 90) |
      (!is.na(x$wc) & x$sex == 1 & x$wc >= 80),
    1L, 0L
  )
  cmrf2 <- ifelse(
    x$dm == 1 |
      (!is.na(x$glu) & x$glu >= 5.6) |
      (!is.na(x$hba1c) & x$hba1c >= 5.7),
    1L, 0L
  )
  cmrf3 <- as.integer(x$htn == 1)
  cmrf4 <- as.integer(x$tg >= 1.7)
  cmrf5 <- as.integer(
    (x$sex == 0 & x$hdl < 1.0) |
      (x$sex == 1 & x$hdl < 1.3)
  )
  data.frame(
    cmrf1 = cmrf1, cmrf2 = cmrf2, cmrf3 = cmrf3,
    cmrf4 = cmrf4, cmrf5 = cmrf5,
    cmrf_count = cmrf1 + cmrf2 + cmrf3 + cmrf4 + cmrf5
  )
}

extract_interaction <- function(fit, interaction_pattern, analysis, outcome) {
  ct <- coef(summary(fit))
  term <- grep(interaction_pattern, rownames(ct), value = TRUE)
  if (!length(term)) stop("Interaction term not found: ", interaction_pattern)
  term <- term[1]
  data.frame(
    analysis = analysis,
    outcome = outcome,
    N = nobs(fit),
    events = sum(model.response(model.frame(fit))),
    term = term,
    interaction_OR = exp(ct[term, "Estimate"]),
    ci_low = exp(ct[term, "Estimate"] - 1.96 * ct[term, "Std. Error"]),
    ci_high = exp(ct[term, "Estimate"] + 1.96 * ct[term, "Std. Error"]),
    p_value = ct[term, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

# BMI interaction is exploratory and distinct from the non-adiposity+BMI estimand.
d$bmi_lt25 <- factor(ifelse(d$bmi < 25, "<25", ">=25"), levels = c(">=25", "<25"))
d$bmi_lt23 <- factor(ifelse(d$bmi < 23, "<23", ">=23"), levels = c(">=23", "<23"))
bmi_rows <- list()
for (outcome in c("f_ge2", "f_ge3")) {
  fit_cont <- glm(
    reformulate(c("cmrf_count * bmi", "age", "sex_label", "center_code"), outcome),
    data = d, family = binomial()
  )
  fit_25 <- glm(
    reformulate(c("cmrf_count * bmi_lt25", "age", "sex_label", "center_code"), outcome),
    data = d, family = binomial()
  )
  fit_23 <- glm(
    reformulate(c("cmrf_count * bmi_lt23", "age", "sex_label", "center_code"), outcome),
    data = d, family = binomial()
  )
  bmi_rows[[length(bmi_rows) + 1]] <- extract_interaction(
    fit_cont, "cmrf_count:bmi", "Continuous BMI interaction", outcome
  )
  bmi_rows[[length(bmi_rows) + 1]] <- extract_interaction(
    fit_25, "cmrf_count:bmi_lt25", "BMI <25 kg/m2 interaction", outcome
  )
  bmi_rows[[length(bmi_rows) + 1]] <- extract_interaction(
    fit_23, "cmrf_count:bmi_lt23", "BMI <23 kg/m2 interaction", outcome
  )
}
bmi_results <- do.call(rbind, bmi_rows)
bmi_results$p_fdr_6 <- p.adjust(bmi_results$p_value, method = "BH")
write.csv(bmi_results, file.path(out_dir, "bmi_interaction_exploratory.csv"),
          row.names = FALSE)

# Pairwise pattern comparisons only within the same CMRF count and for
# patterns represented by at least 30 participants.
pattern_counts <- aggregate(
  id ~ cmrf_count + cmrf_pattern, data = d, FUN = length
)
names(pattern_counts)[3] <- "n"
write.csv(pattern_counts, file.path(out_dir, "pattern_counts_adults.csv"),
          row.names = FALSE)

pairwise_rows <- list()
for (outcome in c("f_ge2", "f_ge3")) {
  for (count in sort(unique(d$cmrf_count))) {
    keep_patterns <- pattern_counts$cmrf_pattern[
      pattern_counts$cmrf_count == count & pattern_counts$n >= 30
    ]
    if (length(keep_patterns) < 2) next
    ds <- droplevels(d[
      d$cmrf_count == count & d$cmrf_pattern %in% keep_patterns,
      , drop = FALSE
    ])
    ds$cmrf_pattern <- factor(ds$cmrf_pattern)
    fit <- glm(
      reformulate(c("cmrf_pattern", "age", "sex_label", "center_code"), outcome),
      data = ds, family = binomial()
    )
    contrast <- as.data.frame(
      summary(pairs(emmeans(fit, ~cmrf_pattern)), infer = TRUE, type = "response")
    )
    if (!nrow(contrast)) next
    pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
      outcome = outcome,
      cmrf_count = count,
      contrast = contrast$contrast,
      N_stratum = nrow(ds),
      odds_ratio = contrast$odds.ratio,
      ci_low = contrast$asymp.LCL,
      ci_high = contrast$asymp.UCL,
      p_value = contrast$p.value,
      stringsAsFactors = FALSE
    )
  }
}
if (length(pairwise_rows)) {
  pairwise_results <- do.call(rbind, pairwise_rows)
  pairwise_results$p_fdr_all <- p.adjust(pairwise_results$p_value, method = "BH")
} else {
  pairwise_results <- data.frame()
}
write.csv(pairwise_results, file.path(out_dir, "pattern_pairwise_exploratory.csv"),
          row.names = FALSE)

# Main-model influence, event/parameter and quasi-separation diagnostics.
diagnostic_rows <- list()
influence_rows <- list()
for (outcome in c("f_ge2", "f_ge3")) {
  fit <- glm(
    reformulate(c("cmrf_count", "age", "sex_label", "center_code"), outcome),
    data = d, family = binomial()
  )
  mf <- model.frame(fit)
  hats <- hatvalues(fit)
  cooks <- cooks.distance(fit)
  std_res <- rstandard(fit)
  prob <- fitted(fit)
  parameter_n <- length(coef(fit))
  events <- sum(model.response(mf))
  nonevents <- nrow(mf) - events
  diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
    outcome = outcome,
    N = nrow(mf),
    events = events,
    nonevents = nonevents,
    parameters = parameter_n,
    events_per_parameter = events / parameter_n,
    nonevents_per_parameter = nonevents / parameter_n,
    max_hat = max(hats),
    hat_above_2p_over_n = sum(hats > 2 * parameter_n / nrow(mf)),
    max_cooks_distance = max(cooks),
    cooks_above_4_over_n = sum(cooks > 4 / nrow(mf)),
    abs_standardized_residual_above_3 = sum(abs(std_res) > 3),
    fitted_probability_below_0_01 = sum(prob < 0.01),
    fitted_probability_above_0_99 = sum(prob > 0.99),
    converged = fit$converged,
    stringsAsFactors = FALSE
  )
  top <- order(cooks, decreasing = TRUE)[seq_len(min(10, length(cooks)))]
  influence_rows[[length(influence_rows) + 1]] <- data.frame(
    outcome = outcome,
    id = d$id[top],
    center_code = d$center_code[top],
    cmrf_count = d$cmrf_count[top],
    fibrosis = d$fibrosis[top],
    hat = hats[top],
    cooks_distance = cooks[top],
    standardized_residual = std_res[top],
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, diagnostic_rows),
          file.path(out_dir, "main_model_diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, influence_rows),
          file.path(out_dir, "main_model_top_influence.csv"), row.names = FALSE)

# Multiple-imputation logged events are preserved, including an explicit zero row.
imp <- readRDS(file.path(out_dir, "multiple_imputation_object.rds"))
if (is.null(imp$loggedEvents) || !nrow(imp$loggedEvents)) {
  logged <- data.frame(
    iteration = NA_integer_, imputation = NA_integer_,
    outcome = NA_character_, method = NA_character_,
    event = "No logged events", stringsAsFactors = FALSE
  )
} else {
  logged <- imp$loggedEvents
}
write.csv(logged, file.path(out_dir, "multiple_imputation_logged_events.csv"),
          row.names = FALSE)

# Observed-versus-imputed distributions and Monte Carlo precision diagnostics.
continuous_mi_vars <- c("wc", "glu", "hba1c", "tg", "hdl")
distribution_rows <- list()
for (v in continuous_mi_vars) {
  observed <- imp$data[[v]][!is.na(imp$data[[v]])]
  imputed <- as.numeric(unlist(imp$imp[[v]], use.names = FALSE))
  distribution_rows[[length(distribution_rows) + 1]] <- data.frame(
    variable = v, source = "Observed", value = observed,
    stringsAsFactors = FALSE
  )
  distribution_rows[[length(distribution_rows) + 1]] <- data.frame(
    variable = v, source = "Imputed", value = imputed,
    stringsAsFactors = FALSE
  )
}
mi_distributions <- do.call(rbind, distribution_rows)
write.csv(
  mi_distributions,
  file.path(out_dir, "multiple_imputation_observed_imputed_distributions.csv"),
  row.names = FALSE
)

mi_precision_rows <- list()
for (outcome in c("f_ge2", "f_ge3")) {
  q <- u <- numeric(imp$m)
  for (i in seq_len(imp$m)) {
    completed <- mice::complete(imp, i)
    for (v in c("sex", "dm", "htn", "alcohol", "hcv")) {
      completed[[v]] <- as.integer(as.character(completed[[v]]))
    }
    completed$sex_label <- factor(
      completed$sex, levels = c(0, 1), labels = c("Male", "Female")
    )
    completed$f_ge2 <- as.integer(completed$fibrosis >= 2)
    completed$f_ge3 <- as.integer(completed$fibrosis >= 3)
    completed <- cbind(completed, derive_cmrf_for_mi(completed))
    fit <- glm(
      reformulate(c("cmrf_count", "age", "sex_label", "center_code"), outcome),
      data = completed, family = binomial()
    )
    q[i] <- unname(coef(fit)["cmrf_count"])
    u[i] <- unname(vcov(fit)["cmrf_count", "cmrf_count"])
  }
  between_variance <- stats::var(q)
  within_variance <- mean(u)
  total_variance <- within_variance + (1 + 1 / imp$m) * between_variance
  lambda <- (1 + 1 / imp$m) * between_variance / total_variance
  monte_carlo_error <- sqrt(between_variance / imp$m)
  mi_precision_rows[[length(mi_precision_rows) + 1]] <- data.frame(
    outcome = outcome,
    m = imp$m,
    pooled_log_or = mean(q),
    pooled_se = sqrt(total_variance),
    within_variance = within_variance,
    between_variance = between_variance,
    fraction_missing_information_approx = lambda,
    monte_carlo_error_log_or = monte_carlo_error,
    monte_carlo_error_pct_pooled_se = 100 * monte_carlo_error / sqrt(total_variance),
    stringsAsFactors = FALSE
  )
}
write.csv(
  do.call(rbind, mi_precision_rows),
  file.path(out_dir, "multiple_imputation_precision_diagnostics.csv"),
  row.names = FALSE
)

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_qc_exploratory.txt"))
message("Completed Paper 1 QC/exploratory analyses: ", out_dir)
