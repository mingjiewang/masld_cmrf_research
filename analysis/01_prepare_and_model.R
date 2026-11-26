#!/usr/bin/env Rscript
# Author: mjwang

# Paper 1 full adult re-analysis.
# Reads the read-only source workbook and writes only to version 20260809_152651.

suppressPackageStartupMessages({
  library(readxl)
  library(digest)
  library(MASS)
  library(ordinal)
  library(mice)
  library(splines)
})

version_id <- "20260809_152651"
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
paper_dir <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))
project_dir <- normalizePath(file.path(paper_dir, ".."))
source_xlsx <- file.path(project_dir, "data", "All_patients_SAFE.xlsx")
out_dir <- file.path(paper_dir, "results", "versions", version_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_xlsx)) stop("Missing source workbook: ", source_xlsx)

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}

or_row <- function(fit, term, outcome, analysis, level = term) {
  z <- coef(summary(fit))[term, ]
  data.frame(
    analysis = analysis,
    outcome = outcome,
    term = term,
    level = level,
    N = nobs(fit),
    events = sum(model.response(model.frame(fit))),
    estimate = unname(z["Estimate"]),
    std_error = unname(z["Std. Error"]),
    OR = exp(unname(z["Estimate"])),
    ci_low = exp(unname(z["Estimate"]) - 1.96 * unname(z["Std. Error"])),
    ci_high = exp(unname(z["Estimate"]) + 1.96 * unname(z["Std. Error"])),
    p_value = unname(z["Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )
}

pool_scalar_rubin <- function(q, u, n, k) {
  m <- length(q)
  qbar <- mean(q)
  ubar <- mean(u)
  b <- stats::var(q)
  total <- ubar + (1 + 1 / m) * b
  if (is.na(b) || b <= .Machine$double.eps) {
    df <- Inf
    fmi <- 0
  } else {
    df_old <- (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2
    lambda <- (1 + 1 / m) * b / total
    df_obs <- ((n - k) + 1) / ((n - k) + 3) * (n - k) * (1 - lambda)
    df <- 1 / (1 / df_old + 1 / df_obs)
    fmi <- lambda + (2 / (df + 3)) * (1 - lambda)
  }
  se <- sqrt(total)
  crit <- if (is.finite(df)) qt(0.975, df) else qnorm(0.975)
  p <- if (is.finite(df)) 2 * pt(abs(qbar / se), df, lower.tail = FALSE) else
    2 * pnorm(abs(qbar / se), lower.tail = FALSE)
  c(
    estimate = qbar, std_error = se, df = df, fmi = fmi, p_value = p,
    OR = exp(qbar), ci_low = exp(qbar - crit * se), ci_high = exp(qbar + crit * se)
  )
}

derive_cmrf <- function(d, prefix = "") {
  cmrf1 <- ifelse(
    d$bmi >= 23 |
      (!is.na(d$wc) & d$sex == 0 & d$wc >= 90) |
      (!is.na(d$wc) & d$sex == 1 & d$wc >= 80),
    1L,
    ifelse(d$bmi < 23 & !is.na(d$wc), 0L, NA_integer_)
  )
  cmrf2 <- ifelse(
    d$dm == 1 |
      (!is.na(d$glu) & d$glu >= 5.6) |
      (!is.na(d$hba1c) & d$hba1c >= 5.7),
    1L,
    ifelse(
      d$dm == 0 & !is.na(d$glu) & d$glu < 5.6 &
        !is.na(d$hba1c) & d$hba1c < 5.7,
      0L, NA_integer_
    )
  )
  cmrf3 <- ifelse(is.na(d$htn), NA_integer_, as.integer(d$htn == 1))
  cmrf4 <- ifelse(is.na(d$tg), NA_integer_, as.integer(d$tg >= 1.7))
  cmrf5 <- ifelse(
    is.na(d$hdl), NA_integer_,
    as.integer(
      (d$sex == 0 & d$hdl < 1.0) |
        (d$sex == 1 & d$hdl < 1.3)
    )
  )
  comps <- data.frame(cmrf1, cmrf2, cmrf3, cmrf4, cmrf5)
  count_min <- rowSums(as.data.frame(lapply(comps, function(x) ifelse(is.na(x), 0L, x))))
  count_max <- rowSums(as.data.frame(lapply(comps, function(x) ifelse(is.na(x), 1L, x))))
  count <- ifelse(complete.cases(comps), rowSums(comps), NA_integer_)
  names(comps) <- paste0(prefix, names(comps))
  cbind(comps,
        setNames(data.frame(count_min), paste0(prefix, "cmrf_count_min")),
        setNames(data.frame(count_max), paste0(prefix, "cmrf_count_max")),
        setNames(data.frame(count), paste0(prefix, "cmrf_count")))
}

message("Reading source workbook: ", source_xlsx)
x <- as.data.frame(read_excel(source_xlsx, sheet = "Sheet1"))
source_names <- names(x)
map <- c(
  id = "患者ID", date = "时间", sex = "性别", age = "年龄",
  height = "身高", weight = "体重", bmi = "BMI", wc = "腰围",
  bri = "BRI", htn = "高血压", dm = "糖尿病", alcohol = "饮酒史",
  glu = "GLU", ALT = "ALT", AST = "AST", alp = "ALP", GGT = "GGT",
  tbil = "TBiL", alb = "ALB", glo = "GLO", bun = "BUN", cr = "Cr",
  ua = "UA", tg = "TG", tc = "Tc", hdl = "HDL", ldl = "LDL",
  cp = "CP", hba1c = "HbA1C", insulin = "Insulin", wbc = "WBC",
  hb = "Hb", plt = "PLT", fer = "Fer", afp = "AFP", ana = "ANA",
  hbsag = "HBsAg", hcv = "HCV", lsm = "LSM", cap = "CAP",
  steatosis = "Steatosis", activity = "Activity", fibrosis = "Fibrosis",
  safe = "SAFE", fast = "FAST", nfs = "NFS", fib4 = "FIB4",
  fib8 = "FIB-8", agile3 = "Agile3+ Score", agile4 = "Agile4 Score",
  cmrf_count_legacy = "CMRF_count"
)
if (length(setdiff(unname(map), names(x)))) {
  stop("Missing mapped columns: ", paste(setdiff(unname(map), names(x)), collapse = ", "))
}
raw <- x[unname(map)]
names(raw) <- names(map)

raw$id <- as.character(raw$id)
raw$source_prefix <- factor(substr(raw$id, 1, 2), levels = c("HZ", "JL", "RJ", "ZS"))
raw$center_code <- factor(
  ifelse(as.character(raw$source_prefix) == "ZS", "RJ", as.character(raw$source_prefix)),
  levels = c("HZ", "JL", "RJ")
)
expected_centers <- c("HZ", "JL", "RJ")
if (any(is.na(raw$center_code))) {
  stop("Unexpected center code")
}
raw$sex <- as.integer(raw$sex)
raw$sex_label <- factor(raw$sex, levels = c(0, 1), labels = c("Male", "Female"))

zero_to_missing <- c("height", "weight", "wc", "glu", "hba1c", "tg", "hdl")
for (v in zero_to_missing) {
  raw[[v]] <- as.numeric(raw[[v]])
  raw[[v]][!is.na(raw[[v]]) & raw[[v]] <= 0] <- NA_real_
}
for (v in setdiff(names(raw), c("id", "date", "source_prefix", "center_code", "sex_label"))) {
  raw[[v]] <- suppressWarnings(as.numeric(raw[[v]]))
}

cmrf <- derive_cmrf(raw)
raw <- cbind(raw, cmrf)
raw$adult <- raw$age >= 18
raw$histology_complete <- !is.na(raw$steatosis) & !is.na(raw$fibrosis)
raw$eligible <- raw$adult & raw$histology_complete & raw$cmrf_count_min >= 1
raw$complete_component <- raw$eligible & !is.na(raw$cmrf_count) & raw$cmrf_count >= 1

eligible <- droplevels(raw[raw$eligible, , drop = FALSE])
canonical <- droplevels(raw[raw$complete_component, , drop = FALSE])
all_age_complete <- droplevels(raw[
  raw$histology_complete & !is.na(raw$cmrf_count) & raw$cmrf_count >= 1,
  , drop = FALSE
])

for (dname in c("eligible", "canonical", "all_age_complete")) {
  d <- get(dname)
  d$f_ge2 <- as.integer(d$fibrosis >= 2)
  d$f_ge3 <- as.integer(d$fibrosis >= 3)
  d$s_ge1 <- as.integer(d$steatosis >= 1)
  d$s_ge2 <- as.integer(d$steatosis >= 2)
  d$s_ge3 <- as.integer(d$steatosis >= 3)
  d$cmrf_count_nonadiposity <- rowSums(d[paste0("cmrf", 2:5)], na.rm = FALSE)
  d$female_hdl_increment <- as.integer(
    d$sex == 1 & !is.na(d$hdl) & d$hdl >= 1.0 & d$hdl < 1.3
  )
  d$cmrf_count_unified_hdl <- d$cmrf_count - d$female_hdl_increment
  d$glycaemia_category <- factor(
    ifelse(d$dm == 1, "Recorded T2DM",
           ifelse(d$cmrf2 == 1, "Dysglycaemia without recorded T2DM", "Normoglycaemia")),
    levels = c("Normoglycaemia", "Dysglycaemia without recorded T2DM", "Recorded T2DM")
  )
  assign(dname, d)
}

stopifnot(
  nrow(raw) == 1100,
  nrow(all_age_complete) == 937,
  nrow(canonical) == 935,
  !anyDuplicated(canonical$id),
  all(canonical$cmrf_count == rowSums(canonical[paste0("cmrf", 1:5)]))
)

write.csv(canonical, file.path(out_dir, "canonical_adult_complete_component.csv"), row.names = FALSE)
write.csv(eligible, file.path(out_dir, "eligible_adult_for_mi.csv"), row.names = FALSE)

# Reproducibility and variable definitions -----------------------------------

manifest <- data.frame(
  relative_path = "data/All_patients_SAFE.xlsx",
  bytes = file.info(source_xlsx)$size,
  modified_time = format(file.info(source_xlsx)$mtime, "%Y-%m-%d %H:%M:%S %z"),
  sha256 = digest(source_xlsx, algo = "sha256", file = TRUE),
  rows = nrow(raw),
  source_columns = length(source_names),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_dir, "input_manifest.csv"), row.names = FALSE)

variable_dictionary <- data.frame(
  variable = c("cmrf1", "cmrf2", "cmrf3", "cmrf4", "cmrf5", "cmrf_count",
               "f_ge2", "f_ge3"),
  label = c("Adiposity", "Glycaemia", "Hypertension",
            "Hypertriglyceridaemia", "Low HDL-C", "CMRF burden",
            "Significant fibrosis", "Advanced fibrosis"),
  rule = c(
    "BMI>=23 kg/m2 OR waist>=90 cm in males/>=80 cm in females",
    "FPG>=5.6 mmol/L OR HbA1c>=5.7% OR recorded T2DM",
    "Recorded hypertension; source contains no treatment field",
    "TG>=1.7 mmol/L; source contains no treatment field",
    "HDL-C<1.0 mmol/L in males/<1.3 mmol/L in females",
    "Sum of five components when all are determinable",
    "Fibrosis stage >=2", "Fibrosis stage >=3"
  ),
  indeterminate_rule = c(
    "BMI<23 with missing waist", "No T2DM but FPG or HbA1c missing",
    "Hypertension missing", "TG missing", "HDL-C missing",
    "Any component indeterminate", "Fibrosis missing", "Fibrosis missing"
  ),
  stringsAsFactors = FALSE
)
write.csv(variable_dictionary, file.path(out_dir, "variable_dictionary.csv"), row.names = FALSE)

flow <- data.frame(
  step_order = 1:7,
  step = c(
    "Source records",
    "Adults aged >=18 years",
    "Adults with steatosis and fibrosis available",
    "Adults with at least one definitely positive CMRF",
    "Adult eligible cohort with >=1 indeterminate component",
    "Adult complete-component analytic cohort",
    "All-age complete-component sensitivity cohort"
  ),
  n = c(
    nrow(raw),
    sum(raw$adult, na.rm = TRUE),
    sum(raw$adult & raw$histology_complete, na.rm = TRUE),
    nrow(eligible),
    sum(!complete.cases(eligible[paste0("cmrf", 1:5)])),
    nrow(canonical),
    nrow(all_age_complete)
  ),
  stringsAsFactors = FALSE
)
write.csv(flow, file.path(out_dir, "analysis_population_flow.csv"), row.names = FALSE)

missing_center <- do.call(rbind, lapply(levels(eligible$center_code), function(cc) {
  d <- eligible[eligible$center_code == cc, , drop = FALSE]
  do.call(rbind, lapply(1:5, function(j) {
    data.frame(
      center_code = cc,
      component = c("Adiposity", "Glycaemia", "Hypertension",
                    "Hypertriglyceridaemia", "Low HDL-C")[j],
      N = nrow(d),
      indeterminate_n = sum(is.na(d[[paste0("cmrf", j)]])),
      indeterminate_pct = 100 * mean(is.na(d[[paste0("cmrf", j)]])),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(missing_center, file.path(out_dir, "cmrf_indeterminate_by_center.csv"), row.names = FALSE)

center_mix <- do.call(rbind, lapply(levels(canonical$center_code), function(cc) {
  d <- canonical[canonical$center_code == cc, , drop = FALSE]
  do.call(rbind, lapply(0:4, function(stage) {
    data.frame(
      center_code = cc, fibrosis_stage = stage, N_center = nrow(d),
      n = sum(d$fibrosis == stage),
      pct = 100 * mean(d$fibrosis == stage),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(center_mix, file.path(out_dir, "center_fibrosis_case_mix.csv"), row.names = FALSE)

transition <- as.data.frame(table(
  legacy = canonical$cmrf_count_legacy,
  consensus = canonical$cmrf_count
))
transition <- transition[transition$Freq > 0, ]
names(transition)[3] <- "n"
write.csv(transition, file.path(out_dir, "legacy_consensus_transition_adults.csv"), row.names = FALSE)

# Descriptive Table 1 --------------------------------------------------------

fmt_mean_sd <- function(x) sprintf("%.1f ± %.1f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_med_iqr <- function(x) {
  q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
  sprintf("%.2f [%.2f–%.2f]", q[2], q[1], q[3])
}
fmt_n_pct <- function(x, value = 1) {
  den <- sum(!is.na(x))
  sprintf("%d/%d (%.1f%%)", sum(x == value, na.rm = TRUE), den,
          100 * mean(x == value, na.rm = TRUE))
}

groups <- c("Overall", paste0("CMRF ", 1:5))
subset_group <- function(g) if (g == "Overall") canonical else
  canonical[canonical$cmrf_count == as.integer(sub("CMRF ", "", g)), , drop = FALSE]
table1_rows <- list()
add_t1 <- function(label, variable, formatter, value = NULL) {
  z <- data.frame(variable = label, stringsAsFactors = FALSE)
  for (g in groups) {
    d <- subset_group(g)
    z[[g]] <- if (is.null(value)) formatter(d[[variable]]) else formatter(d[[variable]], value)
  }
  z$missing_n_overall <- sum(is.na(canonical[[variable]]))
  table1_rows[[length(table1_rows) + 1]] <<- z
}
add_t1("Age, years", "age", fmt_mean_sd)
add_t1("Female sex", "sex", fmt_n_pct, 1)
add_t1("BMI, kg/m2", "bmi", fmt_mean_sd)
add_t1("Waist circumference, cm", "wc", fmt_med_iqr)
add_t1("Fasting glucose, mmol/L", "glu", fmt_med_iqr)
add_t1("HbA1c, %", "hba1c", fmt_med_iqr)
add_t1("Triglycerides, mmol/L", "tg", fmt_med_iqr)
add_t1("HDL-C, mmol/L", "hdl", fmt_med_iqr)
add_t1("ALT, U/L", "ALT", fmt_med_iqr)
add_t1("AST, U/L", "AST", fmt_med_iqr)
add_t1("Platelets, 10^9/L", "plt", fmt_med_iqr)
add_t1("Adiposity", "cmrf1", fmt_n_pct, 1)
add_t1("Glycaemia", "cmrf2", fmt_n_pct, 1)
add_t1("Hypertension", "cmrf3", fmt_n_pct, 1)
add_t1("Hypertriglyceridaemia", "cmrf4", fmt_n_pct, 1)
add_t1("Low HDL-C", "cmrf5", fmt_n_pct, 1)
add_t1("Significant fibrosis (F>=2)", "f_ge2", fmt_n_pct, 1)
add_t1("Advanced fibrosis (F>=3)", "f_ge3", fmt_n_pct, 1)
table1 <- do.call(rbind, table1_rows)
write.csv(table1, file.path(out_dir, "Table1_descriptive.csv"), row.names = FALSE)

# Complete-component primary and key-secondary models -----------------------

main_rows <- list()
functional_rows <- list()
categorical_rows <- list()
canonical$cmrf_factor <- relevel(factor(canonical$cmrf_count), ref = "1")

for (outcome in c("f_ge2", "f_ge3")) {
  fit <- glm(
    reformulate(c("cmrf_count", "age", "sex_label", "center_code"), outcome),
    data = canonical, family = binomial()
  )
  main_rows[[length(main_rows) + 1]] <- or_row(
    fit, "cmrf_count", outcome, "Adult complete-component"
  )
  fit_factor <- glm(
    reformulate(c("cmrf_factor", "age", "sex_label", "center_code"), outcome),
    data = canonical, family = binomial()
  )
  fit_spline <- glm(
    reformulate(c("cmrf_count", "ns(age, 3)", "sex_label", "center_code"), outcome),
    data = canonical, family = binomial()
  )
  fit_interaction <- glm(
    reformulate(c("cmrf_count * center_code", "age", "sex_label"), outcome),
    data = canonical, family = binomial()
  )
  functional_rows[[length(functional_rows) + 1]] <- data.frame(
    outcome = outcome,
    linear_vs_categorical_LRT_p = anova(fit, fit_factor, test = "LRT")$`Pr(>Chi)`[2],
    linear_age_vs_spline_LRT_p = anova(fit, fit_spline, test = "LRT")$`Pr(>Chi)`[2],
    center_interaction_LRT_p = anova(fit, fit_interaction, test = "LRT")$`Pr(>Chi)`[2],
    stringsAsFactors = FALSE
  )
  for (term in grep("^cmrf_factor", names(coef(fit_factor)), value = TRUE)) {
    categorical_rows[[length(categorical_rows) + 1]] <- or_row(
      fit_factor, term, outcome, "CMRF burden categorical",
      level = sub("cmrf_factor", "", term)
    )
  }
}
main_cc <- do.call(rbind, main_rows)
write.csv(main_cc, file.path(out_dir, "main_models_complete_component.csv"), row.names = FALSE)
write.csv(do.call(rbind, functional_rows), file.path(out_dir, "functional_form_tests.csv"), row.names = FALSE)
write.csv(do.call(rbind, categorical_rows), file.path(out_dir, "cmrf_categorical_models.csv"), row.names = FALSE)

# Centre-stratified bootstrap standardized probabilities --------------------

gcomp <- function(fit, d, values = 1:5) {
  vapply(values, function(v) {
    nd <- d
    nd$cmrf_count <- v
    mean(predict(fit, newdata = nd, type = "response"))
  }, numeric(1))
}

set.seed(2026073011)
B <- 1000L
prob_rows <- list()
for (outcome in c("f_ge2", "f_ge3")) {
  fit <- glm(
    reformulate(c("cmrf_count", "age", "sex_label", "center_code"), outcome),
    data = canonical, family = binomial()
  )
  point <- gcomp(fit, canonical)
  boots <- matrix(NA_real_, nrow = B, ncol = 5)
  for (b in seq_len(B)) {
    idx <- unlist(lapply(split(seq_len(nrow(canonical)), canonical$center_code), function(ii) {
      sample(ii, length(ii), replace = TRUE)
    }), use.names = FALSE)
    bd <- droplevels(canonical[idx, , drop = FALSE])
    bf <- try(glm(
      reformulate(c("cmrf_count", "age", "sex_label", "center_code"), outcome),
      data = bd, family = binomial()
    ), silent = TRUE)
    if (!inherits(bf, "try-error")) boots[b, ] <- gcomp(bf, bd)
  }
  ci <- apply(boots, 2, quantile, probs = c(.025, .975), na.rm = TRUE)
  prob_rows[[length(prob_rows) + 1]] <- data.frame(
    outcome = outcome, cmrf_count = 1:5, probability = point,
    ci_low = ci[1, ], ci_high = ci[2, ],
    bootstrap_success = colSums(!is.na(boots)), bootstrap_replicates = B,
    N = nrow(canonical),
    events = sum(canonical[[outcome]]),
    adjustment = "Age, sex, centre; marginal standardization",
    stringsAsFactors = FALSE
  )
}
probability_results <- do.call(rbind, prob_rows)
write.csv(probability_results, file.path(out_dir, "standardized_probabilities.csv"), row.names = FALSE)

fibrosis_distribution <- as.data.frame(table(
  cmrf_count = canonical$cmrf_count,
  fibrosis_stage = canonical$fibrosis
))
names(fibrosis_distribution)[3] <- "n"
fibrosis_distribution$pct <- ave(
  fibrosis_distribution$n, fibrosis_distribution$cmrf_count,
  FUN = function(z) 100 * z / sum(z)
)
write.csv(fibrosis_distribution, file.path(out_dir, "fibrosis_distribution_by_cmrf.csv"), row.names = FALSE)

# Component and glycaemia complete-component models -------------------------

component_rows <- list()
glycaemia_rows <- list()
for (outcome in c("f_ge2", "f_ge3")) {
  fit <- glm(
    reformulate(c(paste0("cmrf", 1:5), "age", "sex_label", "center_code"), outcome),
    data = canonical, family = binomial()
  )
  for (term in paste0("cmrf", 1:5)) {
    component_rows[[length(component_rows) + 1]] <- or_row(
      fit, term, outcome, "Adult complete-component; mutually adjusted components"
    )
  }
  gfit <- glm(
    reformulate(
      c("glycaemia_category", "cmrf1", "cmrf3", "cmrf4", "cmrf5",
        "age", "sex_label", "center_code"),
      outcome
    ),
    data = canonical, family = binomial()
  )
  for (term in grep("^glycaemia_category", names(coef(gfit)), value = TRUE)) {
    glycaemia_rows[[length(glycaemia_rows) + 1]] <- or_row(
      gfit, term, outcome, "Adult complete-component; glycaemia decomposition",
      level = sub("^glycaemia_category", "", term)
    )
  }
}
component_cc <- do.call(rbind, component_rows)
component_cc$p_fdr_all_10 <- p.adjust(component_cc$p_value, method = "BH")
glycaemia_cc <- do.call(rbind, glycaemia_rows)
glycaemia_cc$p_holm_all_4 <- p.adjust(glycaemia_cc$p_value, method = "holm")
write.csv(component_cc, file.path(out_dir, "component_models_complete_component.csv"), row.names = FALSE)
write.csv(glycaemia_cc, file.path(out_dir, "glycaemia_decomposition_complete_component.csv"), row.names = FALSE)

# Raw-variable-informed multiple imputation ---------------------------------

mi_vars <- c(
  "wc", "glu", "hba1c", "tg", "hdl",
  "age", "sex", "center_code", "bmi", "dm", "htn",
  "steatosis", "fibrosis", "ALT", "AST", "plt", "alcohol", "hcv"
)
mi_data <- eligible[mi_vars]
mi_data$center_code <- factor(mi_data$center_code)
mi_data$sex <- factor(mi_data$sex, levels = c(0, 1))
mi_data$dm <- factor(mi_data$dm, levels = c(0, 1))
mi_data$htn <- factor(mi_data$htn, levels = c(0, 1))
mi_data$alcohol <- factor(mi_data$alcohol, levels = c(0, 1))
mi_data$hcv <- factor(mi_data$hcv, levels = c(0, 1))

mi_method <- make.method(mi_data)
mi_method[] <- ""
mi_method[c("wc", "glu", "hba1c", "tg", "hdl")] <- "pmm"
mi_method["alcohol"] <- "logreg"
mi_method["hcv"] <- "logreg"
pred <- make.predictorMatrix(mi_data)
diag(pred) <- 0
pred[c("age", "sex", "center_code", "bmi", "dm", "htn", "steatosis", "fibrosis",
       "ALT", "AST", "plt"), ] <- 0

set.seed(2026073012)
imp <- mice(
  mi_data, m = 50, maxit = 20, method = mi_method,
  predictorMatrix = pred, printFlag = FALSE
)
saveRDS(imp, file.path(out_dir, "multiple_imputation_object.rds"))

mi_missing <- data.frame(
  variable = names(mi_data),
  missing_n = vapply(mi_data, function(z) sum(is.na(z)), integer(1)),
  missing_pct = 100 * vapply(mi_data, function(z) mean(is.na(z)), numeric(1)),
  method = mi_method,
  stringsAsFactors = FALSE
)
write.csv(mi_missing, file.path(out_dir, "multiple_imputation_missingness.csv"), row.names = FALSE)

mi_main_fit <- list(f_ge2 = list(), f_ge3 = list())
mi_component_fit <- list(f_ge2 = setNames(vector("list", 5), paste0("cmrf", 1:5)),
                         f_ge3 = setNames(vector("list", 5), paste0("cmrf", 1:5)))
mi_gly_fit <- list(f_ge2 = list(), f_ge3 = list())

for (i in seq_len(imp$m)) {
  d <- complete(imp, i)
  for (v in c("sex", "dm", "htn", "alcohol", "hcv")) d[[v]] <- as.integer(as.character(d[[v]]))
  d$sex_label <- factor(d$sex, levels = c(0, 1), labels = c("Male", "Female"))
  dc <- derive_cmrf(d)
  d <- cbind(d, dc)
  if (anyNA(d[paste0("cmrf", 1:5)])) stop("MI completed data still contains indeterminate CMRF")
  d$f_ge2 <- as.integer(d$fibrosis >= 2)
  d$f_ge3 <- as.integer(d$fibrosis >= 3)
  d$glycaemia_category <- factor(
    ifelse(d$dm == 1, "Recorded T2DM",
           ifelse(d$cmrf2 == 1, "Dysglycaemia without recorded T2DM", "Normoglycaemia")),
    levels = c("Normoglycaemia", "Dysglycaemia without recorded T2DM", "Recorded T2DM")
  )
  for (outcome in c("f_ge2", "f_ge3")) {
    mf <- glm(
      reformulate(c("cmrf_count", "age", "sex_label", "center_code"), outcome),
      data = d, family = binomial()
    )
    mi_main_fit[[outcome]][[i]] <- c(
      q = unname(coef(mf)["cmrf_count"]),
      u = unname(vcov(mf)["cmrf_count", "cmrf_count"])
    )
    cf <- glm(
      reformulate(c(paste0("cmrf", 1:5), "age", "sex_label", "center_code"), outcome),
      data = d, family = binomial()
    )
    for (term in paste0("cmrf", 1:5)) {
      mi_component_fit[[outcome]][[term]][[i]] <- c(
        q = unname(coef(cf)[term]), u = unname(vcov(cf)[term, term])
      )
    }
    gf <- glm(
      reformulate(c("glycaemia_category", "cmrf1", "cmrf3", "cmrf4", "cmrf5",
                    "age", "sex_label", "center_code"), outcome),
      data = d, family = binomial()
    )
    for (term in grep("^glycaemia_category", names(coef(gf)), value = TRUE)) {
      mi_gly_fit[[outcome]][[term]][[i]] <- c(
        q = unname(coef(gf)[term]), u = unname(vcov(gf)[term, term])
      )
    }
  }
}

pool_fit_list <- function(x, outcome, analysis, term, level = term, n, k) {
  z <- do.call(rbind, x)
  p <- pool_scalar_rubin(z[, "q"], z[, "u"], n = n, k = k)
  data.frame(
    analysis = analysis, outcome = outcome, term = term, level = level,
    N = n, m = nrow(z), events = sum(eligible[[outcome]]),
    estimate = p["estimate"], std_error = p["std_error"], df = p["df"],
    fmi = p["fmi"], OR = p["OR"], ci_low = p["ci_low"],
    ci_high = p["ci_high"], p_value = p["p_value"],
    stringsAsFactors = FALSE
  )
}

mi_main_rows <- list()
mi_component_rows <- list()
mi_gly_rows <- list()
for (outcome in c("f_ge2", "f_ge3")) {
  mi_main_rows[[length(mi_main_rows) + 1]] <- pool_fit_list(
    mi_main_fit[[outcome]], outcome, "Adult eligible; raw-variable MI",
    "cmrf_count", n = nrow(eligible), k = 7
  )
  for (term in paste0("cmrf", 1:5)) {
    mi_component_rows[[length(mi_component_rows) + 1]] <- pool_fit_list(
      mi_component_fit[[outcome]][[term]], outcome,
      "Adult eligible; raw-variable MI; mutually adjusted components",
      term, n = nrow(eligible), k = 11
    )
  }
  for (term in names(mi_gly_fit[[outcome]])) {
    mi_gly_rows[[length(mi_gly_rows) + 1]] <- pool_fit_list(
      mi_gly_fit[[outcome]][[term]], outcome,
      "Adult eligible; raw-variable MI; glycaemia decomposition",
      term, level = sub("^glycaemia_category", "", term),
      n = nrow(eligible), k = 12
    )
  }
}
mi_main <- do.call(rbind, mi_main_rows)
mi_components <- do.call(rbind, mi_component_rows)
mi_components$p_fdr_all_10 <- p.adjust(mi_components$p_value, method = "BH")
mi_gly <- do.call(rbind, mi_gly_rows)
mi_gly$p_holm_all_4 <- p.adjust(mi_gly$p_value, method = "holm")
write.csv(mi_main, file.path(out_dir, "main_models_multiple_imputation.csv"), row.names = FALSE)
write.csv(mi_components, file.path(out_dir, "component_models_multiple_imputation.csv"), row.names = FALSE)
write.csv(mi_gly, file.path(out_dir, "glycaemia_decomposition_multiple_imputation.csv"), row.names = FALSE)

main_combined <- rbind(
  transform(main_cc, m = NA_integer_, df = NA_real_, fmi = NA_real_),
  mi_main[names(transform(main_cc, m = NA_integer_, df = NA_real_, fmi = NA_real_))]
)
write.csv(main_combined, file.path(out_dir, "Table2_main_models.csv"), row.names = FALSE)

# Robustness analyses --------------------------------------------------------

robust_rows <- list()
add_robust <- function(d, outcome, exposure, rhs, label) {
  fit <- glm(reformulate(rhs, outcome), data = droplevels(d), family = binomial())
  z <- or_row(fit, exposure, outcome, label)
  robust_rows[[length(robust_rows) + 1]] <<- z
}

for (outcome in c("f_ge2", "f_ge3")) {
  add_robust(canonical, outcome, "cmrf_count",
             c("cmrf_count", "age", "sex_label", "center_code"), "Adult primary")
  add_robust(all_age_complete, outcome, "cmrf_count",
             c("cmrf_count", "age", "sex_label", "center_code"), "All ages")
  d_s1 <- canonical[canonical$steatosis >= 1, , drop = FALSE]
  add_robust(d_s1, outcome, "cmrf_count",
             c("cmrf_count", "age", "sex_label", "center_code"), "Restricted to S>=1")
  add_robust(canonical, outcome, "cmrf_count_legacy",
             c("cmrf_count_legacy", "age", "sex_label", "center_code"), "Legacy CMRF count")
  add_robust(canonical, outcome, "cmrf_count_unified_hdl",
             c("cmrf_count_unified_hdl", "age", "sex_label", "center_code"),
             "Unified HDL-C <1.0")
  add_robust(canonical, outcome, "cmrf_count_nonadiposity",
             c("cmrf_count_nonadiposity", "bmi", "age", "sex_label", "center_code"),
             "Non-adiposity count + BMI")
  d_alc <- canonical[!is.na(canonical$alcohol), , drop = FALSE]
  add_robust(d_alc, outcome, "cmrf_count",
             c("cmrf_count", "age", "sex_label", "center_code", "factor(alcohol)"),
             "Known alcohol + adjustment")
  d_hcv <- canonical[is.na(canonical$hcv) | canonical$hcv == 0, , drop = FALSE]
  add_robust(d_hcv, outcome, "cmrf_count",
             c("cmrf_count", "age", "sex_label", "center_code"), "Exclude non-zero HCV")
  add_robust(canonical, outcome, "cmrf_count",
             c("cmrf_count", "ns(age, 3)", "sex_label", "center_code"), "Age natural spline")
  for (cc in levels(canonical$center_code)) {
    d <- droplevels(canonical[canonical$center_code != cc, , drop = FALSE])
    add_robust(d, outcome, "cmrf_count",
               c("cmrf_count", "age", "sex_label", "center_code"),
               paste0("Leave out centre ", cc))
  }
  dmin <- eligible
  dmin$cmrf_count_scenario <- dmin$cmrf_count_min
  add_robust(dmin, outcome, "cmrf_count_scenario",
             c("cmrf_count_scenario", "age", "sex_label", "center_code"),
             "All indeterminate negative")
  dmax <- eligible
  dmax$cmrf_count_scenario <- dmax$cmrf_count_max
  add_robust(dmax, outcome, "cmrf_count_scenario",
             c("cmrf_count_scenario", "age", "sex_label", "center_code"),
             "All indeterminate positive")
}
robustness <- do.call(rbind, robust_rows)
mi_robust <- mi_main
mi_robust$analysis <- "Adult raw-variable MI"
robustness <- rbind(
  robustness,
  mi_robust[names(robustness)]
)
scenario_order <- c(
  "Adult primary", "Adult raw-variable MI", "All ages", "Restricted to S>=1",
  "Legacy CMRF count", "Unified HDL-C <1.0", "Non-adiposity count + BMI",
  "Known alcohol + adjustment", "Exclude non-zero HCV", "Age natural spline",
  paste0("Leave out centre ", expected_centers),
  "All indeterminate negative", "All indeterminate positive"
)
robustness$scenario_order <- match(robustness$analysis, scenario_order)
robustness <- robustness[order(robustness$outcome, robustness$scenario_order), ]
write.csv(robustness, file.path(out_dir, "robustness_models.csv"), row.names = FALSE)

# Supportive ordinal and exploratory steatosis -------------------------------

canonical$fibrosis_ord <- ordered(canonical$fibrosis)
ord_fit <- polr(
  fibrosis_ord ~ cmrf_count + age + sex_label + center_code,
  data = canonical, Hess = TRUE, method = "logistic"
)
ord_z <- coef(summary(ord_fit))["cmrf_count", ]
ordinal_result <- data.frame(
  N = nobs(ord_fit), outcome = "F0-F4 ordinal fibrosis",
  OR = exp(ord_z["Value"]),
  ci_low = exp(ord_z["Value"] - 1.96 * ord_z["Std. Error"]),
  ci_high = exp(ord_z["Value"] + 1.96 * ord_z["Std. Error"]),
  p_value = 2 * pnorm(abs(ord_z["t value"]), lower.tail = FALSE),
  stringsAsFactors = FALSE
)
clm_fit <- clm(
  fibrosis_ord ~ cmrf_count + age + sex_label + center_code,
  data = canonical, link = "logit", Hess = TRUE
)
po_test <- as.data.frame(nominal_test(clm_fit))
po_test$term <- rownames(po_test)
rownames(po_test) <- NULL
write.csv(ordinal_result, file.path(out_dir, "ordinal_fibrosis_supportive.csv"), row.names = FALSE)
write.csv(po_test, file.path(out_dir, "ordinal_fibrosis_proportional_odds_test.csv"), row.names = FALSE)

steatosis_rows <- list()
for (outcome in c("s_ge1", "s_ge2", "s_ge3")) {
  fit <- glm(
    reformulate(c("cmrf_count", "age", "sex_label", "center_code"), outcome),
    data = canonical, family = binomial()
  )
  steatosis_rows[[length(steatosis_rows) + 1]] <- or_row(
    fit, "cmrf_count", outcome, "Exploratory steatosis threshold"
  )
}
steatosis_results <- do.call(rbind, steatosis_rows)
steatosis_results$p_fdr_3 <- p.adjust(steatosis_results$p_value, method = "BH")
write.csv(steatosis_results, file.path(out_dir, "steatosis_exploratory.csv"), row.names = FALSE)

# Final machine-readable summaries ------------------------------------------

all_components <- rbind(
  transform(component_cc, dataset = "Complete-component"),
  transform(mi_components, events = events, dataset = "Multiple imputation")[
    names(transform(component_cc, dataset = "Complete-component"))
  ]
)
write.csv(all_components, file.path(out_dir, "Table3_component_models.csv"), row.names = FALSE)

summary_metrics <- data.frame(
  metric = c(
    "source_n", "adult_n", "adult_histology_complete_n", "adult_eligible_n",
    "adult_complete_component_n", "all_age_complete_component_n",
    "adult_primary_f_ge2_events", "adult_key_secondary_f_ge3_events",
    "adult_adiposity_negative_n", "tg_above_20_n", "female_hdl_increment_n"
  ),
  value = c(
    nrow(raw), sum(raw$adult), sum(raw$adult & raw$histology_complete),
    nrow(eligible), nrow(canonical), nrow(all_age_complete),
    sum(canonical$f_ge2), sum(canonical$f_ge3),
    sum(canonical$cmrf1 == 0), sum(canonical$tg > 20, na.rm = TRUE),
    sum(canonical$female_hdl_increment)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_metrics, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Completed Paper 1 re-analysis: ", out_dir)
