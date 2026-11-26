# Revision author: mjwang
# Derived from frozen three-centre figure code; new MI overrides, unchanged aggregate fallback.
#!/usr/bin/env Rscript

# Paper 1 submission-grade figures.
# R is the exclusive plotting and rendering backend.
# All quantitative panels read frozen machine results from 01_reanalysis.R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ragg)
})

version_id <- "20261003_165814"
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
paper_dir <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))
result_dir <- file.path(paper_dir, "results", "versions", version_id)
base_result_dir <- file.path(paper_dir, "results", "versions", "20260809_152651")
figure_inputs <- list()
figure_dir <- file.path(result_dir, "figures")
source_dir <- file.path(result_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(name) {
  path <- file.path(result_dir, name)
  if (!file.exists(path)) path <- file.path(base_result_dir, name)
  if (!file.exists(path)) stop("Missing frozen result: ", path)
  figure_inputs[[name]] <<- data.frame(file = name, path = file.path("paper1_cmrf_pathology", "results", "versions", basename(dirname(path)), name), sha256 = digest::digest(path, file = TRUE, algo = "sha256"))
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

palette <- c(
  ink = "#263238",
  neutral = "#8B9398",
  pale = "#DDE3E7",
  blue = "#2F6B9A",
  blue_light = "#9DC3D5",
  teal = "#3C8D86",
  amber = "#D59B45",
  red = "#B85C55",
  purple = "#7566A6"
)
col <- function(name) unname(palette[name])

theme_nature <- function(base_size = 10, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = palette["ink"]),
      axis.ticks = element_line(linewidth = 0.35, colour = palette["ink"]),
      axis.text = element_text(colour = palette["ink"], size = base_size - 0.3),
      axis.title = element_text(colour = palette["ink"], size = base_size),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "#F4F6F7", colour = NA),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = palette["neutral"]),
      plot.title.position = "plot",
      plot.tag = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      legend.key.height = unit(3.2, "mm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}
theme_set(theme_nature())

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 135) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svg_path <- file.path(figure_dir, paste0(filename, ".svg"))
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  tiff_path <- file.path(figure_dir, paste0(filename, ".tiff"))
  png_path <- file.path(figure_dir, paste0(filename, ".png"))

  svglite::svglite(svg_path, width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(pdf_path, width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(tiff_path, width = w, height = h, units = "in", res = 600,
                 background = "white")
  print(plot)
  dev.off()
  ragg::agg_png(png_path, width = w, height = h, units = "in", res = 300,
                background = "white")
  print(plot)
  dev.off()
}

# Figure 1 -------------------------------------------------------------------

flow <- read_result("analysis_population_flow.csv")
center_mix <- read_result("center_fibrosis_case_mix.csv")
missing_center <- read_result("cmrf_indeterminate_by_center.csv")
transition <- read_result("legacy_consensus_transition_adults.csv")

flow_main <- flow[flow$step_order <= 6, ]
flow_main$short_step <- c(
  "Source records", "Adults", "Histology available",
  "≥1 definite CMRF", "Indeterminate component", "Complete-component"
)
flow_main$short_step <- factor(flow_main$short_step, levels = rev(flow_main$short_step))
flow_main$status <- ifelse(
  flow_main$step_order == 5, "Excluded/uncertain",
  ifelse(flow_main$step_order == 6, "Primary cohort", "Retained")
)
flow_main$label_x <- ifelse(flow_main$n > 300, flow_main$n - 18, flow_main$n + 18)
flow_main$label_hjust <- ifelse(flow_main$n > 300, 1, 0)
flow_main$label_colour <- ifelse(flow_main$n > 300, "white", col("ink"))

p1a <- ggplot(flow_main, aes(x = n, y = short_step, fill = status)) +
  geom_col(width = 0.68) +
  geom_text(
    aes(x = label_x, label = paste0("n=", n), hjust = label_hjust,
        colour = label_colour),
    size = 3.2, family = "Arial"
  ) +
  scale_colour_identity() +
  scale_fill_manual(values = c(
    "Retained" = col("blue_light"),
    "Excluded/uncertain" = col("amber"),
    "Primary cohort" = col("blue")
  )) +
  scale_x_continuous(limits = c(0, 1320), breaks = c(0, 600, 1200),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(labels = c(
    "Source records" = "Source records", "Adults" = "Adults",
    "Histology available" = "Histology\navailable",
    "≥1 definite CMRF" = "≥1 definite CMRF",
    "Indeterminate component" = "Indeterminate\ncomponent",
    "Complete-component" = "Complete\ncomponent"
  )) +
  labs(x = "Participants", y = NULL, title = "Cohort construction") +
  theme(legend.position = "none", axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

center_mix$fibrosis_stage <- factor(
  center_mix$fibrosis_stage, levels = 0:4, labels = paste0("F", 0:4)
)
p1b <- ggplot(center_mix, aes(x = center_code, y = pct, fill = fibrosis_stage)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  scale_x_discrete(labels = setNames(
    paste0(unique(center_mix$center_code), "\nn=",
           center_mix %>% distinct(center_code, N_center) %>% pull(N_center)),
    unique(center_mix$center_code)
  )) +
  scale_fill_manual(values = c(
    F0 = "#E6EEF3", F1 = "#BFD4E1", F2 = "#78ABC2",
    F3 = "#3C7A9D", F4 = "#244C68"
  )) +
  scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Analysis centre", y = "Fibrosis composition",
       fill = "Stage", title = "Centre case mix") +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank())

missing_center$component <- factor(
  missing_center$component,
  levels = c("Adiposity", "Glycaemia", "Hypertension",
             "Hypertriglyceridaemia", "Low HDL-C")
)
p1c <- ggplot(missing_center,
              aes(x = center_code, y = component, fill = indeterminate_pct)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.1f%%", indeterminate_pct)),
            size = 3.2, family = "Arial") +
  scale_fill_gradient(low = "#F4F7F9", high = palette["amber"],
                      name = "Indeterminate (%)") +
  scale_y_discrete(labels = c(
    "Adiposity" = "Adiposity", "Glycaemia" = "Glycaemia",
    "Hypertension" = "Hypertension",
    "Hypertriglyceridaemia" = "High triglycerides", "Low HDL-C" = "Low HDL-C"
  )) +
  labs(x = "Analysis centre", y = NULL,
       title = "Indeterminate components") +
  guides(fill = guide_colourbar(direction = "horizontal", title.position = "top",
                               barwidth = unit(30, "mm"), barheight = unit(3, "mm"))) +
  theme(legend.position = "bottom")

transition$legacy <- as.numeric(as.character(transition$legacy))
transition$consensus <- as.numeric(as.character(transition$consensus))
transition$direction <- ifelse(
  transition$consensus > transition$legacy, "Higher after consensus adjudication",
  ifelse(transition$consensus == transition$legacy, "Unchanged", "Lower")
)
p1d <- ggplot(transition, aes(x = legacy, y = consensus)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              linewidth = 0.4, colour = palette["neutral"]) +
  geom_point(aes(size = n, fill = direction), shape = 21,
             colour = "white", stroke = 0.4, alpha = 0.9) +
  geom_text(aes(label = n), size = 3.2, family = "Arial", colour = palette["ink"]) +
  scale_size_area(max_size = 11, name = "Participants") +
  scale_fill_manual(values = c(
    "Higher after consensus adjudication" = col("amber"),
    "Unchanged" = col("blue_light"),
    "Lower" = col("red")
  ), labels = c("Higher after consensus adjudication" = "Higher",
                "Unchanged" = "Unchanged", "Lower" = "Lower"), name = NULL) +
  scale_x_continuous(breaks = 0:5, limits = c(-0.2, 5.2)) +
  scale_y_continuous(breaks = 1:5, limits = c(0.8, 5.2)) +
  coord_equal() +
  labs(x = "Legacy CMRF count", y = "Reconstructed CMRF count",
       title = "Count reconstruction",
       subtitle = "Bubble area represents participants") +
  guides(fill = guide_legend(order = 1, nrow = 1), size = "none") +
  theme(legend.position = "bottom")

figure1 <- (p1a | p1b) / (p1c | p1d) +
  plot_layout(heights = c(0.92, 1.08)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

write.csv(flow_main, file.path(source_dir, "Figure1a_cohort_flow.csv"), row.names = FALSE)
write.csv(center_mix, file.path(source_dir, "Figure1b_center_case_mix.csv"), row.names = FALSE)
write.csv(missing_center, file.path(source_dir, "Figure1c_component_indeterminacy.csv"), row.names = FALSE)
write.csv(transition, file.path(source_dir, "Figure1d_reclassification.csv"), row.names = FALSE)
save_pub_r(figure1, "Figure1_cohort_and_adjudication", height_mm = 175)

# Figure 2 -------------------------------------------------------------------

fibdist <- read_result("fibrosis_distribution_by_cmrf.csv")
prob <- read_result("standardized_probabilities.csv")
main_cc <- read_result("main_models_complete_component.csv")
main_mi <- read_result("main_models_multiple_imputation.csv")

fibdist$fibrosis_stage <- factor(
  fibdist$fibrosis_stage, levels = 0:4, labels = paste0("F", 0:4)
)
fibdist$cmrf_count <- factor(fibdist$cmrf_count, levels = 1:5)
p2a <- ggplot(fibdist, aes(x = cmrf_count, y = pct, fill = fibrosis_stage)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = c(
    F0 = "#E6EEF3", F1 = "#BFD4E1", F2 = "#78ABC2",
    F3 = "#3C7A9D", F4 = "#244C68"
  )) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "CMRF count", y = "Participants", fill = "Fibrosis stage",
       title = "Observed fibrosis distribution") +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_blank())

make_probability_panel <- function(outcome, title, colour) {
  d <- prob[prob$outcome == outcome, ]
  ggplot(d, aes(x = cmrf_count, y = probability)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
                fill = colour, alpha = 0.16) +
    geom_line(linewidth = 0.75, colour = colour) +
    geom_point(size = 2.2, shape = 21, fill = "white", colour = colour,
               stroke = 0.75) +
    scale_x_continuous(breaks = 1:5) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, max(d$ci_high) * 1.08)) +
    labs(x = "CMRF count", y = "Standardized probability", title = title,
         subtitle = "Age, sex and centre adjusted\n1,000 stratified bootstraps") +
    theme(plot.subtitle = element_text(size = 9))
}
p2b <- make_probability_panel("f_ge2", "Significant fibrosis (F≥2)", palette["blue"])
p2c <- make_probability_panel("f_ge3", "Advanced fibrosis (F≥3)", palette["teal"])

main_compare <- bind_rows(
  main_cc %>% mutate(dataset = "Complete-component"),
  main_mi %>% mutate(dataset = "Multiple imputation")
)
main_compare$outcome_label <- factor(
  main_compare$outcome, levels = c("f_ge3", "f_ge2"),
  labels = c("Advanced fibrosis (F≥3)", "Significant fibrosis (F≥2)")
)
main_compare$dataset <- factor(
  main_compare$dataset, levels = c("Complete-component", "Multiple imputation")
)
p2d <- ggplot(main_compare,
              aes(x = OR, y = outcome_label, colour = dataset, shape = dataset)) +
  geom_vline(xintercept = 1, linewidth = 0.4, linetype = 2,
             colour = palette["neutral"]) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0, linewidth = 0.65,
                position = position_dodge(width = 0.42)) +
  geom_point(size = 2.4, stroke = 0.7,
             position = position_dodge(width = 0.42)) +
  scale_colour_manual(values = c(
    "Complete-component" = col("blue"),
    "Multiple imputation" = col("amber")
  )) +
  scale_shape_manual(values = c("Complete-component" = 16, "Multiple imputation" = 17)) +
  scale_y_discrete(labels = c("Advanced fibrosis (F≥3)" = "F≥3",
                             "Significant fibrosis (F≥2)" = "F≥2")) +
  scale_x_log10(breaks = c(1, 1.2, 1.4, 1.6, 1.8),
                limits = c(0.98, 1.85)) +
  labs(x = "Odds ratio per CMRF\n(log scale)", y = NULL,
       colour = NULL, shape = NULL, title = "Complete-component vs MI") +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  theme(legend.position = "top")

figure2 <- (p2a | p2b) / (p2c | p2d) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

write.csv(fibdist, file.path(source_dir, "Figure2a_fibrosis_distribution.csv"), row.names = FALSE)
write.csv(prob, file.path(source_dir, "Figure2bc_standardized_probabilities.csv"), row.names = FALSE)
write.csv(main_compare, file.path(source_dir, "Figure2d_complete_and_mi.csv"), row.names = FALSE)
save_pub_r(figure2, "Figure2_cmrf_burden_and_fibrosis", height_mm = 175)

# Figure 3 -------------------------------------------------------------------

comp_cc <- read_result("component_models_complete_component.csv")
comp_mi <- read_result("component_models_multiple_imputation.csv")
gly_cc <- read_result("glycaemia_decomposition_complete_component.csv")
gly_mi <- read_result("glycaemia_decomposition_multiple_imputation.csv")

component_labels <- c(
  cmrf1 = "Adiposity", cmrf2 = "Glycaemia", cmrf3 = "Hypertension",
  cmrf4 = "Hypertriglyceridaemia", cmrf5 = "Low HDL-C"
)
comp <- bind_rows(
  comp_cc %>% mutate(dataset = "Complete-component"),
  comp_mi %>% mutate(dataset = "Multiple imputation")
)
comp$component <- factor(
  unname(component_labels[comp$term]),
  levels = rev(unname(component_labels))
)
comp$outcome_label <- factor(
  comp$outcome, levels = c("f_ge2", "f_ge3"),
  labels = c("Significant fibrosis (F≥2)", "Advanced fibrosis (F≥3)")
)
comp$dataset <- factor(comp$dataset, levels = c("Complete-component", "Multiple imputation"))
comp$supported <- comp$p_fdr_all_10 < 0.05

p3a <- ggplot(comp, aes(x = OR, y = component, colour = dataset, shape = supported)) +
  geom_vline(xintercept = 1, linewidth = 0.4, linetype = 2,
             colour = palette["neutral"]) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0, linewidth = 0.55,
                position = position_dodge(width = 0.48)) +
  geom_point(size = 2.1, stroke = 0.7,
             position = position_dodge(width = 0.48)) +
  facet_wrap(~outcome_label, ncol = 2, labeller = as_labeller(c(
    "Significant fibrosis (F≥2)" = "Significant fibrosis\n(F≥2)",
    "Advanced fibrosis (F≥3)" = "Advanced fibrosis\n(F≥3)"
  ))) +
  scale_colour_manual(values = c(
    "Complete-component" = col("blue"),
    "Multiple imputation" = col("amber")
  )) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c(`FALSE` = "BH-FDR ≥0.05", `TRUE` = "BH-FDR <0.05")) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 4),
                labels = c("0.1", "0.5", "1", "2", "4"), limits = c(0.06, 5.2)) +
  scale_y_discrete(labels = c(
    "Adiposity" = "Adiposity", "Glycaemia" = "Glycaemia",
    "Hypertension" = "Hypertension",
    "Hypertriglyceridaemia" = "High triglycerides", "Low HDL-C" = "Low HDL-C"
  )) +
  labs(x = "Mutually adjusted odds ratio (log scale)", y = NULL,
       colour = NULL, shape = NULL,
       title = "Component associations",
       subtitle = "Complete-component and MI analyses") +
  theme(legend.position = "top", legend.box = "vertical")

gly <- bind_rows(
  gly_cc %>% mutate(dataset = "Complete-component"),
  gly_mi %>% mutate(dataset = "Multiple imputation")
)
gly$level <- factor(
  gly$level,
  levels = c("Recorded T2DM", "Dysglycaemia without recorded T2DM")
)
gly$outcome_label <- factor(
  gly$outcome, levels = c("f_ge2", "f_ge3"),
  labels = c("F≥2", "F≥3")
)
gly$dataset <- factor(gly$dataset, levels = c("Complete-component", "Multiple imputation"))
gly$supported <- gly$p_holm_all_4 < 0.05

p3b <- ggplot(gly, aes(x = OR, y = level, colour = outcome_label,
                       shape = dataset, alpha = supported)) +
  geom_vline(xintercept = 1, linewidth = 0.4, linetype = 2,
             colour = palette["neutral"]) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0, linewidth = 0.65,
                position = position_dodge(width = 0.54)) +
  geom_point(size = 2.5, stroke = 0.7,
             position = position_dodge(width = 0.54)) +
  scale_colour_manual(values = c("F≥2" = col("blue"), "F≥3" = col("teal"))) +
  scale_shape_manual(values = c("Complete-component" = 16, "Multiple imputation" = 17)) +
  scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 1), guide = "none") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4, 6), limits = c(0.85, 7.1)) +
  scale_y_discrete(labels = c(
    "Recorded T2DM" = "Recorded T2DM",
    "Dysglycaemia without recorded T2DM" = "Dysglycaemia without\nrecorded T2DM"
  )) +
  labs(
    x = "Odds ratio vs normoglycaemia (log scale)", y = NULL,
    colour = "Outcome", shape = "Dataset",
    title = "Glycaemia decomposition",
    subtitle = "Adjusted for other CMRF components, age, sex and centre"
  ) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, nrow = 1)) +
  theme(legend.position = "top", legend.box = "vertical")

figure3 <- p3a / p3b +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

write.csv(comp, file.path(source_dir, "Figure3a_component_models.csv"), row.names = FALSE)
write.csv(gly, file.path(source_dir, "Figure3b_glycaemia_decomposition.csv"), row.names = FALSE)
save_pub_r(figure3, "Figure3_components_and_glycaemia", height_mm = 190)

# Figure 4 -------------------------------------------------------------------

robust <- read_result("robustness_models.csv")
robust$outcome_label <- factor(
  robust$outcome, levels = c("f_ge2", "f_ge3"),
  labels = c("Significant fibrosis (F≥2)", "Advanced fibrosis (F≥3)")
)
scenario_levels <- robust %>%
  distinct(analysis, scenario_order) %>%
  arrange(scenario_order) %>%
  pull(analysis)
robust$analysis <- factor(robust$analysis, levels = rev(scenario_levels))
robust$class <- case_when(
  grepl("^Leave out", robust$analysis) ~ "Leave-one-centre-out",
  robust$analysis %in% c("Adult primary", "Adult raw-variable MI") ~ "Primary/MI",
  TRUE ~ "Other sensitivity"
)

p4 <- ggplot(robust, aes(x = OR, y = analysis, colour = class)) +
  geom_vline(xintercept = 1, linewidth = 0.45, linetype = 2,
             colour = palette["neutral"]) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 2.0) +
  facet_wrap(~outcome_label, ncol = 2) +
  scale_colour_manual(values = c(
    "Primary/MI" = col("blue"),
    "Other sensitivity" = col("teal"),
    "Leave-one-centre-out" = col("amber")
  )) +
  scale_x_log10(breaks = c(0.8, 1, 1.2, 1.5, 2.0), limits = c(0.78, 2.25)) +
  labs(x = "Odds ratio per additional CMRF (log scale)", y = NULL,
       colour = NULL, title = "Robustness of the CMRF burden association",
       subtitle = "Points and lines show odds ratios and 95% confidence intervals") +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 9)
  )

write.csv(robust, file.path(source_dir, "Figure4_robustness_models.csv"), row.names = FALSE)
save_pub_r(p4, "Figure4_robustness", height_mm = 170)

# Supplementary Figure S1: MI diagnostics -----------------------------------

imp <- readRDS(file.path(base_result_dir, "multiple_imputation_object.rds"))
chain_array <- imp$chainMean
trace_rows <- list()
if (length(dim(chain_array)) == 3) {
  vars <- intersect(c("wc", "glu", "hba1c", "tg", "hdl"), dimnames(chain_array)[[1]])
  for (v in vars) {
    vi <- match(v, dimnames(chain_array)[[1]])
    for (chain in seq_len(dim(chain_array)[3])) {
      trace_rows[[length(trace_rows) + 1]] <- data.frame(
        variable = v,
        iteration = seq_len(dim(chain_array)[2]),
        chain = chain,
        chain_mean = chain_array[vi, , chain],
        stringsAsFactors = FALSE
      )
    }
  }
}
trace <- do.call(rbind, trace_rows)
mi_missing <- read_result("multiple_imputation_missingness.csv")
mi_missing <- mi_missing[mi_missing$method != "", ]
mi_distributions <- read_result("multiple_imputation_observed_imputed_distributions.csv")
mi_precision <- read_result("multiple_imputation_precision_diagnostics.csv")

variable_labels <- c(
  wc = "Waist circumference",
  glu = "Glucose",
  hba1c = "HbA1c",
  tg = "Triglycerides",
  hdl = "HDL-C"
)
trace$variable <- factor(
  trace$variable, levels = names(variable_labels), labels = variable_labels
)
mi_distributions$variable <- factor(
  mi_distributions$variable,
  levels = names(variable_labels), labels = variable_labels
)
mi_distributions$source <- factor(
  mi_distributions$source, levels = c("Observed", "Imputed")
)

pS1a <- ggplot(trace, aes(x = iteration, y = chain_mean, group = chain)) +
  geom_line(alpha = 0.28, linewidth = 0.25, colour = palette["blue"]) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  labs(x = "Iteration", y = "Chain mean", title = "Multiple-imputation trace") +
  theme(strip.text = element_text(face = "bold"))

pS1b <- ggplot(
  mi_distributions,
  aes(x = value, colour = source, fill = source)
) +
  geom_density(alpha = 0.12, linewidth = 0.55, adjust = 1.1) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  scale_colour_manual(values = c("Observed" = col("blue"), "Imputed" = col("amber"))) +
  scale_fill_manual(values = c("Observed" = col("blue"), "Imputed" = col("amber"))) +
  labs(
    x = "Observed or imputed value", y = "Density",
    colour = NULL, fill = NULL,
    title = "Observed and imputed distributions"
  ) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

pS1c <- ggplot(mi_missing, aes(x = missing_pct, y = reorder(variable, missing_pct))) +
  geom_col(width = 0.65, fill = palette["amber"]) +
  geom_text(aes(label = sprintf("%.1f%%", missing_pct)), hjust = -0.12,
            size = 3.2, family = "Arial") +
  scale_x_continuous(limits = c(0, max(mi_missing$missing_pct) * 1.18),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Missing before imputation", y = NULL, title = "Variables imputed") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

mi_precision_long <- rbind(
  data.frame(
    outcome = mi_precision$outcome,
    metric = "Fraction of missing information",
    value = 100 * mi_precision$fraction_missing_information_approx
  ),
  data.frame(
    outcome = mi_precision$outcome,
    metric = "Monte Carlo error / pooled SE",
    value = mi_precision$monte_carlo_error_pct_pooled_se
  )
)
mi_precision_long$outcome <- factor(
  mi_precision_long$outcome,
  levels = c("f_ge2", "f_ge3"), labels = c("F≥2", "F≥3")
)
mi_precision_long$metric <- factor(
  mi_precision_long$metric,
  levels = c("Fraction of missing information", "Monte Carlo error / pooled SE")
)
pS1d <- ggplot(
  mi_precision_long,
  aes(x = value, y = metric, fill = outcome)
) +
  geom_col(width = 0.62, position = position_dodge(width = 0.68)) +
  geom_text(
    aes(label = sprintf("%.1f%%", value)),
    hjust = -0.12, size = 3.2, family = "Arial",
    position = position_dodge(width = 0.68)
  ) +
  scale_fill_manual(values = c("F≥2" = col("blue"), "F≥3" = col("teal"))) +
  scale_x_continuous(
    limits = c(0, max(mi_precision_long$value) * 1.22),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Percentage", y = NULL, fill = "Outcome",
    title = "Imputation precision"
  ) +
  theme(
    legend.position = "top",
    axis.line.y = element_blank(), axis.ticks.y = element_blank()
  )

figS1 <- pS1a / pS1b / (pS1c | pS1d) +
  plot_layout(heights = c(1.05, 1.05, 0.75)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

write.csv(trace, file.path(source_dir, "FigureS1a_mi_trace.csv"), row.names = FALSE)
write.csv(
  mi_distributions,
  file.path(source_dir, "FigureS1b_observed_imputed_distributions.csv"),
  row.names = FALSE
)
write.csv(mi_missing, file.path(source_dir, "FigureS1c_mi_missingness.csv"), row.names = FALSE)
write.csv(mi_precision, file.path(source_dir, "FigureS1d_mi_precision.csv"), row.names = FALSE)
save_pub_r(figS1, "FigureS1_multiple_imputation_diagnostics", height_mm = 220)

# Figure legends and QA notes ------------------------------------------------

legend_text <- c(
  "Fig. 1 | Cohort construction and CMRF adjudication. a Sequential construction of the adult complete-component cohort. b Fibrosis-stage composition across three analysis centres (HZ, JL and RJ; original ZS included in RJ). c Percentage of participants in the adult eligible cohort whose component status was indeterminate, by centre. d Transition from the legacy to reconstructed CMRF count in the adult complete-component cohort; point area is proportional to the number of participants. Source data are provided as a Source Data file.",
  "Fig. 2 | CMRF burden and biopsy-defined fibrosis. a Observed F0–F4 composition by CMRF count in 935 adults. b,c Marginal standardized probabilities of significant (F≥2; 550 events) and advanced fibrosis (F≥3; 246 events) from logistic models adjusted for age, sex and centre; bands show centre-stratified percentile 95% confidence intervals from 1,000 bootstrap samples. d Adjusted odds ratios per additional CMRF in the complete-component cohort (n=935) and raw-variable multiple-imputation cohort (n=1,085; 50 imputations). Source data are provided as a Source Data file.",
  "Fig. 3 | Component-level and glycaemic associations with fibrosis. a Mutually adjusted component odds ratios for F≥2 and F≥3 in complete-component and raw-variable multiple-imputation analyses. Filled points denote associations passing Benjamini–Hochberg correction across ten component–outcome tests. b Associations of dysglycaemia without recorded T2DM and recorded T2DM versus normoglycaemia, adjusted for age, sex, centre and the other four CMRF components; opacity denotes Holm-adjusted support across four comparisons. Lines indicate 95% confidence intervals. Source data are provided as a Source Data file.",
  "Fig. 4 | Robustness of the CMRF burden association. Adjusted odds ratios per additional CMRF across the adult primary analysis, raw-variable multiple imputation, alternative eligibility and exposure definitions, covariate functional forms, extreme missingness scenarios and leave-one-centre-out analyses. Estimates are shown separately for F≥2 and F≥3; lines indicate 95% confidence intervals. These analyses address different sensitivity estimands and should not be interpreted as interchangeable effect estimates. Source data are provided as a Source Data file.",
  "Supplementary Fig. 1 | Multiple-imputation diagnostics. a Chain means across 20 iterations for 50 imputations of the continuous raw measurements used to rederive CMRF status. b Observed and imputed distributions of the continuous raw measurements. c Proportion missing before imputation for continuous measurements and binary auxiliary variables. d Approximate fraction of missing information and Monte Carlo error as a percentage of the pooled standard error for the CMRF-count coefficient. Predictive mean matching was used for continuous variables and logistic regression for binary auxiliary variables."
)
writeLines(legend_text, file.path(result_dir, "figure_legends.txt"))

qa <- data.frame(
  figure = c("Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure S1"),
  core_conclusion = c(
    "Population construction and CMRF adjudication are traceable",
    "Higher CMRF burden is associated with higher adjusted fibrosis probability",
    "Glycaemia, particularly recorded T2DM, shows the strongest component signal",
    "The burden association is directionally robust across planned sensitivities",
    "Raw-variable MI shows stable traces, plausible distributions and low simulation error"
  ),
  n_definition = c(
    "Adult source/eligible/complete-component cohorts",
    "Adult complete-component n=935; MI n=1085",
    "Adult complete-component n=935; MI n=1085",
    "Scenario-specific N shown in source data",
    "Adult eligible n=1085; m=50"
  ),
  interval = c("Not applicable", "95% bootstrap CI and model-based 95% CI",
               "Model-based 95% CI", "Model-based 95% CI",
               "Trace, distribution, FMI and Monte Carlo error diagnostics"),
  adjustment = c("Descriptive", "Age, sex, centre", "Age, sex, centre and other components",
                 "Scenario-specific; source data", "m=50, maxit=20"),
  exports = "SVG, PDF, 600 dpi TIFF, 300 dpi PNG",
  stringsAsFactors = FALSE
)
write.csv(qa, file.path(result_dir, "figure_QA_contract.csv"), row.names = FALSE)

message("Completed Paper 1 figures: ", figure_dir)

write.csv(do.call(rbind, figure_inputs), file.path(result_dir, "figure_input_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo_figures.txt"))

# Rendering integrity: every frozen CI remains a finite error-bar segment.
ci_plots <- list(Figure2d = p2d, Figure3a = p3a, Figure3b = p3b, Figure4 = p4)
ci_sources <- list(Figure2d = main_compare, Figure3a = comp, Figure3b = gly, Figure4 = robust)
ci_checks <- do.call(rbind, lapply(names(ci_plots), function(panel) {
  rendered <- ggplot_build(ci_plots[[panel]])$data[[2]]
  input <- ci_sources[[panel]]
  finite <- is.finite(rendered$xmin) & is.finite(rendered$xmax)
  stopifnot(nrow(rendered) == nrow(input), all(finite))
  data.frame(panel = panel, source_intervals = nrow(input),
    rendered_finite_intervals = sum(finite), minimum_source_CI = min(input$ci_low),
    maximum_source_CI = max(input$ci_high), pass = TRUE, author = "mjwang",
    stringsAsFactors = FALSE)
}))
write.csv(ci_checks, file.path(result_dir, "figure_render_checks.csv"), row.names = FALSE)
