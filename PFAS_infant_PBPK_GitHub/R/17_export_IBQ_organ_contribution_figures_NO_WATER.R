#!/usr/bin/env Rscript
# ================================================================
# NO-WATER manuscript figure export
#
# Figure 3: Six-PFAS contributions to within-organ CTWA composition
# Figure 4: Individual serum BMDL10 ratio distribution
#
# Design goals:
#   - match the earlier journal-style figures
#   - no plot titles
#   - smaller, cleaner text
#   - clearer annotations
#   - consistent colour style across figures
#
# Input:
#   Frozen S1 PBPK workbook only
#   No ODE rerun
#   No drinking-water data
# ================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(tibble)
  library(scales)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", getwd()), winslash = "/", mustWork = FALSE)
VERSION <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
MODE <- tolower(Sys.getenv("PFAS_TRAJECTORY_MODE", "proxy"))
if (MODE != "proxy") VERSION <- paste0(VERSION, "_beta_", MODE)
default_input <- Sys.getenv("PFAS_PBPK_RESULTS", unset = file.path(
  ROOT, "results", VERSION, "tables",
  "PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx"
))

input_xlsx <- if (length(args) >= 1 && nzchar(args[[1]])) {
  args[[1]]
} else {
  default_input
}

if (!file.exists(input_xlsx)) {
  stop(
    "Frozen S1 PBPK workbook not found:\n", input_xlsx,
    "\nRun from the project root or pass the full workbook path."
  )
}

output_dir <- if (length(args) >= 2 && nzchar(args[[2]])) {
  args[[2]]
} else {
  Sys.getenv("PFAS_17_OUTPUT", unset = file.path(
    getwd(), "results", "PFAS_infant_PBPK_S1_renal_totalplasma_v2",
    "IBQ_organ_contribution_figures_NO_WATER"
  ))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("[17-v2] Reading frozen S1 PBPK workbook: ", normalizePath(input_xlsx))

# -----------------------------------------------------------------
# 1. Style settings
# -----------------------------------------------------------------
FONT_FAMILY <- "sans"

compound_order <- c("PFOA", "PFOS", "PFHxS", "PFNA", "PFDA", "PFUdA")
four_compounds <- c("PFOA", "PFNA", "PFHxS", "PFOS")
organ_order <- c("Plasma", "Liver", "Kidney")
scenario_order <- c("Lower", "Central", "Upper")

compound_colours <- c(
  PFOA  = "#0072B2",
  PFOS  = "#D55E00",
  PFHxS = "#56B4E9",
  PFNA  = "#E69F00",
  PFDA  = "#009E73",
  PFUdA = "#CC79A7"
)

scenario_colours <- c(
  Lower   = "#0072B2",
  Central = "#009E73",
  Upper   = "#D55E00"
)

theme_journal <- function(base_size = 11.5) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      text = element_text(colour = "#20252B", family = FONT_FAMILY),
      axis.title = element_text(size = 12, colour = "#20252B"),
      axis.text = element_text(size = 10, colour = "#20252B"),
      axis.line = element_line(colour = "#30363D", linewidth = 0.45),
      axis.ticks = element_line(colour = "#30363D", linewidth = 0.4),
      axis.ticks.length = unit(2.0, "mm"),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(
        size = 11,
        face = "bold",
        colour = "#20252B",
        margin = margin(b = 3, t = 2)
      ),
      legend.position = "bottom",
      legend.title = element_text(size = 10.5, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.width = unit(7, "mm"),
      legend.key.height = unit(4, "mm"),
      panel.grid.major.y = element_line(colour = "#E8EDF2", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(t = 6, r = 8, b = 6, l = 8)
    )
}

save_journal_figure <- function(plot_object, file_stem, width_in, height_in, folder) {
  plot_object <- plot_object + labs(title = NULL, subtitle = NULL, caption = NULL)
  
  pdf_file  <- file.path(folder, paste0(file_stem, ".pdf"))
  tiff_file <- file.path(folder, paste0(file_stem, ".tiff"))
  png_file  <- file.path(folder, paste0(file_stem, ".png"))
  
  if (capabilities("cairo")) {
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot_object,
      device = grDevices::cairo_pdf,
      width = width_in, height = height_in, units = "in",
      bg = "white", family = FONT_FAMILY, limitsize = FALSE
    )
  } else {
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot_object,
      device = "pdf",
      width = width_in, height = height_in, units = "in",
      bg = "white", family = FONT_FAMILY, useDingbats = FALSE,
      limitsize = FALSE
    )
  }
  
  ggplot2::ggsave(
    filename = tiff_file,
    plot = plot_object,
    device = "tiff",
    width = width_in, height = height_in, units = "in",
    dpi = 900, compression = "lzw", bg = "white", limitsize = FALSE
  )
  
  ggplot2::ggsave(
    filename = png_file,
    plot = plot_object,
    device = "png",
    width = width_in, height = height_in, units = "in",
    dpi = 900, bg = "white", limitsize = FALSE
  )
  
  message("[Saved] ", file_stem, " | PDF + 900-dpi TIFF + 900-dpi PNG")
}

# -----------------------------------------------------------------
# 2. Read and validate frozen PBPK output
# -----------------------------------------------------------------
dat0 <- read_excel(input_xlsx, sheet = "PBPK_Results")

required_columns <- c(
  "ID", "Time_days", "Scenario", "Compound",
  "C_Plasma", "C_Liver", "C_Kidney"
)
missing_columns <- setdiff(required_columns, names(dat0))
if (length(missing_columns) > 0) {
  stop("Missing required columns in PBPK_Results: ",
       paste(missing_columns, collapse = ", "))
}

dat <- dat0 %>%
  mutate(
    ID = as.character(ID),
    Time_days = as.numeric(Time_days),
    Scenario = trimws(as.character(Scenario)),
    Compound = trimws(as.character(Compound))
  ) %>%
  filter(toupper(Scenario) == "S1")

target_compounds <- c("PFOA", "PFNA", "PFDA", "PFUdA", "PFHxS", "PFOS")
expected_n <- 175L * 6L * 54L

if (dplyr::n_distinct(dat$ID) != 175L) {
  stop("Expected 175 unique infants; found ", dplyr::n_distinct(dat$ID), ".")
}
if (!setequal(unique(dat$Compound), target_compounds)) {
  stop("Unexpected compounds in S1 output.")
}
if (nrow(dat) != expected_n) {
  stop("Unexpected S1 record count: ", nrow(dat), "; expected ", expected_n, ".")
}
if (anyDuplicated(dat[c("ID", "Compound", "Time_days")])) {
  stop("Duplicate ID-compound-time records were detected.")
}
if (anyNA(dat[c("ID", "Time_days", "Compound", "C_Plasma", "C_Liver", "C_Kidney")])) {
  stop("Missing values detected in required concentration fields.")
}

message("[17-v2] Validated 175 infants x 6 compounds x 54 time points.")

# -----------------------------------------------------------------
# 3. CTWA calculations and within-organ contributions
# -----------------------------------------------------------------
calc_ctwa <- function(x, value_col) {
  x <- x[order(x$Time_days), , drop = FALSE]
  tt <- x$Time_days
  yy <- x[[value_col]]
  
  if (nrow(x) != 54L ||
      abs(tt[1] - 0) > 1e-8 ||
      abs(tail(tt, 1) - 365.25) > 1e-8 ||
      !any(abs(tt - 364) < 1e-8)) {
    stop("Unexpected time grid for an ID-compound trajectory.")
  }
  if (anyNA(yy)) stop("Missing concentration values in ", value_col, ".")
  
  c365 <- approx(x = tt, y = yy, xout = 365, method = "linear")$y
  keep <- tt < 365
  t_int <- c(tt[keep], 365)
  y_int <- c(yy[keep], c365)
  
  auc <- sum(diff(t_int) * (head(y_int, -1) + tail(y_int, -1)) / 2)
  auc / 365
}

ctwa_wide <- dat %>%
  group_by(ID, Compound) %>%
  group_modify(~ tibble(
    CTWA_Plasma_ng_L = calc_ctwa(.x, "C_Plasma"),
    CTWA_Liver_ng_L  = calc_ctwa(.x, "C_Liver"),
    CTWA_Kidney_ng_L = calc_ctwa(.x, "C_Kidney")
  )) %>%
  ungroup()

ctwa_long <- ctwa_wide %>%
  pivot_longer(
    cols = starts_with("CTWA_"),
    names_to = "Organ_variable",
    values_to = "CTWA_0_365_ng_L"
  ) %>%
  mutate(
    Organ = recode(
      Organ_variable,
      CTWA_Plasma_ng_L = "Plasma",
      CTWA_Liver_ng_L = "Liver",
      CTWA_Kidney_ng_L = "Kidney"
    ),
    Organ = factor(Organ, levels = organ_order),
    Compound = factor(Compound, levels = compound_order)
  ) %>%
  group_by(ID, Organ) %>%
  mutate(
    Six_PFAS_sum_CTWA_ng_L = sum(CTWA_0_365_ng_L),
    Contribution_pct = 100 * CTWA_0_365_ng_L / Six_PFAS_sum_CTWA_ng_L
  ) %>%
  ungroup()

share_qc <- ctwa_long %>%
  group_by(ID, Organ) %>%
  summarise(Share_sum_pct = sum(Contribution_pct), .groups = "drop")

if (any(abs(share_qc$Share_sum_pct - 100) > 1e-8)) {
  stop("Within-organ contribution percentages failed the 100% compositional QC.")
}

organ_summary <- ctwa_long %>%
  group_by(Organ, Compound) %>%
  summarise(
    N = n_distinct(ID),
    Median_pct = median(Contribution_pct),
    Q25_pct = unname(quantile(Contribution_pct, 0.25)),
    Q75_pct = unname(quantile(Contribution_pct, 0.75)),
    P90_pct = unname(quantile(Contribution_pct, 0.90)),
    .groups = "drop"
  )

pfoa_pfos_combined <- ctwa_long %>%
  filter(Compound %in% c("PFOA", "PFOS")) %>%
  group_by(ID, Organ) %>%
  summarise(Combined_share_pct = sum(Contribution_pct), .groups = "drop") %>%
  group_by(Organ) %>%
  summarise(
    N = n(),
    Median_pct = median(Combined_share_pct),
    Q25_pct = unname(quantile(Combined_share_pct, 0.25)),
    Q75_pct = unname(quantile(Combined_share_pct, 0.75)),
    P90_pct = unname(quantile(Combined_share_pct, 0.90)),
    .groups = "drop"
  )

# -----------------------------------------------------------------
# 4. Day-365 serum BMDL10 ratio calculations
# -----------------------------------------------------------------
ratio_table <- tribble(
  ~Scenario, ~Compound, ~R_serum_plasma,
  "Lower",   "PFOA", 1.1,
  "Lower",   "PFNA", 1.2,
  "Lower",   "PFHxS", 1.0,
  "Lower",   "PFOS", 1.1,
  "Central", "PFOA", 1.2,
  "Central", "PFNA", 1.3,
  "Central", "PFHxS", 1.1,
  "Central", "PFOS", 1.1,
  "Upper",   "PFOA", 1.2,
  "Upper",   "PFNA", 1.4,
  "Upper",   "PFHxS", 1.1,
  "Upper",   "PFOS", 1.2
)

day365_plasma <- dat %>%
  filter(Compound %in% four_compounds) %>%
  group_by(ID, Compound) %>%
  group_modify(~ {
    x <- arrange(.x, Time_days)
    c365_ng_L <- approx(
      x = x$Time_days,
      y = x$C_Plasma,
      xout = 365,
      method = "linear"
    )$y
    tibble(C_plasma_day365_ng_mL = c365_ng_L / 1000)
  }) %>%
  ungroup()

ibq_compound <- day365_plasma %>%
  inner_join(
    ratio_table,
    by = "Compound",
    relationship = "many-to-many"
  ) %>%
  mutate(
    Scenario = factor(Scenario, levels = scenario_order),
    C_serum_equivalent_ng_mL = C_plasma_day365_ng_mL * R_serum_plasma
  )

stopifnot(nrow(ibq_compound) == 175 * 4 * 3)
stopifnot(all(ibq_compound %>% count(ID, Scenario) %>% pull(n) == 4))

ibq_individual <- ibq_compound %>%
  group_by(ID, Scenario) %>%
  summarise(
    SerumEq_Sum4_ng_mL = sum(C_serum_equivalent_ng_mL),
    R_BMDL10 = SerumEq_Sum4_ng_mL / 17.5,
    .groups = "drop"
  )

ibq_summary <- ibq_individual %>%
  group_by(Scenario) %>%
  summarise(
    N = n_distinct(ID),
    SerumEq_Median_ng_mL = median(SerumEq_Sum4_ng_mL),
    SerumEq_Q25_ng_mL = unname(quantile(SerumEq_Sum4_ng_mL, 0.25)),
    SerumEq_Q75_ng_mL = unname(quantile(SerumEq_Sum4_ng_mL, 0.75)),
    R_BMDL10_Median = median(R_BMDL10),
    R_BMDL10_Q25 = unname(quantile(R_BMDL10, 0.25)),
    R_BMDL10_Q75 = unname(quantile(R_BMDL10, 0.75)),
    R_BMDL10_P90 = unname(quantile(R_BMDL10, 0.90)),
    R_BMDL10_P95 = unname(quantile(R_BMDL10, 0.95)),
    N_R_BMDL10_ge_1 = sum(R_BMDL10 >= 1),
    Pct_R_BMDL10_ge_1 = mean(R_BMDL10 >= 1) * 100,
    .groups = "drop"
  )

central_component <- ibq_compound %>%
  filter(Scenario == "Central") %>%
  group_by(ID) %>%
  mutate(
    Share_pct = 100 * C_serum_equivalent_ng_mL / sum(C_serum_equivalent_ng_mL)
  ) %>%
  ungroup() %>%
  group_by(Compound) %>%
  summarise(
    N = n_distinct(ID),
    Median_pct = median(Share_pct),
    Q25_pct = unname(quantile(Share_pct, 0.25)),
    Q75_pct = unname(quantile(Share_pct, 0.75)),
    .groups = "drop"
  )

# -----------------------------------------------------------------
# 5. Figure 3. Within-organ composition of six PFAS
# -----------------------------------------------------------------
p3 <- ggplot(
  ctwa_long,
  aes(x = Compound, y = Contribution_pct, fill = Compound)
) +
  geom_boxplot(
    width = 0.62,
    linewidth = 0.42,
    colour = "#30363D",
    outlier.shape = NA,
    alpha = 0.95
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.7,
    fill = "white",
    colour = "#20252B",
    stroke = 0.35
  ) +
  facet_wrap(~ Organ, nrow = 1) +
  scale_fill_manual(values = compound_colours, guide = "none") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = NULL,
    y = expression("Contribution to six-PFAS " * C[TWA] * " (%)")
  ) +
  theme_journal() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(
      size = 9.5,
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 9.8),
    axis.title.y = element_text(size = 12),
    panel.spacing = unit(6.5, "mm"),
    legend.position = "none"
  )

save_journal_figure(
  p3,
  file_stem = "Figure_3_organ_contribution_boxplots",
  width_in = 7.2,
  height_in = 4.6,
  folder = output_dir
)

# -----------------------------------------------------------------
# 6. Figure 4. Serum BMDL10 ratio distributions across conversion scenarios
# -----------------------------------------------------------------
ibq_labels <- ibq_summary %>%
  mutate(
    Scenario = factor(Scenario, levels = scenario_order),
    Label = sprintf("%d (%.1f%%) >= 1", N_R_BMDL10_ge_1, Pct_R_BMDL10_ge_1)
  )

ymax_plot <- max(ibq_individual$R_BMDL10, na.rm = TRUE)
label_y <- ymax_plot * 1.04
upper_lim <- ymax_plot * 1.12

p4 <- ggplot(
  ibq_individual,
  aes(x = Scenario, y = R_BMDL10, fill = Scenario)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "#5B6770"
  ) +
  geom_boxplot(
    width = 0.58,
    linewidth = 0.42,
    colour = "#30363D",
    outlier.shape = 21,
    outlier.size = 1.2,
    outlier.stroke = 0.25,
    outlier.fill = "white"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.9,
    fill = "white",
    colour = "#20252B",
    stroke = 0.35
  ) +
  geom_text(
    data = ibq_labels,
    aes(x = Scenario, y = label_y, label = Label),
    inherit.aes = FALSE,
    size = 3.0,
    lineheight = 0.95
  ) +
  scale_fill_manual(values = scenario_colours, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(0, upper_lim),
    expand = expansion(mult = c(0, 0)),
    labels = label_number(accuracy = 0.1)
  ) +
  labs(
    x = "Serum-to-plasma ratio scenario",
    y = expression("Serum benchmark ratio (" * R[BMDL10] * ")")
  ) +
  theme_journal() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9.8),
    axis.title = element_text(size = 12)
  )

save_journal_figure(
  p4,
  file_stem = "Figure_4_R_BMDL10_distribution",
  width_in = 5.6,
  height_in = 4.8,
  folder = output_dir
)

# -----------------------------------------------------------------
# 7. Export numerical results and QC
# -----------------------------------------------------------------
ratio_display <- ratio_table %>%
  group_by(Scenario) %>%
  summarise(
    Ratio_PFOA_PFNA_PFHxS_PFOS = paste(
      R_serum_plasma[match(c("PFOA", "PFNA", "PFHxS", "PFOS"), Compound)],
      collapse = " / "
    ),
    .groups = "drop"
  )

ibq_summary_export <- ibq_summary %>%
  mutate(Scenario = as.character(Scenario)) %>%
  left_join(ratio_display, by = "Scenario")

qc <- tibble(
  QC_item = c(
    "Source workbook",
    "Scenario retained",
    "Unique infants",
    "Compounds",
    "Records",
    "Expected records",
    "Time points per ID-compound",
    "Duplicate ID-compound-time rows",
    "AUC interval",
    "Day-365 handling",
    "CTWA definition",
    "Contribution denominator",
    "Compositional QC",
    "R_BMDL10 join structure",
    "R_BMDL10 expansion rows",
    "EFSA reference point",
    "PBPK ODE rerun"
  ),
  Result = c(
    basename(input_xlsx),
    "S1 only",
    dplyr::n_distinct(dat$ID),
    paste(target_compounds, collapse = ", "),
    nrow(dat),
    expected_n,
    54,
    sum(duplicated(dat[c("ID", "Compound", "Time_days")])),
    "0-365 days",
    "Linear interpolation between days 364 and 365.25",
    "Trapezoidal AUC0-365 divided by 365 days",
    "Sum of six compound CTWA values in the same infant and organ",
    "Each infant-organ set sums to 100%",
    "Intentional many-to-many by compound x 3 scenarios",
    nrow(ibq_compound),
    "17.5 ng/mL serum BMDL10 for the four-PFAS mixture",
    "No"
  )
)

write.xlsx(
  list(
    Organ_Contribution_Summary = organ_summary,
    PFOA_PFOS_Combined = pfoa_pfos_combined,
    R_BMDL10_Summary = ibq_summary_export,
    R_BMDL10_Individual = ibq_individual,
    SerumEq_Component_Shares = central_component,
    Organ_Contribution_Individual = ctwa_long,
    QC = qc
  ),
  file = file.path(output_dir, "PFAS_R_BMDL10_organ_contribution_results.xlsx"),
  overwrite = TRUE
)

message("[17-v2] Completed. Outputs saved to: ", normalizePath(output_dir))
message("[17-v2] Figures updated to journal style (no titles, smaller text, consistent colours).")

