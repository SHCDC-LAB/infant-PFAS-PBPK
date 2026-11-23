#!/usr/bin/env Rscript
# =============================================================================
# R/16_export_main_figures_1_2_NO_WATER.R
# Exports manuscript Figures 1 and 2 only. Figures 3 and 4 are exported by 17.
# No ODE rerun and no environmental-water inputs.
# =============================================================================
options(stringsAsFactors = FALSE)
pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "scales")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", unset = getwd()), winslash = "/", mustWork = FALSE)
VERSION <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
MODE <- tolower(Sys.getenv("PFAS_TRAJECTORY_MODE", "proxy"))
if (MODE != "proxy") VERSION <- paste0(VERSION, "_beta_", MODE)
PBPK <- Sys.getenv("PFAS_PBPK_RESULTS", unset = file.path(ROOT, "results", VERSION, "tables", "PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx"))
RISK <- Sys.getenv("PFAS_RISK_OUTPUT", unset = file.path(ROOT, "results", VERSION, "risk_characterization_NO_WATER", "infant_PFAS_risk_characterization_NO_WATER.xlsx"))
if (!file.exists(PBPK)) stop("Missing PBPK workbook: ", PBPK)
if (!file.exists(RISK)) stop("Missing risk workbook; run script 15 first: ", RISK)
OUT <- file.path(ROOT, "results", VERSION, "journal_figures_NO_WATER", "main")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

milk <- readxl::read_excel(RISK, sheet = "Milk_EDI_Individual") %>%
  mutate(Compound = factor(Compound, levels = c("PFOA","PFNA","PFDA","PFUdA","PFHxS","PFOS")),
         Period = factor(Period, levels = c("0-180 days","180-365 days","0-365 days")))
pbpk <- readxl::read_excel(PBPK, sheet = "PBPK_Results") %>%
  filter(Scenario == "S1") %>%
  mutate(Time_days = as.numeric(Time_days), ID = as.character(ID),
         Compound = as.character(Compound)) %>%
  filter(Time_days <= 365)

theme_journal <- function() theme_classic(base_size = 11, base_family = "sans") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 9.5),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom", panel.grid = element_blank())

p1 <- ggplot(milk, aes(x = Period, y = Mean_EDI_ng_kg_day, fill = Compound)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.45, linewidth = 0.35) +
  facet_wrap(~ Compound, scales = "free_y", ncol = 3) +
  scale_y_log10(labels = scales::label_number()) +
  labs(x = NULL, y = "Milk-derived EDI (ng/kg bw/day)") +
  theme_journal() + theme(legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1, size = 8.5))
ggsave(file.path(OUT, "Figure_1_milk_EDI_distribution.pdf"), p1, width = 7.2, height = 5.8)
ggsave(file.path(OUT, "Figure_1_milk_EDI_distribution.tiff"), p1, width = 7.2, height = 5.8, dpi = 900, compression = "lzw")

long <- pbpk %>%
  select(ID, Compound, Time_days, C_Plasma, C_Liver, C_Kidney) %>%
  pivot_longer(cols = starts_with("C_"), names_to = "Organ", values_to = "Concentration_ng_L") %>%
  mutate(Organ = recode(Organ, C_Plasma = "Plasma", C_Liver = "Liver", C_Kidney = "Kidney"))
summ <- long %>%
  group_by(Compound, Organ, Time_days) %>%
  summarise(Median = median(Concentration_ng_L, na.rm = TRUE),
            Q25 = quantile(Concentration_ng_L, .25, na.rm = TRUE),
            Q75 = quantile(Concentration_ng_L, .75, na.rm = TRUE), .groups = "drop")
summ$Compound <- factor(summ$Compound, levels = c("PFOA","PFNA","PFDA","PFUdA","PFHxS","PFOS"))
summ$Organ <- factor(summ$Organ, levels = c("Plasma","Liver","Kidney"))
org_cols <- c(Plasma = "#0072B2", Liver = "#D55E00", Kidney = "#009E73")
p2 <- ggplot(summ, aes(Time_days, Median, colour = Organ, fill = Organ)) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = .12, colour = NA) +
  geom_line(linewidth = .8) + facet_wrap(~Compound, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = org_cols) + scale_fill_manual(values = org_cols) +
  coord_cartesian(xlim = c(0,365)) +
  scale_x_continuous(breaks = c(0,90,180,270,365)) +
  labs(x = "Postnatal age (days)", y = "PBPK-predicted concentration (ng/L)", colour = "Compartment") +
  theme_journal()
ggsave(file.path(OUT, "Figure_2_organ_concentration_trajectories.pdf"), p2, width = 7.2, height = 6.2)
ggsave(file.path(OUT, "Figure_2_organ_concentration_trajectories.tiff"), p2, width = 7.2, height = 6.2, dpi = 900, compression = "lzw")
message("Figures 1 and 2 exported to: ", OUT)
