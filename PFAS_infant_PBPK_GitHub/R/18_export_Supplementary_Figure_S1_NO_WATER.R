#!/usr/bin/env Rscript
# =============================================================================
# R/18_export_Supplementary_Figure_S1_NO_WATER.R
# Derives compound ranks from the actual EDI and frozen PBPK outputs.
# No manually transcribed rank matrix; no ODE rerun; no water data.
# =============================================================================
options(stringsAsFactors = FALSE)
pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "tibble", "openxlsx")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", unset = getwd()), winslash = "/", mustWork = FALSE)
VERSION <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
MODE <- tolower(Sys.getenv("PFAS_TRAJECTORY_MODE", "proxy"))
if (MODE != "proxy") VERSION <- paste0(VERSION, "_beta_", MODE)
PBPK <- Sys.getenv("PFAS_PBPK_RESULTS", unset = file.path(ROOT, "results", VERSION, "tables", "PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx"))
RISK <- Sys.getenv("PFAS_RISK_OUTPUT", unset = file.path(ROOT, "results", VERSION, "risk_characterization_NO_WATER", "infant_PFAS_risk_characterization_NO_WATER.xlsx"))
CONTRIB <- Sys.getenv("PFAS_17_OUTPUT", unset = file.path(ROOT, "results", VERSION, "IBQ_organ_contribution_figures_NO_WATER", "PFAS_R_BMDL10_organ_contribution_results.xlsx"))
for (f in c(PBPK,RISK,CONTRIB)) if (!file.exists(f)) stop("Required input missing: ", f)
OUT <- file.path(ROOT, "results", VERSION, "supplementary_figure_S1_NO_WATER")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
chem <- c("PFOA","PFNA","PFDA","PFUdA","PFHxS","PFOS")
edi <- readxl::read_excel(RISK, sheet = "Milk_EDI_Individual") %>%
  filter(Period == "0-365 days") %>%
  group_by(Compound) %>% summarise(Value = median(Mean_EDI_ng_kg_day), .groups = "drop") %>%
  mutate(Metric = "Milk EDI")
pbpk <- readxl::read_excel(PBPK, sheet = "PBPK_Results") %>% filter(Scenario == "S1")
day365 <- pbpk %>% filter(Compound %in% chem) %>%
  group_by(ID, Compound) %>%
  group_modify(~ tibble(Value = approx(as.numeric(.x$Time_days), .x$C_Plasma, xout = 365)$y)) %>%
  ungroup() %>% group_by(Compound) %>%
  summarise(Value = median(Value), .groups = "drop") %>%
  mutate(Metric = "Plasma day 365")
org <- readxl::read_excel(CONTRIB, sheet = "Organ_Contribution_Individual") %>%
  group_by(Compound, Organ) %>% summarise(Value = median(CTWA_0_365_ng_L), .groups = "drop") %>%
  mutate(Metric = paste0(Organ, " CTWA"))
values <- bind_rows(edi, day365, org) %>% filter(Compound %in% chem)
ranked <- values %>% group_by(Metric) %>%
  mutate(Rank = rank(-Value, ties.method = "min")) %>% ungroup()
ranked$Compound <- factor(ranked$Compound, levels = rev(chem))
ranked$Metric <- factor(ranked$Metric, levels = c("Milk EDI","Plasma day 365","Plasma CTWA","Liver CTWA","Kidney CTWA"))
p <- ggplot(ranked, aes(Metric, Compound, fill = Rank)) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = Rank), size = 3.4, fontface = "bold") +
  scale_fill_gradientn(colours = c("#084594","#2171B5","#4292C6","#9ECAE1","#DEEBF7","#F7FBFF"),
                       limits = c(1,6), breaks = 1:6, name = "Rank\n(1 = highest)") +
  labs(x = NULL, y = NULL) + theme_classic(base_size = 10.5) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9.5), panel.grid = element_blank(),
        legend.position = "right")
ggsave(file.path(OUT,"Supplementary_Figure_S1_exposure_rankings.pdf"),p,width=7.2,height=4.5)
ggsave(file.path(OUT,"Supplementary_Figure_S1_exposure_rankings.tiff"),p,width=7.2,height=4.5,dpi=900,compression="lzw")
openxlsx::write.xlsx(list(Rank_Values = ranked), file.path(OUT,"Supplementary_Figure_S1_rank_values.xlsx"), overwrite=TRUE)
message("Supplementary Figure S1 exported: ", OUT)
