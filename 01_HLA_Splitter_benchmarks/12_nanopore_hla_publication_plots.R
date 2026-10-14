# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Final nanopore typing and reference-concordance figures.
# Original source: 14_revision_nanopore_hla_publication_plots.R
# Figure/table scope: FigS1c; FigS1d; FigS1e; FigS1f
# Data cohort: 30PBMC. Inputs retain the original Data/Result hierarchy.
# Run with run_analysis.R; configure HLA_SPLITTER_DATA_ROOT before execution.
# Input objects, metadata and source tables are read from the configured cohort.
# Outputs: files written by this script under the cohort Figures/Result directories.
# Preserve named color mappings, annotation columns and statistical units.
# See README.md and manuscript_consistency_review.md for release limitations.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggbeeswarm)
  library(readr)
library(gtools)
})

`%||%` <- function(x, y) if (length(x) == 0 || is.na(x) || is.null(x)) y else x

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", file_arg[1] %||% file.path(getwd(), "14_revision_nanopore_hla_publication_plots.R"))
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
pbmc_dir <- getOption("hla.cohort_dir")
nanopore_dir <- file.path(pbmc_dir, "Data", "HLA_list", "nanopore")

hla_levels <- c("A", "B", "C", "DRB1", "DQB1", "DPB1")
hla_colors <- c(
  "A" = "#F8766D",
  "B" = "#B79F00",
  "C" = "#00BA38",
  "DRB1" = "#00BFC4",
  "DQB1" = "#619CFF",
  "DPB1" = "#F564E3"
)
hla_light_colors <- c(
  "A" = "#F9A8A3",
  "B" = "#D8CB64",
  "C" = "#6ED88C",
  "DRB1" = "#6ED9DC",
  "DQB1" = "#9FC2FF",
  "DPB1" = "#F6A5E8"
)

truncate_hla <- function(allele, field = 2) {
  allele <- as.character(allele)
  ifelse(
    is.na(allele) | allele == "",
    NA_character_,
    str_extract(allele, if (field == 1) "^[A-Z0-9]+\\*[0-9]+" else "^[A-Z0-9]+\\*[0-9]+:[0-9]+")
  )
}

theme_hla <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(colour = "black", size = base_size),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.45),
      panel.grid.minor = element_line(colour = "grey94", linewidth = 0.25),
      panel.border = element_rect(fill = NA, colour = "grey35", linewidth = 0.55),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = base_size + 0.8),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.42, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

save_plot <- function(plot, file_base, width = 5.2, height = 4.0) {
  ggsave(file.path(nanopore_dir, paste0(file_base, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(nanopore_dir, paste0(file_base, ".png")), plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
}

# Read the saved object or source table used by this analysis.
typing <- read_csv(file.path(nanopore_dir, "HLAtyping result.csv"), show_col_types = FALSE)

typing_long <- typing %>%
  pivot_longer(
    cols = -c(`Sample ID`, Method),
    names_to = c("HLA_allele", "Allele"),
    names_sep = "_",
    values_to = "HLA_type"
  ) %>%
  filter(HLA_allele %in% hla_levels) %>%
  mutate(
    HLA_2field = truncate_hla(HLA_type, 2),
    HLA_1field = truncate_hla(HLA_type, 1)
  )

accuracy_call_df <- typing_long %>%
  group_by(`Sample ID`, HLA_allele, Allele) %>%
  summarise(
    SBT_2field = HLA_2field[Method == "SBT"][1],
    ONT_2field = HLA_2field[Method == "ONT-Rapid"][1],
    SBT_1field = HLA_1field[Method == "SBT"][1],
    ONT_1field = HLA_1field[Method == "ONT-Rapid"][1],
    .groups = "drop"
  ) %>%
  mutate(
    Concordance = case_when(
      is.na(SBT_2field) | is.na(ONT_2field) ~ "Missing",
      SBT_2field == ONT_2field ~ "2-field match",
      SBT_1field == ONT_1field ~ "1-field match",
      TRUE ~ "Mismatch"
    ),
    Is_1field_match = !is.na(SBT_1field) & !is.na(ONT_1field) & SBT_1field == ONT_1field,
    Is_2field_match = !is.na(SBT_2field) & !is.na(ONT_2field) & SBT_2field == ONT_2field
  ) %>%
  mutate(HLA_allele = factor(HLA_allele, levels = hla_levels))

accuracy_df <- accuracy_call_df %>%
  group_by(HLA_allele) %>%
  summarise(
    `1-field` = mean(Is_1field_match, na.rm = TRUE) * 100,
    `2-field` = mean(Is_2field_match, na.rm = TRUE) * 100,
    Tested_alleles = n(),
    .groups = "drop"
  ) %>%
  pivot_longer(c(`1-field`, `2-field`), names_to = "Resolution", values_to = "Accuracy_pct") %>%
  mutate(
    HLA_allele = factor(HLA_allele, levels = hla_levels),
    Resolution = factor(Resolution, levels = c("1-field", "2-field"))
  )

p_accuracy <- ggplot(accuracy_df, aes(x = HLA_allele, y = Accuracy_pct, fill = Resolution)) +
  geom_col(
    position = position_dodge(width = 0.70),
    width = 0.58,
    linewidth = 0.55,
    colour = "black",
    alpha = 0.72
  ) +
  scale_fill_manual(values = c("1-field" = "#5DA5DA", "2-field" = "#F17CB0"), name = "Resolution") +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Nanopore HLA typing accuracy", x = "HLA allele", y = "Typing concordance") +
  theme_hla(12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

coverage <- read_csv(file.path(nanopore_dir, "coverage.csv"), show_col_types = FALSE)

coverage_long <- coverage %>%
  pivot_longer(
    cols = -SampleID,
    names_to = c("HLA_allele", "Allele"),
    names_sep = "_",
    values_to = "HLA_count"
  ) %>%
  filter(HLA_allele %in% hla_levels) %>%
  mutate(
    HLA_allele = factor(HLA_allele, levels = hla_levels),
    log10_count = log10(as.numeric(HLA_count) + 1)
  )

p_coverage <- ggplot(coverage_long, aes(x = HLA_allele, y = log10_count, fill = HLA_allele, colour = HLA_allele)) +
  geom_violin(width = 0.72, alpha = 0.44, trim = FALSE, linewidth = 0.55, colour = "black") +
  ggbeeswarm::geom_quasirandom(width = 0.22, size = 1.45, alpha = 0.70) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.1, stroke = 0.35, fill = "white", colour = "black") +
  scale_fill_manual(values = hla_light_colors, name = "HLA allele") +
  scale_colour_manual(values = hla_light_colors, name = "HLA allele") +
  labs(title = NULL, x = "HLA allele", y = expression(Log[10]~"(HLA count +1)")) +
  theme_hla(12) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.55, colour = "black")),
    colour = "none"
  )

# Export the analysis source data for figure reproduction.
write_csv(accuracy_call_df, file.path(nanopore_dir, "HLA_accuracy_call_source.csv"))
write_csv(accuracy_df, file.path(nanopore_dir, "HLA_accuracy_1field_2field_summary_source.csv"))
write_csv(coverage_long, file.path(nanopore_dir, "HLA_locus_coverage_source.csv"))

save_plot(p_accuracy, "HLA_accuracy", width = 5.0, height = 3.4)
save_plot(p_coverage, "HLA_locus coverage2", width = 5.2, height = 4.0)

message("Updated nanopore HLA plots in: ", nanopore_dir)
