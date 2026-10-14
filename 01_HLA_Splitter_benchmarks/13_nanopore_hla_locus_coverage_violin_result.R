# Manuscript section: HLA-Splitter enables accurate demultiplexing across low-input settings
# Purpose: Final nanopore typing and reference-concordance figures.
# Original source: 14b_revision_nanopore_hla_locus_coverage_violin_result.R
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
  library(ggplot2)
  library(ggbeeswarm)
})

project_dir <- file.path(getOption("hla.data_root"), "30PBMC")
project_dir <- gsub("\\\\", "/", project_dir)
nanopore_dir <- file.path(project_dir, "Result", "nanopore")
dir.create(nanopore_dir, recursive = TRUE, showWarnings = FALSE)

hla_levels <- c("A", "B", "C", "DRB1", "DQB1", "DPB1")
hla_light_colors <- c(
  "A" = "#F9A8A3",
  "B" = "#D8CB64",
  "C" = "#6ED88C",
  "DRB1" = "#6ED9DC",
  "DQB1" = "#9FC2FF",
  "DPB1" = "#F6A5E8"
)

theme_hla <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(colour = "black", size = base_size),
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
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

save_plot_safe <- function(plot, filename, width, height) {
  pdf_path <- file.path(nanopore_dir, paste0(filename, ".pdf"))
  png_path <- file.path(nanopore_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  }
  file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
  file.copy(tmp_png, png_path, overwrite = TRUE)
  invisible(c(pdf = pdf_path, png = png_path))
}

message("Reading nanopore coverage from Result/nanopore...")
# Read the saved object or source table used by this analysis.
coverage <- read.csv(file.path(nanopore_dir, "coverage.csv"), check.names = FALSE)

coverage_long <- coverage %>%
  select(any_of(c("SampleName", "SampleID")), matches("^(A|B|C|DRB1|DQB1|DPB1)_")) %>%
  pivot_longer(
    cols = matches("^(A|B|C|DRB1|DQB1|DPB1)_"),
    names_to = c("HLA_allele", "Allele"),
    names_sep = "_",
    values_to = "HLA_count"
  ) %>%
  filter(HLA_allele %in% hla_levels) %>%
  mutate(
    HLA_count = suppressWarnings(as.numeric(HLA_count)),
    HLA_allele = factor(HLA_allele, levels = hla_levels),
    Genome_coverage = HLA_count
  ) %>%
  filter(!is.na(Genome_coverage), is.finite(Genome_coverage), Genome_coverage > 0)

# Export the analysis source data for figure reproduction.
write.csv(coverage_long, file.path(nanopore_dir, "HLA_locus_coverage_source.csv"), row.names = FALSE)

p_coverage <- ggplot(coverage_long, aes(x = HLA_allele, y = Genome_coverage, fill = HLA_allele, colour = HLA_allele)) +
  geom_violin(width = 0.68, alpha = 0.42, trim = FALSE, linewidth = 0.55, colour = "black") +
  ggbeeswarm::geom_quasirandom(width = 0.20, size = 1.45, alpha = 0.70) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.1, stroke = 0.35, fill = "white", colour = "black") +
  scale_fill_manual(values = hla_light_colors, name = "HLA allele") +
  scale_colour_manual(values = hla_light_colors, name = "HLA allele") +
  scale_y_log10(
    breaks = c(100, 200, 500, 1000, 2000, 5000),
    minor_breaks = c(seq(100, 900, by = 100), seq(1000, 9000, by = 1000)),
    labels = scales::label_comma(),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  labs(x = "HLA allele", y = "Genome coverage (X)") +
  theme_hla(11) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.55, colour = "black")),
    colour = "none"
  )

save_plot_safe(p_coverage, "HLA_locus coverage2", width = 4.20, height = 3.60)
message("Saved HLA_locus coverage2.pdf and .png in: ", nanopore_dir)
