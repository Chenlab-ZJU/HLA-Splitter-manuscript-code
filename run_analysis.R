#!/usr/bin/env Rscript
# Configure manuscript paths and run one retained analysis in a fresh R process.
# The runner never installs packages or reruns upstream analyses automatically.

main <- function() {
  # Use UTF-8 when the parent project directory contains non-ASCII characters.
  if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
    suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
  }
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- sub("^--file=", "", file_arg[[1]])
  if (!grepl("^(/|[A-Za-z]:)", script)) script <- file.path(getwd(), script)
  code_root <- gsub("\\\\", "/", dirname(script))
  args <- commandArgs(TRUE)
  manifest_path <- file.path(code_root, "manuscript_code_manifest.csv")
  if (file.exists(manifest_path)) {
    manifest_text <- readLines(manifest_path, encoding = "UTF-8", warn = FALSE)
    manifest_text[[1]] <- sub("^\ufeff", "", manifest_text[[1]])
    manifest <- read.csv(text = paste(manifest_text, collapse = "\n"),
                         stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  } else {
    # Discover scripts from their metadata headers when no external index is present.
    files <- list.files(code_root, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
    relative <- substring(gsub("\\\\", "/", files), nchar(code_root) + 2L)
    keep <- grepl("^0[0-5]_[^/]+/", relative)
    files <- files[keep]
    relative <- relative[keep]
    rows <- lapply(seq_along(files), function(i) {
      header <- readLines(files[[i]], n = 30L, encoding = "UTF-8", warn = FALSE)
      cohort <- grep("^# Data cohort:", header, value = TRUE)
      role <- grep("^# Purpose:", header, value = TRUE)
      if (!length(cohort)) stop("Missing Data cohort header: ", relative[[i]])
      data.frame(file = relative[[i]],
                 cohort = sub("^# Data cohort: (.*)\\. Inputs.*$", "\\1", cohort[[1]]),
                 role = if (length(role)) sub("^# Purpose: ", "", role[[1]]) else "Analysis script",
                 stringsAsFactors = FALSE)
    })
    if (!length(rows)) stop("No analysis scripts found under ", code_root)
    manifest <- do.call(rbind, rows)
  }
  if (!length(args) || identical(args[[1]], "--list")) {
    cat(paste(manifest$file, manifest$role, sep = "  :  "), sep = "\n")
    return(invisible(NULL))
  }
  mode <- if (args[[1]] %in% c("--check", "--self-test")) args[[1]] else "run"
  if (mode == "--self-test") {
    options(hla.code_root = code_root)
    algorithm <- manifest$file[grepl("_HLA_demultiplex_v2.R$", manifest$file)]
    env <- new.env(parent = globalenv())
    sys.source(file.path(code_root, algorithm), envir = env)
    stopifnot(is.function(env$HLA_demultiplex))
    cat("PASS: manuscript manifest and retained HLA-Splitter function load successfully.\n")
    return(invisible(NULL))
  }
  if (mode == "--check" && length(args) < 2L) {
    stop("Usage: Rscript run_analysis.R --check <relative-script-path>")
  }
  target <- if (mode == "--check") args[[2]] else args[[1]]
  target <- gsub("\\\\", "/", target)
  row <- manifest[manifest$file == target, , drop = FALSE]
  if (nrow(row) != 1L) stop("Select a script listed by --list: ", target)
  revision_root <- Sys.getenv("HLA_SPLITTER_REVISION_ROOT", dirname(code_root))
  project_root <- Sys.getenv("HLA_SPLITTER_PROJECT_ROOT", dirname(revision_root))
  bundled_data <- file.path(revision_root, "Data")
  default_data <- if (dir.exists(bundled_data)) bundled_data else file.path(project_root, "Data")
  data_root <- Sys.getenv("HLA_SPLITTER_DATA_ROOT", default_data)
  data_root <- gsub("\\\\", "/", data_root)
  if (!dir.exists(data_root)) stop("Missing HLA_SPLITTER_DATA_ROOT: ", data_root)
  cohort_dir <- file.path(data_root, row$cohort[[1]])
  if (row$cohort[[1]] == "shared") cohort_dir <- data_root
  if (!dir.exists(cohort_dir)) stop("Missing cohort directory: ", cohort_dir)
  options(hla.code_root = code_root, hla.data_root = data_root,
          hla.project_root = project_root, hla.revision_root = revision_root,
          hla.cohort_dir = cohort_dir, hla.script = target)
  # Environment names preserve the input configuration used by the original scripts.
  vars <- c(SYN_PROJECT_DIR = "synthetic", CSF_PROJECT_DIR = "csf",
            LAVAGE_PROJECT_DIR = "balf", URINE_PROJECT_DIR = "urine",
            GSE189536_PROJECT_DIR = "gse189536", GSE189536_DATA_DIR = "gse189536",
            PBMC_PROJECT_DIR = "pbmc30")
  do.call(Sys.setenv, as.list(setNames(file.path(data_root, vars), names(vars))))
  Sys.setenv(REVISION_PROJECT_DIR = revision_root)
  file <- file.path(code_root, target)
  parse(file, encoding = "UTF-8")
  if (mode == "--check") {
    cat("PASS: selected script parses; cohort directory and path configuration exist.\n")
    cat("This preflight does not validate every dynamic input or execute the analysis.\n")
    return(invisible(NULL))
  }
  # Input-only bundles do not contain the output folders assumed by older scripts.
  for (folder in c("Figures", "Result")) {
    dir.create(file.path(cohort_dir, folder), recursive = TRUE, showWarnings = FALSE)
  }
  setwd(cohort_dir)
  set.seed(2026)
  message("Running ", target, "\nInput/output cohort: ", cohort_dir)
  # A separate script environment contains legacy rm(list = ls()) calls. Persistent
  # path options survive cleanup; source() helpers retain their original behavior.
  analysis_env <- new.env(parent = globalenv())
  sys.source(file, envir = analysis_env)
  message("Completed ", target)
}
main()
