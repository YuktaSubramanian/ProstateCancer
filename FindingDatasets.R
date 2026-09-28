# FindingDatasets.R
#
# Purpose:
#   For each GEO dataset listed below, this script:
#     1. Downloads the GSE (expression matrix + sample metadata) via GEOquery
#     2. Saves the raw ExpressionSet to DownloadedData/
#     3. Checks for time-to-event + event-status columns in the clinical data
#     4. Reports the data processing level (raw counts, log2, normalised, etc.)
#     5. Writes a per-dataset report to the console and to a combined text file
#
# Datasets:
#   GSE21034  GSE46691  GSE10645  GSE25136
#   GSE18916  GSE68907  GSE141551
#
# Output:
#   DownloadedData/dataset_report.txt

library(GEOquery)

# ── Paths ────────────────────────────────────────────────────────────────────

project_dir  <- "/home/yukta/Yukta!/Codes"
data_dir     <- file.path(
  project_dir, "ProstateCancer", "DownloadedData"
)
report_file  <- file.path(data_dir, "dataset_report.txt")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# ── Datasets to process ───────────────────────────────────────────────────────

gse_ids <- c(
  "GSE21034",
  "GSE46691",
  "GSE10645",
  "GSE25136",
  "GSE18916",
  "GSE68907",
  "GSE141551"
)

# ── Survival-column detection helpers ────────────────────────────────────────
# Keywords that suggest a time-to-event column
TIME_PATTERNS <- c(
  "time", "days", "months", "years", "survival",
  "follow", "duration", "os_", "_os", "pfs_", "_pfs",
  "rfs_", "_rfs", "dfs_", "_dfs", "bcr_", "_bcr",
  "recur", "relapse", "progression"
)

# Keywords that suggest an event / status column
EVENT_PATTERNS <- c(
  "status", "event", "death", "deceased", "vital",
  "censor", "outcome", "dead", "alive",
  "os_event", "pfs_event", "rfs_event", "dfs_event",
  "bcr_event", "recur", "relapse", "progression"
)

match_cols <- function(col_names, patterns) {
  hits <- character(0)
  for (cn in col_names) {
    cn_lower <- tolower(cn)
    if (any(vapply(patterns, function(p) grepl(p, cn_lower), logical(1)))) {
      hits <- c(hits, cn)
    }
  }
  hits
}

detect_survival_cols <- function(clinical_df) {
  cols <- colnames(clinical_df)
  list(
    time_cols  = match_cols(cols, TIME_PATTERNS),
    event_cols = match_cols(cols, EVENT_PATTERNS)
  )
}

# ── Processing level detector ─────────────────────────────────────────────────
# Guesses the level of pre-processing from value range and dataset metadata.
detect_processing_level <- function(eset) {
  expr    <- exprs(eset)
  val_min <- min(expr, na.rm = TRUE)
  val_max <- max(expr, na.rm = TRUE)

  level <- "Unknown"

  if (val_min >= 0 && val_max > 100) {
    level <- paste0(
      "Likely raw counts or RPKM/FPKM/TPM ",
      "(range: ", round(val_min, 2), " - ", round(val_max, 2), ")"
    )
  } else if (val_min >= 0 && val_max <= 30) {
    level <- paste0(
      "Likely log2-transformed (range: ",
      round(val_min, 2), " - ", round(val_max, 2), ")"
    )
  } else if (val_min < 0) {
    level <- paste0(
      "Likely log2 or z-score normalised (negative values present); ",
      "range: ", round(val_min, 2), " - ", round(val_max, 2)
    )
  }

  level
}

# ── Per-dataset reporting ─────────────────────────────────────────────────────

emit <- function(..., con = NULL) {
  msg <- paste0(...)
  message(msg)
  if (!is.null(con)) cat(msg, "\n", file = con, sep = "", append = TRUE)
}

report_dataset <- function(gse_id, con) {
  sep_line <- strrep("=", 70)
  emit(sep_line, con = con)
  emit("DATASET: ", gse_id, con = con)
  emit(sep_line, con = con)

  # -- Download ---------------------------------------------------------------
  emit("Downloading from GEO ...", con = con)
  gse_list <- tryCatch(
    getGEO(
      gse_id,
      GSEMatrix  = TRUE,
      destdir    = data_dir,
      AnnotGPL   = FALSE,
      getGPL     = FALSE
    ),
    error = function(e) {
      emit("  ERROR downloading ", gse_id, ": ", conditionMessage(e),
           con = con)
      return(NULL)
    }
  )

  if (is.null(gse_list) || length(gse_list) == 0) {
    emit("  SKIPPED (download failed)", con = con)
    return(invisible(NULL))
  }

  n_platforms <- length(gse_list)
  emit("  Platforms found: ", n_platforms, con = con)

  # Loop over platforms (usually just one)
  for (p_idx in seq_len(n_platforms)) {
    eset     <- gse_list[[p_idx]]
    pl_name  <- annotation(eset)
    if (nchar(pl_name) == 0) pl_name <- paste0("platform_", p_idx)

    emit("", con = con)
    emit("  --- Platform ", p_idx, ": ", pl_name, " ---", con = con)

    # -- Expression matrix summary -------------------------------------------
    expr        <- exprs(eset)
    n_genes     <- nrow(expr)
    n_samples   <- ncol(expr)
    emit("  Genes x Samples : ", n_genes, " x ", n_samples, con = con)
    emit("  Value range     : [",
         round(min(expr, na.rm = TRUE), 3), ", ",
         round(max(expr, na.rm = TRUE), 3), "]",
         con = con)
    emit("  Any NA          : ", anyNA(expr), con = con)

    # -- Processing level -----------------------------------------------------
    proc_level <- detect_processing_level(eset)
    emit("  Processing level: ", proc_level, con = con)

    # -- Clinical data --------------------------------------------------------
    clinical   <- pData(eset)
    all_cols   <- colnames(clinical)
    n_cols     <- length(all_cols)
    emit("", con = con)
    emit("  Clinical columns (", n_cols, " total):", con = con)

    # Print all column names (wrapped at 70 chars for readability)
    line_buf <- "    "
    for (cn in all_cols) {
      candidate <- paste0(cn, "  ")
      if (nchar(line_buf) + nchar(candidate) > 70) {
        emit(line_buf, con = con)
        line_buf <- "    "
      }
      line_buf <- paste0(line_buf, candidate)
    }
    if (nchar(trimws(line_buf)) > 0) emit(line_buf, con = con)

    # -- Survival detection ---------------------------------------------------
    surv_hits  <- detect_survival_cols(clinical)
    time_cols  <- surv_hits$time_cols
    event_cols <- surv_hits$event_cols

    emit("", con = con)
    emit("  >> Survival data check:", con = con)

    if (length(time_cols) > 0) {
      emit("     Time-to-event columns found (", length(time_cols), "):",
           con = con)
      for (tc in time_cols) {
        vals    <- suppressWarnings(as.numeric(clinical[[tc]]))
        n_valid <- sum(!is.na(vals))
        emit("       - ", tc,
             "  [non-NA: ", n_valid, "/", nrow(clinical), "]",
             con = con)
      }
    } else {
      emit("     No time-to-event columns detected", con = con)
    }

    if (length(event_cols) > 0) {
      emit("     Event/status columns found (", length(event_cols), "):",
           con = con)
      for (ec in event_cols) {
        uniq_vals <- unique(clinical[[ec]])
        uniq_vals <- head(uniq_vals[!is.na(uniq_vals)], 10)
        emit("       - ", ec,
             "  [unique values: ",
             paste(uniq_vals, collapse = " | "), "]",
             con = con)
      }
    } else {
      emit("     No event/status columns detected", con = con)
    }

    # -- Verdict ------------------------------------------------------------
    has_survival <- length(time_cols) > 0 && length(event_cols) > 0
    emit("", con = con)
    emit("  >> Has usable survival data: ",
         ifelse(has_survival, "YES", "NO (missing time or event column)"),
         con = con)
  }

  emit("", con = con)
}

# ── Run ───────────────────────────────────────────────────────────────────────

# Open report file (overwrite if exists)
con <- file(report_file, open = "w")
on.exit(close(con), add = TRUE)

cat(
  "Dataset Suitability Report\n",
  "Generated: ", format(Sys.time()), "\n",
  strrep("=", 70), "\n\n",
  file = con, sep = ""
)
message("Writing report to: ", report_file)

for (gse_id in gse_ids) {
  report_dataset(gse_id, con)
}

emit(strrep("=", 70), con = con)
emit("Done. All datasets processed.", con = con)
message("Report saved to: ", report_file)
