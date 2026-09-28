# This script runs dorothea library to identify genes that 
# are regulated by each transcription factor of interest.
# We can then use this list as a signature for future analyses 
# of transcriptomic data.

# dorothea was installed using bioconda already
# the libraries

# Working directoy
proj_dir <- "/home/yukta/Yukta!/Codes/ProstateCancer"

library(dorothea)
library(dplyr)
library(tidyr)
library(readr)

# Create output directory with timestamp
output_dir <-  paste0(proj_dir, "DownloadedData/Dorothea")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created output directory:", output_dir, "\n")
}

# Load human DoRothEA interactions
data("dorothea_hs", package = "dorothea")

# List of TFs of interest
# The coloquial name used in the study may differ from the proper gene name 
# So create a dictionary
tf_dict <- list(
  "OCT4" = "POU5F1",
  "LIN28" = "LIN28A",  # Also consider LIN28B if needed
  "P63" = "TP63",
  "NKX3-1" = "NKX3-1",  # Prostate-specific homeobox gene
  "NANOG" = "NANOG",
  "SOX2" = "SOX2",
  "AR" = "AR",  # Androgen receptor
  "GATA2" = "GATA2",
  "FOXA1" = "FOXA1"
)

# Function to call for each TF gene
get_tf_signatures <- function(tf_name) {
  tf_network <- dorothea_hs |>
    dplyr::filter(
      tf == tf_name,
      confidence %in% c("A", "B", "C", "D", "E")  # Include all confidence levels
    )
  
  if (nrow(tf_network) == 0) {
    warning("No interactions found for ", tf_name, " in dorothea database.")
    return(list(
      TF = tf_name,
      positive_signature = character(0),
      negative_signature = character(0),
      total_interactions = 0
    ))
  }
  
  positive_signature <- tf_network |>
    dplyr::filter(mor > 0) |>
    dplyr::distinct(target) |>
    dplyr::pull(target) |>
    sort()
  
  negative_signature <- tf_network |>
    dplyr::filter(mor < 0) |>
    dplyr::distinct(target) |>
    dplyr::pull(target) |>
    sort()
  
  list(
    TF = tf_name,
    positive_signature = positive_signature,
    negative_signature = negative_signature,
    total_interactions = nrow(tf_network),
    confidence_summary = table(tf_network$confidence)
  )
}

# save the data separately for each TF in a readable format (csv)
# Iterate over the tf_dict and pass corresponding gene name in the dictionary
all_signatures <- list()

for (tf_name in names(tf_dict)) {
  gene_name <- tf_dict[[tf_name]]
  cat("Processing", tf_name, "(", gene_name, ")\n")
  
  signatures <- get_tf_signatures(gene_name)
  all_signatures[[tf_name]] <- signatures
  
  # Save positive signatures
  if (length(signatures$positive_signature) > 0) {
    positive_df <- data.frame(
      TF = gene_name,
      gene = signatures$positive_signature,
      regulation = "positive"
    )
    write.csv(
      positive_df,
      file.path(output_dir, paste0(gene_name, "_positive_signatures.csv")),
      row.names = FALSE
    )
  }
  
  # Save negative signatures
  if (length(signatures$negative_signature) > 0) {
    negative_df <- data.frame(
      TF = gene_name,
      gene = signatures$negative_signature,
      regulation = "negative"
    )
    write.csv(
      negative_df,
      file.path(output_dir, paste0(gene_name, "_negative_signatures.csv")),
      row.names = FALSE
    )
  }
  
  # Save complete table with all interactions
  complete_table <- dorothea_hs |>
    dplyr::filter(
      tf == gene_name,
      confidence %in% c("A", "B")
    ) |>
    dplyr::transmute(
      TF = tf,
      gene = target,
      regulation = dplyr::if_else(mor > 0, "positive", "negative"),
      confidence,
      mor
    ) |>
    dplyr::distinct()
  
  write.csv(
    complete_table,
    file.path(output_dir, paste0(gene_name, "_DoRothEA_AB_signatures.csv")),
    row.names = FALSE
  )
}

# Print summary
cat("\n=== Summary of extracted signatures ===\n")

# Create summary dataframe
summary_df <- data.frame(
  TF = character(),
  Gene_Symbol = character(),
  Positive_Genes = integer(),
  Negative_Genes = integer(),
  Total_Interactions = integer(),
  Has_Data = logical(),
  stringsAsFactors = FALSE
)

for (tf_name in names(all_signatures)) {
  sig <- all_signatures[[tf_name]]
  
  summary_row <- data.frame(
    TF = tf_name,
    Gene_Symbol = sig$TF,
    Positive_Genes = length(sig$positive_signature),
    Negative_Genes = length(sig$negative_signature),
    Total_Interactions = ifelse(is.null(sig$total_interactions), 0, sig$total_interactions),
    Has_Data = length(sig$positive_signature) > 0 || length(sig$negative_signature) > 0,
    stringsAsFactors = FALSE
  )
  
  summary_df <- rbind(summary_df, summary_row)
  
  cat(tf_name, "(", sig$TF, "):\n")
  cat("  Positive signatures:", length(sig$positive_signature), "genes\n")
  cat("  Negative signatures:", length(sig$negative_signature), "genes\n")
  cat("  Total interactions:", ifelse(is.null(sig$total_interactions), 0, sig$total_interactions), "\n")
  
  # Print confidence levels if available
  if (!is.null(sig$confidence_summary)) {
    cat("  Confidence levels:", paste(names(sig$confidence_summary), "=", sig$confidence_summary, collapse = ", "), "\n")
  }
  
  # Print first few genes if available
  if (length(sig$positive_signature) > 0) {
    first_pos <- head(sig$positive_signature, 3)
    cat("  First positive genes:", paste(first_pos, collapse = ", "))
    if (length(sig$positive_signature) > 3) cat(" ...\n") else cat("\n")
  }
  if (length(sig$negative_signature) > 0) {
    first_neg <- head(sig$negative_signature, 3)
    cat("  First negative genes:", paste(first_neg, collapse = ", "))
    if (length(sig$negative_signature) > 3) cat(" ...\n") else cat("\n")
  }
  
  if (length(sig$positive_signature) == 0 && length(sig$negative_signature) == 0) {
    cat("  WARNING: No interactions found in dorothea database\n")
  }
  cat("\n")
}

# Save summary to CSV
write.csv(summary_df, file.path(output_dir, "TF_signatures_summary.csv"), row.names = FALSE)
cat("Summary saved to", file.path(output_dir, "TF_signatures_summary.csv"), "\n")

# Save the complete results as RData for later use
save(all_signatures, file = file.path(output_dir, "TF_signatures_results.RData"))
cat("Results saved to", file.path(output_dir, "TF_signatures_results.RData"), "\n")

# Also save a combined CSV with all signatures
all_signatures_df <- data.frame()
for (tf_name in names(all_signatures)) {
  sig <- all_signatures[[tf_name]]
  if (length(sig$positive_signature) > 0) {
    pos_df <- data.frame(
      TF = sig$TF,
      TF_colloquial = tf_name,
      gene = sig$positive_signature,
      regulation = "positive",
      stringsAsFactors = FALSE
    )
    all_signatures_df <- rbind(all_signatures_df, pos_df)
  }
  if (length(sig$negative_signature) > 0) {
    neg_df <- data.frame(
      TF = sig$TF,
      TF_colloquial = tf_name,
      gene = sig$negative_signature,
      regulation = "negative",
      stringsAsFactors = FALSE
    )
    all_signatures_df <- rbind(all_signatures_df, neg_df)
  }
}

if (nrow(all_signatures_df) > 0) {
  write.csv(all_signatures_df, file.path(output_dir, "all_TF_signatures_combined.csv"), row.names = FALSE)
  cat("Combined signatures saved to", file.path(output_dir, "all_TF_signatures_combined.csv"), "\n")
}

cat("\n=== Script Execution Complete ===\n")
cat("Files created:\n")
cat("1. Individual CSVs for each TF (positive/negative/complete)\n")
cat("2. TF_signatures_summary.csv - Summary of all TFs\n")
cat("3. TF_signatures_results.RData - Complete R object\n")
cat("4. all_TF_signatures_combined.csv - All signatures in one file\n")
cat("\nTo use these signatures in transcriptomic analysis:\n")
cat("- Positive signatures: genes upregulated by the TF\n")
cat("- Negative signatures: genes downregulated by the TF\n")
cat("- Use in GSEA, GSVA, or other gene set enrichment methods\n")