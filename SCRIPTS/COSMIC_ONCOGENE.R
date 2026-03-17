# COSMIC Oncogene Enrichment Analysis

# Download COSMIC Cancer Gene Census if not present
 
download_cosmic_data <- function(cosmic_dir = "COSMIC_data") {
  
  # Create directory if it doesn't exist
  if (!dir.exists(cosmic_dir)) {
    dir.create(cosmic_dir, recursive = TRUE)
    message("Created COSMIC data directory: ", cosmic_dir)
  }
  
  # Define file path
  cosmic_file <- file.path(cosmic_dir, "cancer_gene_census.csv")
  
  # Check if file already exists and is valid
  if (file.exists(cosmic_file)) {
    # Check if file is valid
    first_line <- readLines(cosmic_file, n = 1)
    if (grepl("<!DOCTYPE|<html|<HTML", first_line)) {
      message("Existing file appears to be an HTML error page. Deleting and re-downloading...")
      file.remove(cosmic_file)
    } else {
      message("COSMIC Cancer Gene Census file already exists at: ", cosmic_file)
      return(cosmic_file)
    }
  }
  
  # Provide manual download instructions
  message("\n=======================================================")
  message("COSMIC requires registration for bulk downloads.")
  message("Please manually download the Cancer Gene Census:")
  message("1. Go to: https://cancer.sanger.ac.uk/census")
  message("2. Register/Login (free for academic use)")
  message("3. Download 'cancer_gene_census.csv'")
  message("4. Save it to: ", cosmic_file)
  message("=======================================================\n")
  
  stop("Please download COSMIC Cancer Gene Census manually and save to: ", cosmic_file)
}

# Load and process COSMIC Cancer Gene Census data

load_cosmic_data <- function(cosmic_file) {
  
  if (!file.exists(cosmic_file)) {
    stop("COSMIC file not found at: ", cosmic_file)
  }
  
  message("Loading COSMIC Cancer Gene Census data...")
  
  # Check if file is valid
  first_line <- readLines(cosmic_file, n = 1)
  if (grepl("<!DOCTYPE|<html|<HTML", first_line)) {
    stop("COSMIC file appears to be an HTML error page, not actual data. Please download manually.")
  }
  
  # Read the file
  cosmic_data <- tryCatch({
    read.csv(cosmic_file, header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    # Try tab-separated if comma-separated fails
    read.delim(cosmic_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  })
  
  message("Loaded ", nrow(cosmic_data), " cancer genes from COSMIC")
  message("COSMIC data columns: ", paste(names(cosmic_data), collapse = ", "))
  
  return(cosmic_data)
}

# Perform Kruskal-Wallis test with Dunn's post-hoc for a single gene

perform_kw_dunn_test <- function(gene_counts, groups) {
  
  # Kruskal-Wallis test
  kw_result <- kruskal.test(gene_counts ~ groups)
  
  # If KW is significant, perform Dunn's test
  dunn_result <- NULL
  dunn_24_mock_padj <- NA
  dunn_48_mock_padj <- NA
  sig_24_vs_mock <- FALSE
  sig_48_vs_mock <- FALSE
  
  if (kw_result$p.value < 0.05) {
    dunn_result <- dunnTest(gene_counts ~ groups, method = "bonferroni")
    
    # Extract comparisons
    dunn_table <- dunn_result$res
    
    # Find 24hr vs Mock comparison
    mock_24_comparison <- dunn_table[grepl("Mock.*24hr|24hr.*Mock", dunn_table$Comparison), ]
    if (nrow(mock_24_comparison) > 0) {
      dunn_24_mock_padj <- mock_24_comparison$P.adj[1]
      sig_24_vs_mock <- dunn_24_mock_padj < 0.05
    }
    
    # Find 48hr vs Mock comparison
    mock_48_comparison <- dunn_table[grepl("Mock.*48hr|48hr.*Mock", dunn_table$Comparison), ]
    if (nrow(mock_48_comparison) > 0) {
      dunn_48_mock_padj <- mock_48_comparison$P.adj[1]
      sig_48_vs_mock <- dunn_48_mock_padj < 0.05
    }
  }
  
  return(list(
    kw_pvalue = kw_result$p.value,
    dunn_24_mock_padj = dunn_24_mock_padj,
    dunn_48_mock_padj = dunn_48_mock_padj,
    sig_24_vs_mock = sig_24_vs_mock,
    sig_48_vs_mock = sig_48_vs_mock
  ))
}

# Main COSMIC Oncogene Enrichment Analysis Function

cosmic_oncogene_analysis <- function(dgelist,
                                     results_24hr,
                                     results_48hr,
                                     rpkm_matrix,
                                     cosmic_dir = "COSMIC_data",
                                     output_file = "COSMIC_oncogene_analysis.tsv",
                                     fdr_threshold = 0.05) {
  
  message("\n========== Starting COSMIC Oncogene Enrichment Analysis ==========\n")
  
  # 1. Download/load COSMIC data
  cosmic_file <- download_cosmic_data(cosmic_dir)
  cosmic_data <- load_cosmic_data(cosmic_file)
  
  # 2. Prepare gene lists
  # Extract gene symbols from COSMIC
  cosmic_gene_col <- if("Gene.Symbol" %in% names(cosmic_data)) {
    "Gene.Symbol"
  } else if("Gene Symbol" %in% names(cosmic_data)) {
    "Gene Symbol"
  } else if("GENE_SYMBOL" %in% names(cosmic_data)) {
    "GENE_SYMBOL"
  } else if("Gene_Symbol" %in% names(cosmic_data)) {
    "Gene_Symbol"
  } else if("SYMBOL" %in% names(cosmic_data)) {
    "SYMBOL"
  } else {
    # Print available columns
    message("Available COSMIC columns: ", paste(names(cosmic_data), collapse = ", "))
    stop("Could not find gene symbol column in COSMIC data. Please check column names.")
  }
  
  message("Using COSMIC column: ", cosmic_gene_col)
  cosmic_genes <- unique(cosmic_data[[cosmic_gene_col]])
  cosmic_genes <- cosmic_genes[!is.na(cosmic_genes) & cosmic_genes != ""]
  message("Found ", length(cosmic_genes), " unique cancer genes in COSMIC")
  
  # Print sample of COSMIC genes for debugging
  message("Sample COSMIC genes: ", paste(head(cosmic_genes, 10), collapse = ", "))
  
  # 3. Match genes - try multiple matching strategies
  dataset_genes <- rownames(rpkm_matrix)
  message("Your dataset has ", length(dataset_genes), " genes")
  message("Sample dataset genes: ", paste(head(dataset_genes, 10), collapse = ", "))
  
  # Strategy 1: Direct match
  analyzed_genes <- intersect(dataset_genes, cosmic_genes)
  message("Direct match found ", length(analyzed_genes), " overlapping genes")
  
  # Strategy 2: If no direct matches, try extracting base gene symbols
  # (in case your genes have "_ENSG..." suffix)
  if (length(analyzed_genes) == 0) {
    message("No direct matches. Trying to extract base gene symbols...")
    
    # Extract base symbols
    dataset_base_symbols <- sapply(strsplit(dataset_genes, "_"), function(x) x[1])
    
    # Find matches
    base_matches <- intersect(dataset_base_symbols, cosmic_genes)
    message("Found ", length(base_matches), " matches using base symbols")
    
    if (length(base_matches) > 0) {
      # Map back to full gene names in dataset
      analyzed_genes <- dataset_genes[dataset_base_symbols %in% base_matches]
      message("Mapped to ", length(analyzed_genes), " genes in your dataset")
    }
  }
  
  # Strategy 3: Try case-insensitive matching
  if (length(analyzed_genes) == 0) {
    message("Trying case-insensitive matching...")
    cosmic_genes_upper <- toupper(cosmic_genes)
    dataset_genes_upper <- toupper(dataset_genes)
    matches_upper <- intersect(dataset_genes_upper, cosmic_genes_upper)
    
    if (length(matches_upper) > 0) {
      analyzed_genes <- dataset_genes[toupper(dataset_genes) %in% matches_upper]
      message("Found ", length(analyzed_genes), " case-insensitive matches")
    }
  }
  
  if (length(analyzed_genes) == 0) {
    message("\n=== DIAGNOSTIC INFORMATION ===")
    message("COSMIC genes (first 20): ", paste(head(cosmic_genes, 20), collapse = ", "))
    message("Your genes (first 20): ", paste(head(dataset_genes, 20), collapse = ", "))
    message("==============================\n")
    stop("No overlap between COSMIC genes and your dataset. Check gene naming conventions above.")
  }
  
  message("Successfully matched ", length(analyzed_genes), " cancer genes for analysis\n")
  
  # 4. Prepare group assignments from dgelist
  groups <- factor(dgelist$samples$group)
  
  # 5. Perform statistical tests for each cancer gene
  message("\nPerforming Kruskal-Wallis tests with Dunn's post-hoc...")
  
  stat_results <- lapply(analyzed_genes, function(gene) {
    gene_counts <- as.numeric(rpkm_matrix[gene, ])
    test_result <- perform_kw_dunn_test(gene_counts, groups)
    
    return(data.frame(
      Gene = gene,
      KW_pvalue = test_result$kw_pvalue,
      Dunn_24vsMock_padj = test_result$dunn_24_mock_padj,
      Dunn_48vsMock_padj = test_result$dunn_48_mock_padj,
      Sig_24hr_vs_Mock = test_result$sig_24_vs_mock,
      Sig_48hr_vs_Mock = test_result$sig_48_vs_mock,
      stringsAsFactors = FALSE
    ))
  })
  
  stat_results_df <- do.call(rbind, stat_results)
  
  # 6. Extract sample-level RPKM values
  sample_names <- colnames(rpkm_matrix)
  mock_samples <- grep("Mock", sample_names, value = TRUE)
  h5n1_24_samples <- grep("24h", sample_names, value = TRUE)
  h5n1_48_samples <- grep("48h", sample_names, value = TRUE)
  
  rpkm_sample_df <- as.data.frame(rpkm_matrix[analyzed_genes, ])
  rpkm_sample_df$Gene <- rownames(rpkm_sample_df)
  
  # 7. Calculate group averages
  rpkm_sample_df$Mock_avg <- rowMeans(rpkm_matrix[analyzed_genes, mock_samples, drop = FALSE])
  rpkm_sample_df$H5N1_24hr_avg <- rowMeans(rpkm_matrix[analyzed_genes, h5n1_24_samples, drop = FALSE])
  rpkm_sample_df$H5N1_48hr_avg <- rowMeans(rpkm_matrix[analyzed_genes, h5n1_48_samples, drop = FALSE])
  
  # 8. Merge with edgeR results
  # Prepare 24hr results
  results_24hr_select <- results_24hr %>%
    dplyr::select(logFC, PValue, FDR) %>%
    rename(logFC_24hr = logFC,
           pvalue_24hr = PValue,
           padj_24hr = FDR)
  results_24hr_select$Gene <- rownames(results_24hr)
  
  # Prepare 48hr results
  results_48hr_select <- results_48hr %>%
    dplyr::select(logFC, PValue, FDR) %>%
    rename(logFC_48hr = logFC,
           pvalue_48hr = PValue,
           padj_48hr = FDR)
  results_48hr_select$Gene <- rownames(results_48hr)
  
  # 9. Merge all data
  final_table <- rpkm_sample_df %>%
    left_join(stat_results_df, by = "Gene") %>%
    left_join(results_24hr_select, by = "Gene") %>%
    left_join(results_48hr_select, by = "Gene")
  
  # 10. Add gene annotation from dgelist
  gene_info <- dgelist$genes %>%
    dplyr::select(SYMBOL_UNIQUE, ENSEMBL, CHR, START, END, Length) %>%
    rename(Gene = SYMBOL_UNIQUE,
           Gene_ID = ENSEMBL,
           gene_chr = CHR,
           gene_start = START,
           gene_end = END,
           gene_length = Length)
  
  final_table <- final_table %>%
    left_join(gene_info, by = "Gene")
  
  # 11. Add COSMIC annotations
  cosmic_annotations <- cosmic_data %>%
    dplyr::select(any_of(c(cosmic_gene_col, "Role.in.Cancer", "Tumour.Types.Somatic.", 
                    "Tumour.Types..Germline.", "Molecular.Genetics"))) %>%
    rename_with(~"Gene", all_of(cosmic_gene_col))
  
  # Column names
  if ("Role.in.Cancer" %in% names(cosmic_annotations)) {
    cosmic_annotations <- cosmic_annotations %>%
      rename(COSMIC_Role = Role.in.Cancer)
  }
  
  if ("Tumour.Types.Somatic." %in% names(cosmic_annotations)) {
    cosmic_annotations <- cosmic_annotations %>%
      rename(tumour_types = Tumour.Types.Somatic.)
  } else if ("Tumour.Types..Somatic." %in% names(cosmic_annotations)) {
    cosmic_annotations <- cosmic_annotations %>%
      rename(tumour_types = Tumour.Types..Somatic.)
  }
  
  final_table <- final_table %>%
    left_join(cosmic_annotations, by = "Gene")
  
  # 12. Reorder columns for final output
  # Individual sample columns
  sample_cols <- setdiff(colnames(final_table), 
                        c("Gene", "Gene_ID", "gene_chr", "gene_start", "gene_end", 
                          "gene_length", "Mock_avg", "H5N1_24hr_avg", "H5N1_48hr_avg",
                          "logFC_24hr", "pvalue_24hr", "padj_24hr",
                          "logFC_48hr", "pvalue_48hr", "padj_48hr",
                          "KW_pvalue", "Dunn_24vsMock_padj", "Dunn_48vsMock_padj",
                          "Sig_24hr_vs_Mock", "Sig_48hr_vs_Mock",
                          "tumour_types", "COSMIC_Role", "Molecular.Genetics"))
  
  final_column_order <- c(
    "Gene_ID",
    "Gene",
    "tumour_types",
    sample_cols,  # All individual sample RPKM values
    "Mock_avg",
    "H5N1_24hr_avg", 
    "H5N1_48hr_avg",
    "logFC_24hr",
    "pvalue_24hr",
    "padj_24hr",
    "logFC_48hr",
    "pvalue_48hr", 
    "padj_48hr",
    "KW_pvalue",
    "Dunn_24vsMock_padj",
    "Dunn_48vsMock_padj",
    "Sig_24hr_vs_Mock",
    "Sig_48hr_vs_Mock",
    "gene_chr",
    "gene_start",
    "gene_end",
    "gene_length",
    "COSMIC_Role",
    "Molecular.Genetics"
  )
  
  # Keep only columns that exist
  final_column_order <- final_column_order[final_column_order %in% names(final_table)]
  
  final_table <- final_table %>%
    dplyr::select(all_of(final_column_order))
  
  # 13. Sort by significance
  final_table <- final_table %>%
    arrange(KW_pvalue, padj_24hr, padj_48hr)
  
  # 14. Write output as TSV
  write.table(final_table, 
              file = output_file, 
              sep = "\t",
              row.names = FALSE, 
              quote = FALSE)
  message("\nResults written to: ", output_file)
  
  # 15. Print summary statistics
  message("\n========== Analysis Summary ==========")
  message("Total cancer genes analyzed: ", nrow(final_table))
  message("Significant by Kruskal-Wallis (p < 0.05): ", 
          sum(final_table$KW_pvalue < 0.05, na.rm = TRUE))
  message("Significant 24hr vs Mock (Dunn's test): ", 
          sum(final_table$Sig_24hr_vs_Mock, na.rm = TRUE))
  message("Significant 48hr vs Mock (Dunn's test): ", 
          sum(final_table$Sig_48hr_vs_Mock, na.rm = TRUE))
  message("Significant in both timepoints: ", 
          sum(final_table$Sig_24hr_vs_Mock & final_table$Sig_48hr_vs_Mock, na.rm = TRUE))
  message("======================================\n")
  
  return(final_table)
}
