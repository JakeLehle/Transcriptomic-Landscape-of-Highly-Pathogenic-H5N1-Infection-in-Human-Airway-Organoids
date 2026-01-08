# KEGG Pathway Gene Table

get_kegg_pathway_genes_keggrest_only <- function(kegg_pathway_id, max_retries = 5) {
  require(KEGGREST)
  require(dplyr)

  cat("Retrieving genes for pathway:", kegg_pathway_id, "using KEGGREST API\n")

  for (attempt in 1:max_retries) {
    cat("KEGGREST attempt", attempt, "of", max_retries, "...\n")

    result <- try({
      # Get pathway info using KEGGREST
      pathway_info <- KEGGREST::keggGet(kegg_pathway_id)

      if (length(pathway_info) == 0) {
        stop("Pathway not found: ", kegg_pathway_id)
      }

      # Extract genes from pathway information
      pathway_genes <- pathway_info[[1]]$GENE

      if (is.null(pathway_genes) || length(pathway_genes) == 0) {
        stop("No genes found in pathway: ", kegg_pathway_id)
      }

      # Extract gene symbols - format is "gene_symbol; gene_description"
      gene_symbols <- sapply(strsplit(pathway_genes, ";"), function(x) trimws(x[1]))

      # Remove any empty strings and duplicates
      gene_symbols <- unique(gene_symbols[gene_symbols != "" & !is.na(gene_symbols)])

      if (length(gene_symbols) == 0) {
        stop("No valid gene symbols extracted")
      }

      cat("Successfully retrieved", length(gene_symbols), "genes via KEGGREST\n")
      return(gene_symbols)

    }, silent = TRUE)

    if (!inherits(result, "try-error") && length(result) > 0) {
      return(result)
    }

    # Wait before retry
    if (attempt < max_retries) {
      wait_time <- attempt * 2  # Increasing wait: 2, 4, 6, 8 seconds
      cat("Retrying in", wait_time, "seconds...\n")
      Sys.sleep(wait_time)
    }
  }

  # If all retries fail
  stop("Failed to retrieve genes for pathway: ", kegg_pathway_id, " after ", max_retries, " attempts")
}

# Get Pathway Name using KEGGREST
get_kegg_pathway_name <- function(kegg_pathway_id) {
  require(KEGGREST)

  pathway_name <- tryCatch({
    pathway_info <- KEGGREST::keggGet(kegg_pathway_id)
    if (length(pathway_info) > 0) {
      name <- pathway_info[[1]]$NAME
      if (!is.null(name)) {
        # Clean up the name - remove " - Homo sapiens (human)" part
        name <- gsub(" - .*", "", name)
        return(name)
      }
    }
    return(kegg_pathway_id)  # Return ID if name not found
  }, error = function(e) {
    return(kegg_pathway_id)  # Return ID on error
  })

  return(pathway_name)
}

# Create Color Gradient Function for logFC Values
create_logFC_color_gradient <- function(values, 
                                       palette = c("#2166AC", "#F7F7F7", "#B2182B")) {
  # Handle NA values
  valid_values <- values[!is.na(values)]
  if (length(valid_values) == 0) {
    return(rep("#FFFFFF", length(values)))
  }
  
  # Use fixed range for consistent color mapping across all pathways
  # This ensures -5 to +5 shows full saturation
  value_range <- c(-5, 5)
  max_abs <- max(abs(valid_values), na.rm = TRUE)
  
  # If values exceed our fixed range, extend the range
  if (max_abs > 5) {
    value_range <- c(-ceiling(max_abs), ceiling(max_abs))
  }
  
  # Create color mapping function with more colors for smoother gradient
  color_fn <- colorRampPalette(palette)
  color_scale <- color_fn(201)  # 201 colors for smoother gradient
  
  # Map values to colors
  colors <- sapply(values, function(x) {
    if (is.na(x)) return("#FFFFFF")
    
    # Scale value to 0-200 index
    scaled <- round((x - value_range[1]) / (value_range[2] - value_range[1]) * 200)
    scaled <- max(0, min(200, scaled))  # Ensure within bounds
    
    return(color_scale[scaled + 1])
  })
  
  return(colors)
}

# Draw Improved Color Legend for logFC Gradient
draw_color_legend <- function(values, color_palette = c("#2166AC", "#F7F7F7", "#B2182B")) {
  if (length(values) == 0 || all(is.na(values))) {
    return()
  }
  
  # Calculate range using the same logic as color function
  max_abs <- max(abs(values), na.rm = TRUE)
  value_range <- c(-5, 5)
  if (max_abs > 5) {
    value_range <- c(-ceiling(max_abs), ceiling(max_abs))
  }
  
  # Create gradient colors
  color_fn <- colorRampPalette(color_palette)
  gradient_colors <- color_fn(201)
  
  # Legend dimensions
  gradient_width <- 0.6
  gradient_height <- 0.1
  legend_y <- -10
  
  # Draw gradient rectangle
  x_positions <- seq(0.2, 0.2 + gradient_width, length.out = 201)
  
  for (i in 1:200) {
    grid.rect(x = x_positions[i], y = legend_y, 
              width = gradient_width/200, height = gradient_height,
              just = c("left", "center"),
              gp = gpar(fill = gradient_colors[i], col = NA))
  }
  
  # Add border around gradient
  grid.rect(x = 0.2, y = legend_y, 
            width = gradient_width, height = gradient_height,
            just = c("left", "center"),
            gp = gpar(fill = NA, col = "black"))
  
  # Add legend title - centered
  grid.text("log2 Fold Change", x = 0.5, y = legend_y + 0.15,
            gp = gpar(fontsize = 10, fontface = "bold"),
            hjust = 0.5, vjust = 0.5)
  
  # Add scale labels - centered below gradient
  grid.text(sprintf("%.1f", value_range[1]), x = 0.2, y = legend_y - 0.12,
            gp = gpar(fontsize = 9),
            hjust = 0.5, vjust = 0.5)

  grid.text("0.0", x = 0.5, y = legend_y - 0.12,
            gp = gpar(fontsize = 9),
            hjust = 0.5, vjust = 0.5)

  grid.text(sprintf("%.1f", value_range[2]), x = 0.8, y = legend_y - 0.12,
            gp = gpar(fontsize = 9),
            hjust = 0.5, vjust = 0.5)
  
  # Add direction labels
  grid.text("Downregulated", x = 0.35, y = legend_y - 0.2,
            gp = gpar(fontsize = 8, col = "blue", fontface = "bold"),
            hjust = 0.5, vjust = 0.5)

  grid.text("Upregulated", x = 0.65, y = legend_y - 0.2,
            gp = gpar(fontsize = 8, col = "red", fontface = "bold"),
            hjust = 0.5, vjust = 0.5)
}

# Generate Colored PDF Table
generate_colored_pdf_table <- function(table_data, pathway_id, pathway_name, output_dir, logFC_threshold) {
  
  pdf_file <- file.path(output_dir, paste0("Pathway_Table_Colored_", pathway_id, ".pdf"))
  
  pdf_table <- table_data %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::select(Gene, logFC_24hr, logFC_48hr) %>%
    rename(`24hr_vs_Mock` = logFC_24hr,
           `48hr_vs_Mock` = logFC_48hr)
  
  # Calculate dynamic height with more space for legend
  table_height <- max(18, nrow(pdf_table) * 0.25 + 4)
  pdf(pdf_file, width = 12, height = table_height)
  
  tryCatch({
    grid.newpage()
    
    # Layout
    table_width <- 0.3
    cell_height <- 0.02
    cell_width <- 0.1  # Reduced column width
    
    # Centering Calc
    start_x <- 0.4
    start_y <- 0.95
    
    # Title - Centered
    grid.text(paste(pathway_name, "(", pathway_id, ")"),
              x = 0.5, y = 0.99, 
              gp = gpar(fontsize = 16, fontface = "bold"),
              hjust = 0.5, vjust = 0.5)
    
    grid.text(paste0("Genes with |logFC| > ", logFC_threshold, " in at least one timepoint | Total genes:", nrow(pdf_table)),
              x = 0.5, y = 0.97, 
              gp = gpar(fontsize = 12),
              hjust = 0.5, vjust = 0.5)
    
    if (nrow(pdf_table) > 0 && nrow(pdf_table) <= 50) {
      # Column headers
      headers <- c("Gene", "24hr vs Mock", "48hr vs Mock")
      
      # Header background
      for (col in 1:3) {
        grid.rect(x = start_x + (col-1)*cell_width, y = start_y,
                  width = cell_width, height = cell_height,
                  gp = gpar(fill = "lightgray", col = "black"))
        
        # x and y centering
        grid.text(headers[col], 
                  x = start_x + (col-1)*cell_width + cell_width/100,
                  y = start_y + cell_height/100,
                  gp = gpar(fontsize = 11, fontface = "bold"),
                  hjust = 0.5, vjust = 0.5)
      }
      
      # Get color ranges
      all_values <- c(pdf_table$`24hr_vs_Mock`, pdf_table$`48hr_vs_Mock`)
      color_palette <- c("#2166AC", "#F7F7F7", "#B2182B")
      
      # Data rows
      for (row in 1:nrow(pdf_table)) {
        y_pos <- start_y - row * cell_height
        
        # Gene name cell
        grid.rect(x = start_x, y = y_pos,
                  width = cell_width, height = cell_height,
                  gp = gpar(fill = "white", col = "black"))
        
        grid.text(pdf_table$Gene[row],
                  x = start_x + cell_width/100, 
                  y = y_pos + cell_height/100,
                  gp = gpar(fontsize = 12, fontface = "bold"),  # Increased font
                  hjust = 0.5, vjust = 0.5)
        
        # 24hr logFC cell
        color_24hr <- create_logFC_color_gradient(pdf_table$`24hr_vs_Mock`[row], color_palette)
        grid.rect(x = start_x + cell_width, y = y_pos,
                  width = cell_width, height = cell_height,
                  gp = gpar(fill = color_24hr, col = "black"))
        
        grid.text(sprintf("%.2f", pdf_table$`24hr_vs_Mock`[row]),
                  x = start_x + cell_width + cell_width/100, 
                  y = y_pos + cell_height/100,
                  gp = gpar(fontsize = 14, fontface = "bold"),
                  just = c("center", "center"))  # Fixed centering
        
        # 48hr logFC cell
        color_48hr <- create_logFC_color_gradient(pdf_table$`48hr_vs_Mock`[row], color_palette)
        grid.rect(x = start_x + 2*cell_width, y = y_pos,
                  width = cell_width, height = cell_height,
                  gp = gpar(fill = color_48hr, col = "black"))
        
        grid.text(sprintf("%.2f", pdf_table$`48hr_vs_Mock`[row]),
                  x = start_x + 2*cell_width + cell_width/40,
                  y = y_pos + cell_height/40,
                  gp = gpar(fontsize = 14, fontface = "bold"),
                  hjust = 0.5, vjust = 0.5)
      }
      
      # Ledgend position
      legend_y <- 0.02  # Further down to avoid overlap
      pushViewport(viewport(x = 0.5, y = legend_y, width = 0.7, height = 0.15, 
                           just = c("center", "bottom")))
      draw_color_legend(all_values, color_palette)
      popViewport()
      
    } else if (nrow(pdf_table) > 50) {
      grid.text("Table too large for colored display (>50 genes)\nSee CSV file for full data",
                x = 0.5, y = 0.5, 
                gp = gpar(fontsize = 12, col = "red"),
                hjust = 0.5, vjust = 0.5)
    }
    
  }, error = function(e) {
    grid.text(paste("Error generating colored table:", e$message), 
              x = 0.5, y = 0.5, 
              gp = gpar(fontsize = 10, col = "red"),
              hjust = 0.5, vjust = 0.5)
  })
  
  dev.off()
  return(pdf_file)
}


# Pathway Comparison table
create_pathway_comparison_table_keggrest_only <- function(deg_24hr, deg_48hr, kegg_pathway_id,
                                                         logFC_threshold = 1.5, output_dir = getwd()) {

  # Load required libraries
  require(KEGGREST)
  require(dplyr)
  require(tibble)
  require(ggplot2)
  require(gridExtra)
  require(scales)
  require(grid)

  # Validate inputs
  if (!all(c("logFC", "FDR") %in% colnames(deg_24hr))) {
    stop("deg_24hr must contain 'logFC' and 'FDR' columns")
  }
  if (!all(c("logFC", "FDR") %in% colnames(deg_48hr))) {
    stop("deg_48hr must contain 'logFC' and 'FDR' columns")
  }

  cat("=== Starting Pathway Analysis ===\n")
  pathway_name <- get_kegg_pathway_name(kegg_pathway_id)
  cat("KEGG Pathway:", pathway_name, "(", kegg_pathway_id, ")\n")
  cat("LogFC Threshold:", logFC_threshold, "\n")

  # Step 1: Get KEGG pathway genes
  pathway_genes <- get_kegg_pathway_genes_keggrest_only(kegg_pathway_id)

  cat("Total genes in pathway:", length(pathway_genes), "\n")

  if (length(pathway_genes) == 0) {
    stop("No genes found for pathway: ", kegg_pathway_id)
  }

  # Step 2: Prepare comparative data
  comparative_data <- data.frame(
    GeneSymbol = rownames(deg_24hr),
    logFC_24hr = deg_24hr$logFC,
    logFC_48hr = deg_48hr$logFC,
    FDR_24hr = deg_24hr$FDR,
    FDR_48hr = deg_48hr$FDR,
    stringsAsFactors = FALSE
  ) %>%
    filter(GeneSymbol %in% pathway_genes) %>%
    mutate(
      avg_logFC = (logFC_24hr + logFC_48hr) / 2,
      max_abs_logFC = pmax(abs(logFC_24hr), abs(logFC_48hr))
    )

  cat("Genes from pathway found in dataset:", nrow(comparative_data), "\n")

  # Step 3: Apply filtering criteria
  filtered_data <- comparative_data %>%
    filter(max_abs_logFC > logFC_threshold) %>%
    filter(!(abs(logFC_24hr) < logFC_threshold & abs(logFC_48hr) < logFC_threshold))

  cat("Genes after filtering:", nrow(filtered_data), "\n")

  if (nrow(filtered_data) == 0) {
    warning("No genes passed the filtering criteria for pathway: ", kegg_pathway_id)
    return(list(
      summary = list(
        pathway_id = kegg_pathway_id,
        pathway_name = pathway_name,
        total_pathway_genes = length(pathway_genes),
        filtered_genes = 0,
        message = "No genes passed filtering criteria"
      ),
      data = NULL
    ))
  }

  # Step 4: Create final table for output
  final_table <- filtered_data %>%
    dplyr::select(GeneSymbol, logFC_24hr, logFC_48hr, avg_logFC) %>%
    arrange(desc(avg_logFC)) %>%
    column_to_rownames("GeneSymbol")

  # Step 5: Save CSV file
  csv_file <- file.path(output_dir, paste0("Pathway_Data_", kegg_pathway_id, ".csv"))
  write.csv(final_table, csv_file, row.names = TRUE)

  # Step 6: Generate enhanced PDF with color coding
  pdf_file <- generate_colored_pdf_table(final_table, kegg_pathway_id, pathway_name, output_dir, logFC_threshold)

  # Step 7: Create comprehensive summary statistics
  summary_stats <- list(
    pathway_id = kegg_pathway_id,
    pathway_name = pathway_name,
    total_pathway_genes = length(pathway_genes),
    genes_in_dataset = nrow(comparative_data),
    filtered_genes = nrow(final_table),
    up_regulated_24hr = sum(final_table$logFC_24hr > logFC_threshold),
    down_regulated_24hr = sum(final_table$logFC_24hr < -logFC_threshold),
    up_regulated_48hr = sum(final_table$logFC_48hr > logFC_threshold),
    down_regulated_48hr = sum(final_table$logFC_48hr < -logFC_threshold),
    csv_file = csv_file,
    pdf_file = pdf_file
  )

  # Print enhanced summary
  cat("\n=== Analysis Summary ===\n")
  cat("Pathway:", summary_stats$pathway_name, "(", summary_stats$pathway_id, ")\n")
  cat("Total pathway genes:", summary_stats$total_pathway_genes, "\n")
  cat("Genes found in dataset:", summary_stats$genes_in_dataset, "\n")
  cat("Genes passing filter:", summary_stats$filtered_genes, "\n")
  cat("Up-regulated in 24hr:", summary_stats$up_regulated_24hr, "\n")
  cat("Down-regulated in 24hr:", summary_stats$down_regulated_24hr, "\n")
  cat("Up-regulated in 48hr:", summary_stats$up_regulated_48hr, "\n")
  cat("Down-regulated in 48hr:", summary_stats$down_regulated_48hr, "\n")
  cat("Output files saved to:", output_dir, "\n")

  # Return comprehensive results
  results <- list(
    summary = summary_stats,
    data = final_table,
    pathway_genes = pathway_genes,
    pathway_name = pathway_name
  )

  return(results)
}



# Main Batch Analysis Function
batch_pathway_analysis_keggrest_only <- function(deg_24hr, deg_48hr, pathway_list,
                                                logFC_threshold = 1.5, output_dir = getwd()) {

  # Check if KEGGREST is available
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("KEGGREST package is required. Install with: BiocManager::install('KEGGREST')")
  }

  results <- list()
  successful_pathways <- character()
  failed_pathways <- character()

  cat("Starting batch pathway analysis using KEGGREST API only\n")
  cat("Number of pathways to process:", length(pathway_list), "\n")
  cat("Maximum retries per pathway: 5\n")

  for (pathway_id in pathway_list) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("Processing pathway:", pathway_id, "\n")

    pathway_result <- tryCatch({
      create_pathway_comparison_table_keggrest_only(
        deg_24hr = deg_24hr,
        deg_48hr = deg_48hr,
        kegg_pathway_id = pathway_id,
        logFC_threshold = logFC_threshold,
        output_dir = output_dir
      )
    }, error = function(e) {
      cat("ERROR processing", pathway_id, ":", e$message, "\n")
      failed_pathways <<- c(failed_pathways, pathway_id)
      NULL
    })

    if (!is.null(pathway_result) && !is.null(pathway_result$data)) {
      results[[pathway_id]] = pathway_result
      successful_pathways = c(successful_pathways, pathway_id)
      cat("Successfully processed:", pathway_id, "\n")
    } else {
      failed_pathways = c(failed_pathways, pathway_id)
      cat("Failed to process:", pathway_id, "\n")
    }
  }

  # Generate comprehensive summary report
  generate_summary_report(results, successful_pathways, failed_pathways, output_dir)

  return(results)
}

# Generate Summary Report
generate_summary_report <- function(results, successful_pathways, failed_pathways, output_dir) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("BATCH ANALYSIS SUMMARY REPORT\n")
  cat("Method: KEGGREST API Only (5 retries per pathway)\n")
  cat(rep("=", 70), "\n", sep = "")

  total_pathways <- length(successful_pathways) + length(failed_pathways)
  success_rate <- round(length(successful_pathways) / total_pathways * 100, 1)

  cat("Successful pathways:", length(successful_pathways), "/", total_pathways,
      paste0("(", success_rate, "%)"), "\n\n")

  if (length(successful_pathways) > 0) {
    cat("SUCCESSFUL PATHWAYS:\n")
    for (pathway in successful_pathways) {
      summary <- results[[pathway]]$summary
      cat(sprintf("  %-12s: %-40s | Genes: %3d | Filtered: %3d\n",
                 pathway, substr(summary$pathway_name, 1, 40),
                 summary$genes_in_dataset, summary$filtered_genes))
    }
    cat("\n")
  }

  if (length(failed_pathways) > 0) {
    cat("FAILED PATHWAYS:", length(failed_pathways), "\n")
    cat("  ", paste(failed_pathways, collapse = ", "), "\n\n")
    cat("Note: Failed pathways may be due to KEGG API limitations or network issues.\n")
    cat("Consider running failed pathways again later.\n")
  }

  # Save summary to file
  summary_file <- file.path(output_dir, "Pathway_Analysis_Summary_KEGGREST.txt")
  sink(summary_file)
  cat("Pathway Analysis Summary - KEGGREST API Only\n")
  cat("Generated:", date(), "\n")
  cat("Total pathways processed:", total_pathways, "\n")
  cat("Successful pathways:", length(successful_pathways), "\n")
  cat("Failed pathways:", length(failed_pathways), "\n\n")

  if (length(successful_pathways) > 0) {
    cat("DETAILED RESULTS:\n")
    for (pathway in successful_pathways) {
      summary <- results[[pathway]]$summary
      cat(sprintf("%s: %s\n", pathway, summary$pathway_name))
      cat(sprintf("  Total genes: %d | In dataset: %d | Filtered: %d\n",
                 summary$total_pathway_genes, summary$genes_in_dataset,
                 summary$filtered_genes))
      cat(sprintf("  Up 24hr: %d | Down 24hr: %d | Up 48hr: %d | Down 48hr: %d\n\n",
                 summary$up_regulated_24hr, summary$down_regulated_24hr,
                 summary$up_regulated_48hr, summary$down_regulated_48hr))
    }
  }
  sink()

  cat("Detailed summary saved to:", summary_file, "\n")
}

# Test function with color functionality
test_keggrest_with_colors <- function() {
  cat("Testing KEGGREST with color functionality...\n")

  # Test with a known pathway
  test_pathway <- "hsa04060"
  cat("Testing pathway:", test_pathway, "\n")

  genes <- get_kegg_pathway_genes_keggrest_only(test_pathway)
  cat("Retrieved", length(genes), "genes\n")

  name <- get_kegg_pathway_name(test_pathway)
  cat("Pathway name:", name, "\n")

  # Test color gradient
  test_values <- c(-5, -2.5, -1, 0, 1, 2.5, 5)
  colors <- create_logFC_color_gradient(test_values)
  cat("Color gradient test successful\n")
  
  cat("KEGGREST with color functionality working correctly!\n")
}

# Initialize and test
cat("KEGG Pathway Analysis Pipeline\n")
test_keggrest_with_colors()
