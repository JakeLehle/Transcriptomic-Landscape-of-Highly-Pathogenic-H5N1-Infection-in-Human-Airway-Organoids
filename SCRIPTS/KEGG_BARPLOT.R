#' KEGG Pathway Bar Plot Generator
#' 
#' Generates side-by-side bar plots for KEGG pathway genes showing logFC values
#' at 24hr and 48hr timepoints with color-coded bars

#' Create Bar Plot for Pathway Comparison
#' 
#' @param table_data Data frame with genes as rownames and logFC columns
#' @param pathway_id KEGG pathway ID
#' @param pathway_name Full pathway name
#' @param output_dir Directory to save the plot
#' @param logFC_threshold Threshold used for filtering (for subtitle)
#' @return Path to saved PDF file
generate_pathway_barplot <- function(table_data, pathway_id, pathway_name, 
                                    output_dir, logFC_threshold = 3) {
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  
  # Prepare data for plotting
  plot_data <- table_data %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::select(Gene, logFC_24hr, logFC_48hr) %>%
    tidyr::pivot_longer(cols = c(logFC_24hr, logFC_48hr),
                       names_to = "Timepoint",
                       values_to = "logFC") %>%
    mutate(Timepoint = factor(Timepoint, 
                             levels = c("logFC_24hr", "logFC_48hr"),
                             labels = c("24hr vs Mock", "48hr vs Mock")))
  
  # Sort genes by average logFC (descending)
  gene_order <- table_data %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(avg_logFC = (logFC_24hr + logFC_48hr) / 2) %>%
    arrange(desc(avg_logFC)) %>%
    pull(Gene)
  
  plot_data$Gene <- factor(plot_data$Gene, levels = gene_order)
  
  # Create color mapping function
  color_scale <- scale_fill_gradient2(
    low = "#2166AC",      # Dark blue for downregulated
    mid = "#F7F7F7",      # White for neutral
    high = "#B2182B",     # Dark red for upregulated
    midpoint = 0,
    limits = c(-max(abs(plot_data$logFC)), max(abs(plot_data$logFC))),
    name = "log2 Fold Change"
  )
  
  # Calculate dynamic plot dimensions
  n_genes <- length(unique(plot_data$Gene))
  plot_width <- max(12, n_genes * 0.4)
  plot_height <- max(8, n_genes * 0.25)
  
  # Create the bar plot
  p <- ggplot(plot_data, aes(x = Gene, y = logFC, fill = logFC)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             width = 0.7, color = "black", linewidth = 0.3) +
    # Add logFC values as text labels on bars
    geom_text(aes(label = sprintf("%.2f", logFC)),
              hjust = ifelse(plot_data$logFC > 0, -0.1, 1.1),
              size = 3.5, fontface = "bold") +
    color_scale +
    facet_wrap(~ Timepoint, ncol = 1, scales = "free_y") +
    coord_flip() +  # Horizontal bars for better gene name readability
    labs(
      title = paste0(pathway_name, " (", pathway_id, ")"),
      subtitle = paste0("Genes with |logFC| > ", logFC_threshold, 
                       " in at least one timepoint | Total genes: ", n_genes),
      x = "Gene Symbol",
      y = "log2 Fold Change"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      axis.text.y = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 12, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    # Expand x-axis limits to accommodate text labels
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))
  
  # Save the plot
  pdf_file <- file.path(paste0(output_dir, "Pathway_Barplot_", pathway_id, ".pdf"))
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height, 
         units = "in", dpi = 300)
  
  cat("✓ Bar plot saved:", pdf_file, "\n")
  return(pdf_file)
}

#' Generate Alternative Side-by-Side Bar Plot (Both timepoints on same panel)
#' 
#' @param table_data Data frame with genes as rownames and logFC columns
#' @param pathway_id KEGG pathway ID
#' @param pathway_name Full pathway name
#' @param output_dir Directory to save the plot
#' @param logFC_threshold Threshold used for filtering (for subtitle)
#' @return Path to saved PDF file
generate_pathway_barplot_sidebyside <- function(table_data, pathway_id, pathway_name, 
                                               output_dir, logFC_threshold = 3) {
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  
  # Prepare data for plotting
  plot_data <- table_data %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::select(Gene, logFC_24hr, logFC_48hr) %>%
    tidyr::pivot_longer(cols = c(logFC_24hr, logFC_48hr),
                       names_to = "Timepoint",
                       values_to = "logFC") %>%
    mutate(Timepoint = factor(Timepoint, 
                             levels = c("logFC_24hr", "logFC_48hr"),
                             labels = c("24hr", "48hr")))
  
  # Sort genes by average logFC (descending)
  gene_order <- table_data %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(avg_logFC = (logFC_24hr + logFC_48hr) / 2) %>%
    arrange(desc(avg_logFC)) %>%
    pull(Gene)
  
  plot_data$Gene <- factor(plot_data$Gene, levels = gene_order)
  
  # Create color mapping function
  color_scale <- scale_fill_gradient2(
    low = "#2166AC",      # Dark blue for downregulated
    mid = "#F7F7F7",      # White for neutral
    high = "#B2182B",     # Dark red for upregulated
    midpoint = 0,
    limits = c(-max(abs(plot_data$logFC)), max(abs(plot_data$logFC))),
    name = "log2FC"
  )
  
  # Calculate dynamic plot dimensions
  n_genes <- length(unique(plot_data$Gene))
  plot_height <- max(8, n_genes * 0.3)
  plot_width <- 12
  
  # Create the side-by-side bar plot
  p <- ggplot(plot_data, aes(x = Gene, y = logFC, fill = logFC)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
             aes(group = Timepoint), width = 0.8, color = "black", linewidth = 0.3) +
    # Add logFC values as text labels on bars
    geom_text(aes(label = sprintf("%.2f", logFC), group = Timepoint),
              position = position_dodge(width = 0.9),
              hjust = ifelse(plot_data$logFC > 0, -0.1, 1.1),
              size = 3, fontface = "bold") +
    color_scale +
    coord_flip() +
    labs(
      title = paste0(pathway_name, " (", pathway_id, ")"),
      subtitle = paste0("Genes with |logFC| > ", logFC_threshold, 
                       " in at least one timepoint | Total genes: ", n_genes),
      x = "Gene Symbol",
      y = "log2 Fold Change"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      axis.text.y = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    # Add timepoint labels
    facet_grid(. ~ Timepoint) +
    # Expand x-axis limits to accommodate text labels
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))
  
  # Save the plot
  pdf_file <- file.path(output_dir, paste0("Pathway_Barplot_SideBySide_", pathway_id, ".pdf"))
  ggsave(pdf_file, plot = p, width = plot_width, height = plot_height, 
         units = "in", dpi = 300)
  
  cat("✓ Side-by-side bar plot saved:", pdf_file, "\n")
  return(pdf_file)
}

#' Create Pathway Bar Plot - Main Wrapper Function
#' 
#' This function integrates with your existing pipeline and generates bar plots
#' 
#' @param deg_24hr DEG results for 24hr timepoint
#' @param deg_48hr DEG results for 48hr timepoint
#' @param kegg_pathway_id KEGG pathway ID
#' @param logFC_threshold Threshold for filtering genes
#' @param output_dir Output directory for plots
#' @param plot_style Either "faceted" (separate panels) or "sidebyside" (grouped bars)
#' @return List with summary and plot file path
create_pathway_barplot <- function(deg_24hr, deg_48hr, kegg_pathway_id,
                                   logFC_threshold = 3, output_dir = getwd(),
                                   plot_style = "faceted") {
  
  require(KEGGREST)
  require(dplyr)
  
  cat("\n=== Creating Bar Plot for Pathway ===\n")
  
  # Get pathway name
  pathway_name <- get_kegg_pathway_name(kegg_pathway_id)
  cat("Pathway:", pathway_name, "(", kegg_pathway_id, ")\n")
  
  # Get pathway genes
  pathway_genes <- get_kegg_pathway_genes_keggrest_only(kegg_pathway_id)
  
  # Prepare data (same as table function)
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
  
  # Apply filtering
  filtered_data <- comparative_data %>%
    filter(max_abs_logFC > logFC_threshold) %>%
    filter(!(abs(logFC_24hr) < logFC_threshold & abs(logFC_48hr) < logFC_threshold))
  
  cat("Genes passing filter:", nrow(filtered_data), "\n")
  
  if (nrow(filtered_data) == 0) {
    warning("No genes passed the filtering criteria for pathway: ", kegg_pathway_id)
    return(NULL)
  }
  
  # Create table data
  final_table <- filtered_data %>%
    dplyr::select(GeneSymbol, logFC_24hr, logFC_48hr, avg_logFC) %>%
    arrange(desc(avg_logFC)) %>%
    tibble::column_to_rownames("GeneSymbol")
  
  # Generate bar plot based on style
  if (plot_style == "sidebyside") {
    plot_file <- generate_pathway_barplot_sidebyside(
      final_table, kegg_pathway_id, pathway_name, output_dir, logFC_threshold
    )
  } else {
    plot_file <- generate_pathway_barplot(
      final_table, kegg_pathway_id, pathway_name, output_dir, logFC_threshold
    )
  }
  
  return(list(
    plot_file = plot_file,
    n_genes = nrow(final_table),
    pathway_name = pathway_name
  ))
}

#' Batch Bar Plot Generation
#' 
#' Generate bar plots for multiple pathways
#' 
#' @param deg_24hr DEG results for 24hr timepoint
#' @param deg_48hr DEG results for 48hr timepoint
#' @param pathway_list Vector of KEGG pathway IDs
#' @param logFC_threshold Threshold for filtering genes
#' @param output_dir Output directory for plots
#' @param plot_style Either "faceted" or "sidebyside"
#' @return List of results for each pathway
batch_pathway_barplots <- function(deg_24hr, deg_48hr, pathway_list,
                                   logFC_threshold = 3, output_dir = getwd(),
                                   plot_style = "faceted") {
  
  results <- list()
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("BATCH BAR PLOT GENERATION\n")
  cat(rep("=", 70), "\n", sep = "")
  cat("Number of pathways:", length(pathway_list), "\n")
  cat("Plot style:", plot_style, "\n\n")
  
  for (pathway_id in pathway_list) {
    cat("Processing:", pathway_id, "...\n")
    
    result <- tryCatch({
      create_pathway_barplot(
        deg_24hr = deg_24hr,
        deg_48hr = deg_48hr,
        kegg_pathway_id = pathway_id,
        logFC_threshold = logFC_threshold,
        output_dir = output_dir,
        plot_style = plot_style
      )
    }, error = function(e) {
      cat("❌ Error:", e$message, "\n")
      NULL
    })
    
    if (!is.null(result)) {
      results[[pathway_id]] <- result
      cat("✅ Completed:", pathway_id, "\n\n")
    }
  }
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("SUMMARY: Generated", length(results), "bar plots\n")
  cat(rep("=", 70), "\n", sep = "")
  
  return(results)
}

cat("✓ KEGG Pathway Bar Plot Generator loaded successfully\n")
cat("Available functions:\n")
cat("  - create_pathway_barplot(): Generate bar plot for single pathway\n")
cat("  - batch_pathway_barplots(): Generate bar plots for multiple pathways\n")
cat("Plot styles: 'faceted' (separate panels) or 'sidebyside' (grouped bars)\n")
