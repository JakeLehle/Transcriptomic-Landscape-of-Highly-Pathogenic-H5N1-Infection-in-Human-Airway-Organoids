# Heatmap for KEGG pathway genes

# Function to create heatmap for genes in a specific KEGG pathway
create_kegg_pathway_heatmap <- function(pathway_id = "hsa04060", # Cytokine-cytokine receptor interaction
                                        pathway_name = NULL,
                                        expression_matrix = NULL,
                                        y_all_object = NULL,
                                        sample_groups = NULL,
                                        output_filename = NULL) {
  
  # If expression_matrix is not provided, create it from y_all_object
  if (is.null(expression_matrix)) {
    if (is.null(y_all_object)) {
      stop("Either expression_matrix or y_all_object must be provided")
    }
    expression_matrix <- cpm(y_all_object, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  }
  
  # If sample_groups is not provided, try to get it from y_all_object
  if (is.null(sample_groups)) {
    if (!is.null(y_all_object)) {
      sample_groups <- y_all_object$samples$group
    } else {
      stop("Sample groups must be provided either directly or through y_all_object")
    }
  }
  
  # If pathway_name is not provided, try to get it from KEGG
  if (is.null(pathway_name)) {
    # Get pathway info using KEGGREST
    pathway_info <- KEGGREST::keggGet(pathway_id)
    pathway_name <- pathway_info[[1]]$NAME
    if (is.null(pathway_name)) {
      pathway_name <- pathway_id
    } else {
      pathway_name <- gsub(" - .*", "", pathway_name)  # Clean up the name
    }
  }
  
  # Get genes in the KEGG pathway using KEGGREST
  pathway_info <- KEGGREST::keggGet(pathway_id)
  pathway_genes <- pathway_info[[1]]$GENE
  
  # Extract gene IDs - they're in the format "gene_symbol; gene_name"
  gene_ids <- names(pathway_genes)
  gene_symbols <- sapply(strsplit(pathway_genes, ";"), function(x) x[1])
  
  # Filter out genes not in our expression matrix
  genes_in_matrix <- gene_symbols[gene_symbols %in% rownames(expression_matrix)]
  
  if (length(genes_in_matrix) == 0) {
    warning(paste("No genes from pathway", pathway_id, "found in expression matrix"))
    return(NULL)
  }
  
  # Subset the expression matrix for these genes
  pathway_expression <- expression_matrix[genes_in_matrix, ]
  
  # Create annotation for samples
  sample_annotation <- data.frame(Group = sample_groups)
  rownames(sample_annotation) <- colnames(expression_matrix)
  
  # Create color palette for groups
  unique_groups <- unique(sample_groups)
  group_colors <- brewer.pal(length(unique_groups), "Set1")[1:length(unique_groups)]
  names(group_colors) <- unique_groups
  annotation_colors <- list(Group = group_colors)
  
  # Set output filename if not provided
  if (is.null(output_filename)) {
    output_filename <- paste0("Heatmap_", pathway_id, ".pdf")
  }
  
  # Create the heatmap
  pheatmap(pathway_expression,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           border_color = NA,
           scale = "row",  # Scale by row (gene) to show relative expression
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           show_rownames = TRUE,
           show_colnames = TRUE,
	   display_numbers = round(pathway_expression, 1),
           annotation_col = sample_annotation,
           annotation_colors = annotation_colors,
           main = paste(pathway_name),
           fontsize_row = 10,
           fontsize_col = 12,
           filename = output_filename,
           width = 10,
           height = max(6, length(genes_in_matrix) * 0.2 + 4))
  
  cat(paste("Heatmap saved as:", output_filename, "\n"))
  
  # Return the list of genes in the pathway that were plotted
  return(list(genes = genes_in_matrix, pathway_name = pathway_name))
}

