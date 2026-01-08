run_comprehensive_enrichment <- function(deg_results, 
                                        comparison_name, 
                                        y_all_object,
                                        pval_cutoff = 0.05, 
                                        logfc_cutoff = 1.5,
                                        fdr_cutoff = 0.05) {
  
  original_wd <- "/work/sdz852/WORKING/RNA-seq/H5N1/JAKE/aligned/"
  
  # Create a directory for this comparison
  dir_name <- paste0("Enrichment_Results_", comparison_name)
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  # Set working directory for output files
  setwd(dir_name)
  
  cat(paste("\n=== Running enrichment analysis for:", comparison_name, "===\n"))
  
  # Prepare the gene list for enrichment analysis
  significant_genes <- deg_results %>%
    as.data.frame() %>%
    filter(FDR < fdr_cutoff & abs(logFC) > logfc_cutoff) %>%
    pull(SYMBOL_UNIQUE)
  
  if (length(significant_genes) == 0) {
    cat("No significant genes found for enrichment analysis.\n")
    setwd(original_wd)
    return(NULL)
  }
  
  # Extract ENTREZIDs from y_all$genes for significant genes
  sig_gene_annotation <- y_all_object$genes[significant_genes, ]
  ensembl_ids <- sig_gene_annotation$ENSEMBL
  
  # Map ENSEMBL to ENTREZID
  entrez_mapping <- bitr(ensembl_ids, 
                         fromType = "ENSEMBL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  
  significant_entrez <- unique(entrez_mapping$ENTREZID)
  
  # Prepare background gene set (all genes in analysis)
  all_entrez <- bitr(y_all_object$genes$ENSEMBL,
                     fromType = "ENSEMBL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)$ENTREZID
  
  cat(paste("Number of significant genes:", length(significant_genes), "\n"))
  cat(paste("Number with ENTREZID mapping:", length(significant_entrez), "\n"))
  
  # Define theme for all plots
  large_text_theme <- theme(
    text = element_text(size = 30),           # Increase general text size by 100% (from ~8 to 16)
    axis.text = element_text(size = 24),      # Axis text size
    axis.text.y = element_text(size = 28),
    axis.text.x = element_text(size = 24),
    axis.title = element_text(size = 32),     # Axis title size
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),  # Plot title
    legend.text = element_text(size = 24),    # Legend text
    legend.title = element_text(size = 30)    # Legend title
  )
  
  # ---------------------------------------------------------------------------
  # 1. GENE ONTOLOGY (GO) ENRICHMENT
  # ---------------------------------------------------------------------------
  cat("Running GO enrichment...\n")
  
  go_ontologies <- c("BP", "CC", "MF")
  go_results <- list()
  
  for (ont in go_ontologies) {
    go_enrich <- enrichGO(gene          = significant_entrez,
                          universe      = all_entrez,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "ENTREZID",
                          ont           = ont,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = pval_cutoff,
                          qvalueCutoff  = 0.2,
                          readable      = TRUE)
    
    go_results[[ont]] <- go_enrich
    
    # Save results
    if (nrow(go_enrich) > 0) {
      filename <- paste0(comparison_name, "_GO_", ont, ".csv")
      write.csv(go_enrich@result, filename)
      
      # Create dot plot
      p <- dotplot(go_enrich, showCategory = 10, 
                   title = paste("GO", ont, "-", comparison_name)) +
        large_text_theme
      plot_filename <- paste0(comparison_name, "_GO_", ont, "_Dotplot.pdf")
      ggsave(plot_filename, plot = p, width = 12, height = 10, dpi = 300)  # Increased size to accommodate larger text
      
      # Create enrichment map
      if (nrow(go_enrich) >= 5) {
        p2 <- emapplot(pairwise_termsim(go_enrich), showCategory = 10, size_category = 3, node_label = "none", size_edge = 2) + 
          geom_cnet_label(node_label = "all", size = 13) +
          theme(legend.text = element_text(size = 24),
                legend.title = element_text(size = 30))
        map_filename <- paste0(comparison_name, "_GO_", ont, "_EnrichmentMap.pdf")
        ggsave(map_filename, plot = p2, width = 14, height = 12, dpi = 300)  # Increased size
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # 2. DISEASE ENRICHMENT
  # ---------------------------------------------------------------------------
  cat("Running disease enrichment analysis...\n")
  
  disease_enrich <- tryCatch({
    enrichDO(gene          = significant_entrez,
             pvalueCutoff  = pval_cutoff,
             pAdjustMethod = "BH",
             universe      = all_entrez,
             minGSSize     = 5,
             maxGSSize     = 500,
             qvalueCutoff  = 0.2,
             readable      = TRUE)
  }, error = function(e) {
    cat("DO enrichment failed, trying DisGeNET instead...\n")
    enrichDGN(gene          = significant_entrez,
              pvalueCutoff  = pval_cutoff,
              pAdjustMethod = "BH",
              universe      = all_entrez,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.2,
              readable      = TRUE)
  })
  
  if (!is.null(disease_enrich) && nrow(disease_enrich) > 0) {
    filename <- paste0(comparison_name, "_Disease_Enrichment.csv")
    write.csv(disease_enrich@result, filename)
    
    p <- dotplot(disease_enrich, showCategory = 10, 
                 title = paste("Disease Enrichment -", comparison_name)) +
      large_text_theme
    plot_filename <- paste0(comparison_name, "_Disease_Dotplot.pdf")
    ggsave(plot_filename, plot = p, width = 12, height = 10, dpi = 300)
  }
  
  # ---------------------------------------------------------------------------
  # 3. KEGG PATHWAY ENRICHMENT
  # ---------------------------------------------------------------------------
  cat("Running KEGG pathway enrichment...\n")
  
  kegg_enrich <- enrichKEGG(gene         = significant_entrez,
                            organism     = "hsa",
                            keyType      = "kegg",
                            pvalueCutoff = pval_cutoff,
                            pAdjustMethod = "BH",
                            universe     = all_entrez)
  
  if (!is.null(kegg_enrich) && nrow(kegg_enrich) > 0) {
    kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    filename <- paste0(comparison_name, "_KEGG_Pathways.csv")
    write.csv(kegg_enrich@result, filename)
    
    p <- dotplot(kegg_enrich, showCategory = 10, 
                 title = paste("KEGG Pathways -", comparison_name)) +
      large_text_theme
    plot_filename <- paste0(comparison_name, "_KEGG_Dotplot.pdf")
    ggsave(plot_filename, plot = p, width = 12, height = 10, dpi = 300)
    
    # Save top pathway visualization
    if (nrow(kegg_enrich) > 0) {
      top_pathway <- kegg_enrich@result$ID[1]
      pathway_filename <- paste0(comparison_name, "_KEGG_Pathway_", top_pathway, ".pdf")
      pathview(gene.data  = significant_entrez,
               pathway.id = top_pathway,
               species    = "hsa",
               limit      = list(gene = 2, cpd = 1),
               kegg.dir = getwd())
      
      # Rename the output file to include comparison name
      default_file <- paste0("hsa", top_pathway, ".pathview.pdf")
      if (file.exists(default_file)) {
        file.rename(default_file, pathway_filename)
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # 4. REACTOME PATHWAY ENRICHMENT
  # ---------------------------------------------------------------------------
  cat("Running Reactome pathway enrichment...\n")
  
  reactome_enrich <- enrichPathway(gene          = significant_entrez,
                                   organism      = "human",
                                   pvalueCutoff  = pval_cutoff,
                                   pAdjustMethod = "BH",
                                   readable      = TRUE)
  
  if (!is.null(reactome_enrich) && nrow(reactome_enrich) > 0) {
    filename <- paste0(comparison_name, "_Reactome_Pathways.csv")
    write.csv(reactome_enrich@result, filename)
    
    p <- dotplot(reactome_enrich, showCategory = 10, 
                 title = paste("Reactome Pathways -", comparison_name)) +
      large_text_theme
    plot_filename <- paste0(comparison_name, "_Reactome_Dotplot.pdf")
    ggsave(plot_filename, plot = p, width = 12, height = 10, dpi = 300)
  }


  # ---------------------------------------------------------------------------
  # 5. GENE SET ENRICHMENT ANALYSIS (GSEA) - INDIVIDUAL PLOTS VERSION
  # ---------------------------------------------------------------------------
  cat("Running GSEA...\n")
  
  # Prepare ranked gene list
  ranked_genes <- deg_results$logFC
  names(ranked_genes) <- rownames(deg_results)
  ranked_genes <- sort(ranked_genes, decreasing = TRUE, na.last = NA)
  
  # Convert to ENTREZID for GSEA
  ranked_entrez <- bitr(names(ranked_genes),
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)
  
  ranked_entrez_df <- data.frame(ENTREZID = ranked_entrez$ENTREZID,
                                 logFC = ranked_genes[ranked_entrez$SYMBOL])
  ranked_entrez_vector <- ranked_entrez_df$logFC
  names(ranked_entrez_vector) <- ranked_entrez_df$ENTREZID
  
  # GSEA for KEGG
  gsea_kegg <- tryCatch({
    gsea_result <- gseKEGG(geneList     = ranked_entrez_vector,
                           organism     = "hsa",
                           minGSSize    = 10,
                           maxGSSize    = 500,
                           eps = 0,
                           pvalueCutoff = pval_cutoff,
                           verbose      = FALSE)
    cat(paste("GSEA completed successfully with", nrow(gsea_result), "pathways\n"))
    gsea_result
  }, error = function(e) {
    cat("GSEA KEGG failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
    filename <- paste0(comparison_name, "_GSEA_KEGG.csv")
    write.csv(gsea_kegg@result, filename)
    cat(paste("GSEA CSV saved:", filename, "\n"))
  } else {
    cat("Skipping GSEA plotting - no results or empty results\n")
  }


  # ---------------------------------------------------------------------------
  # SAVE ALL RESULTS
  # ---------------------------------------------------------------------------
  save(go_results, disease_enrich, kegg_enrich, reactome_enrich, gsea_kegg,
       file = paste0(comparison_name, "_All_Enrichment_Results.RData"))
  
  cat(paste("Enrichment analysis completed for:", comparison_name, "\n"))
  cat(paste("Results saved in directory:", dir_name, "\n"))
  
  # Return to original working directory
  setwd(original_wd)
  
  return(list(GO = go_results, 
              Disease = disease_enrich, 
              KEGG = kegg_enrich, 
              Reactome = reactome_enrich, 
              GSEA = gsea_kegg))
}
