# Loading screen section

# Define all required packages
cran_packages <- c(
  "rlang", "devtools", "stats", "statmod", "tibble", "dplyr", 
  "ggplot2", "pheatmap", "RColorBrewer", "ggtangle", "ggrepel", 
  "grid", "gridExtra", "tidyr", "FSA"
)

bioc_packages <- c(
  "BiocManager", "edgeR", "DESeq2", "Rsubread", "org.Hs.eg.db",
  "clusterProfiler", "DOSE", "ReactomePA", "pathview", "enrichplot",
  "KEGGREST", "biomaRt"
)

# Function to install and load CRAN packages
install_if_missing_cran <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

# Function to install and load Bioconductor packages
install_if_missing_bioc <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) BiocManager::install(new_packages)
}

# Install missing packages
install_if_missing_cran(cran_packages)
install_if_missing_bioc(bioc_packages)

# Load all packages at once
invisible(lapply(c(cran_packages, bioc_packages), library, character.only = TRUE))

# Print loaded packages for verification
cat("Successfully loaded", length(c(cran_packages, bioc_packages)), "packages\n")

# --------------------------------------------------------------
# Main bluk RNA-seq work flow
# --------------------------------------------------------------

OUTPUT_DIR <- "/path/to/your/output/dir/must/contain/raw/fastqs"
SOURCE_DIR <- "path/to/dir/with/scripts/"

# Step 1: prepare your own reference genome index.
setwd(OUTPUT_DIR)

bam_files <- list.files(pattern = "\\.bam$")

print(c('Here is the top 10 items', head(bam_files, n=10)))

source(paste0(SOURCE_DIR, "FEATURE_COUNTS.R"))

fc_obj <- compact_featureCounts(bam_files)

# You can check that your fc object was set correctly using head() and specificng to return the information in the counts 
head(fc_obj$counts)
ncol(fc_obj$counts)
fc_obj$targets
fc_obj$annotation$Length
fc_obj$annotation$GeneID

# edgeR analysis
group_all <- ifelse(grepl("24h", bam_files), "H5N1_24hr",
             ifelse(grepl("48h", bam_files), "H5N1_48hr",
             ifelse(grepl("Mock", bam_files), "Mock", NA)))
if(any(is.na(group_all))) stop("Some files couldn't be assigned to a group!")
sample_metadata <- data.frame(sample_name = bam_files, group = group_all) 

# Check the metadata
print(sample_metadata)

# Make the DEGList object
y_all <- DGEList(counts=fc_obj$counts, group=sample_metadata$group)

# Get gene lengths for human genes
# First map your gene symbols to ENTREZ IDs
# Add gene lengths to y_all$genes before filtering
# First extract the geneID and lnegth from the featureCount output
gene_annotation <- data.frame(
  ENSEMBL = fc_obj$annotation$GeneID,
  Length = fc_obj$annotation$Length
)
# Now get the common gene symbols to swap that out with the ENSEMBL ID they are currently using
GeneSymbol <- mapIds(org.Hs.eg.db, 
                     keys = rownames(y_all), 
                     keytype = "ENSEMBL", 
                     column = "SYMBOL")


gene_annotation$SYMBOL <- GeneSymbol

# Get chromosome locations using biomaRt
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get chromosome, start, and end positions
gene_locations <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = gene_annotation$ENSEMBL,
  mart = mart
)

# Sort data to make sure they match up 
gene_locations <- gene_locations %>%
  mutate(ensembl_gene_id = as.character(ensembl_gene_id)) %>%
  arrange(match(ensembl_gene_id, gene_annotation$ENSEMBL))

# Merge with gene annotation
gene_annotation <- gene_annotation %>%
  left_join(gene_locations, by = c("ENSEMBL" = "ensembl_gene_id")) %>%
  rename(
    CHR = chromosome_name,
    START = start_position,
    END = end_position
  )

gene_annotation$SYMBOL_FINAL <- ifelse(is.na(gene_annotation$SYMBOL),
                                       gene_annotation$ENSEMBL,
                                       gene_annotation$SYMBOL)

symbol_freq <- table(gene_annotation$SYMBOL_FINAL)
duplicate_symbols <- names(symbol_freq[symbol_freq > 1])

gene_annotation$SYMBOL_UNIQUE <- gene_annotation$SYMBOL_FINAL

if(length(duplicate_symbols) > 0) {
  warning(paste("Found", length(duplicate_symbols), "duplicate identifiers. Appending Ensembl ID to make them unique."))
  
  # For all these duplicate rows add in a custom symbol by combining the gene symbo and the ensembl ID
  is_duplicate <- gene_annotation$SYMBOL_FINAL %in% duplicate_symbols
  
  gene_annotation$SYMBOL_UNIQUE[is_duplicate] <- paste(gene_annotation$SYMBOL_FINAL[is_duplicate],
                                         gene_annotation$ENSEMBL[is_duplicate],
                                         sep = "_")
}

symbol_freq <- table(gene_annotation$SYMBOL_UNIQUE)
duplicate_symbols <- names(symbol_freq[symbol_freq > 1])

# Assign the comprehensive data frame to y_all$genes
if(all(rownames(y_all$counts) == gene_annotation$ENSEMBL)) {
  print("Everything matched up! Adding Gene symbols and gene lengths.")
  y_all$genes <- gene_annotation
} else {
  warning("Oops looks like something wasn't in order I'll try to fix that now. Double check this data!")
  gene_annotation <- gene_annotation[match(rownames(y_all$counts), gene_annotation$ENSEMBL), ]
  y_all$genes <- gene_annotation
}

# Set rownames of the DGEList to the unique identifiers
rownames(y_all$counts) <- y_all$genes$SYMBOL_UNIQUE
rownames(y_all$genes) <- y_all$genes$SYMBOL_UNIQUE

# Inspect the result to see how NAs and duplicates were handled
head(y_all$genes, n = 15)
cat("\nSample of genes where original SYMBOL was NA:\n")
print(head(y_all$genes[is.na(y_all$genes$SYMBOL), ]))
cat("\nSample of genes with duplicate SYMBOL_FINAL:\n")
dupe_examples <- y_all$genes[y_all$genes$SYMBOL_FINAL %in% head(duplicate_symbols), ]
print(dupe_examples[order(dupe_examples$SYMBOL_FINAL), ])

#Filter out genes that have no expression
keepy_all <- filterByExpr(y_all)
y_all <- y_all[keepy_all, , keep.lib.sizes=FALSE]

#Normalize the library size using TMM normalization this is log normalization of the mean (M)
y_all <- calcNormFactors(y_all, method = "TMM")

# To make a standard cpm adjusted count matrix that does not take into account the gene length for the normalization
CPM <- cpm(y_all, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

# To make a rpkm count matrix that does take into account the gene length you need to find the gene length 
RPKM <- rpkm(y_all, gene.length = y_all$genes$Length, normalized.lib.sizes = TRUE, log =TRUE, prior.count = 1)
head(RPKM)

# Okay now I'm going to write out both files
logCPM_with_symbols <- as.data.frame(CPM) %>%
  rownames_to_column(var = "GeneSymbol")

write.csv(logCPM_with_symbols,
          file = "logCPM_with_gene_symbols.csv",
          row.names = FALSE,
          quote = FALSE)

# Same thing for RPKM
logRPKM_with_symbols <- as.data.frame(RPKM) %>%
  rownames_to_column(var = "GeneSymbol")

write.csv(logRPKM_with_symbols,
          file = "logRPKM_with_gene_symbols.csv",
          row.names = FALSE,
          quote = FALSE)

# Okay now lets calculate DEG
design <- model.matrix(~ 0 + group, data = y_all$samples)

# Give the columns meaningful names based on the group levels
colnames(design) <- levels(y_all$samples$group)
print("Design matrix:")
print(design)

# 2. Estimate dispersions
#    This critical step estimates the biological variability for each gene
y_all <- estimateDisp(y_all, design, robust = TRUE)

#    Visualize the dispersion estimates (highly recommended)
plotBCV(y_all, main = "Biological Coefficient of Variation")

# 3. Fit the negative binomial generalized linear model (GLM)
fit <- glmQLFit(y_all, design, robust = TRUE)

#    Visualize the quasi-likelihood dispersions (optional but useful)
plotQLDisp(fit, main = "QL Dispersions")

# 4. Define contrasts and perform statistical testing
#    Contrast: H5N1_24hr vs Mock
con_24hr_vs_mock <- makeContrasts(H5N1_24hr_vs_Mock = H5N1_24hr - Mock, levels = design)
qlf_24hr_vs_mock <- glmQLFTest(fit, contrast = con_24hr_vs_mock)

#    Contrast: H5N1_48hr vs Mock 
con_48hr_vs_mock <- makeContrasts(H5N1_48hr_vs_Mock = H5N1_48hr - Mock, levels = design)
qlf_48hr_vs_mock <- glmQLFTest(fit, contrast = con_48hr_vs_mock)

#    Contrast: H5N1_48hr vs H5N1_24hr (to see time-dependent changes)
con_48hr_vs_24hr <- makeContrasts(H5N1_48hr_vs_24hr = H5N1_48hr - H5N1_24hr, levels = design)
qlf_48hr_vs_24hr <- glmQLFTest(fit, contrast = con_48hr_vs_24hr)

# 5. Extract and examine results
#    For H5N1_24hr vs Mock
print("Top 10 DEGs for H5N1_24hr vs Mock:")
print(topTags(qlf_24hr_vs_mock, n = 10))

# Get full results with multiple testing correction
results_24hr_vs_mock <- topTags(qlf_24hr_vs_mock, n = nrow(y_all$counts), 
                               adjust.method = "BH", sort.by = "P")$table

# Add a column indicating significance (e.g., FDR < 0.05 and |logFC| > 1)
results_24hr_vs_mock$significant <- ifelse(results_24hr_vs_mock$FDR < 0.05 & 
                                          abs(results_24hr_vs_mock$logFC) > 1.5, 
                                         "yes", "no")

# Check the results
print(paste("Number of significant DEGs (FDR < 0.05, |logFC| > 1):", 
            sum(results_24hr_vs_mock$significant == "yes")))
print("Summary of results:")
print(summary(results_24hr_vs_mock))

# Repeat for other comparisons
results_48hr_vs_mock <- topTags(qlf_48hr_vs_mock, n = nrow(y_all$counts), 
                               adjust.method = "BH", sort.by = "P")$table
results_48hr_vs_mock$significant <- ifelse(results_48hr_vs_mock$FDR < 0.05 & 
                                          abs(results_48hr_vs_mock$logFC) > 1.5, 
                                         "yes", "no")

results_48hr_vs_24hr <- topTags(qlf_48hr_vs_24hr, n = nrow(y_all$counts), 
                               adjust.method = "BH", sort.by = "P")$table
results_48hr_vs_24hr$significant <- ifelse(results_48hr_vs_24hr$FDR < 0.05 & 
                                          abs(results_48hr_vs_24hr$logFC) > 1.5, 
                                         "yes", "no")

# Save all results
write.csv(results_24hr_vs_mock, file = "H5N1_24hr_vs_Mock_DEGs.csv")
write.csv(results_48hr_vs_mock, file = "H5N1_48hr_vs_Mock_DEGs.csv")
write.csv(results_48hr_vs_24hr, file = "H5N1_48hr_vs_24hr_DEGs.csv")

# -----------------------------------------------------------------------------
# VOLCANO PLOTS (FIGURE 2)
# -----------------------------------------------------------------------------

source(paste0(SOURCE_DIR, "VOLCANO_PLOT.R"))


# Create the volcano plot for 24hr_vs_Mock
volcano_24hr_vs_mock <- create_enhanced_volcano(
  deg_results = results_24hr_vs_mock,
  comparison_name = "H5N1_24hr_vs_Mock"
)

# Create the volcano plot for 24hr_vs_Mock
volcano_48hr_vs_mock <- create_enhanced_volcano(
  deg_results = results_48hr_vs_mock,
  comparison_name = "H5N1_48hr_vs_Mock"
)

volcano_48hr_vs_24hr <- create_enhanced_volcano(
  deg_results = results_48hr_vs_24hr,
  comparison_name = "H5N1_48hr_vs_24hr"
)

# -----------------------------------------------------------------------------
# MASTER ENRICHMENT ANALYSIS FUNCTION (FIGURES 3-4)
# -----------------------------------------------------------------------------

source(paste0(SOURCE_DIR, "ENRICHMENT.R"))

# -----------------------------------------------------------------------------
# RUN ENRICHMENT FOR ALL COMPARISONS
# -----------------------------------------------------------------------------

# Run enrichment for each comparison
enrichment_24hr_vs_mock <- run_comprehensive_enrichment(
  deg_results = results_24hr_vs_mock,
  comparison_name = "H5N1_24hr_vs_Mock",
  y_all_object = y_all
)

enrichment_48hr_vs_mock <- run_comprehensive_enrichment(
  deg_results = results_48hr_vs_mock,
  comparison_name = "H5N1_48hr_vs_Mock",
  y_all_object = y_all
)

enrichment_48hr_vs_24hr <- run_comprehensive_enrichment(
  deg_results = results_48hr_vs_24hr,
  comparison_name = "H5N1_48hr_vs_24hr",
  y_all_object = y_all
)

# -----------------------------------------------------------------------------
# KEGG HEATMAP PLOTS (SUPPLEMENTAL FIGURES)
# -----------------------------------------------------------------------------

source(paste0(SOURCE_DIR, "HEATMAP.R"))

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04060",
  y_all_object = y_all,
  output_filename = "Heatmap_Cytokine_Pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04668",
  y_all_object = y_all,
  output_filename = "Heatmap_TNF_Pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04622",
  y_all_object = y_all,
  output_filename = "Heatmap_RIG-I-like_receptor_Pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04620",
  y_all_object = y_all,
  output_filename = "Heatmap_Toll-like_receptor_Pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04512",
  y_all_object = y_all,
  output_filename = "Heatmap_ECM-receptor_interaction_Pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04064",
  y_all_object = y_all,
  output_filename = "Heatmap_NF-kappa_B_signaling_pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04217",
  y_all_object = y_all,
  output_filename = "Heatmap_Necroptosis_pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04630",
  y_all_object = y_all,
  output_filename = "Heatmap_JAK-STAT_signaling_pathway_pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04010",
  y_all_object = y_all,
  output_filename = "Heatmap_MAPK_signaling_pathway_pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04151",
  y_all_object = y_all,
  output_filename = "Heatmap_PI3K-Akt_signaling_pathway_pathway.pdf"
)

########### NON-Significant Sections
cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04310",
  y_all_object = y_all,
  output_filename = "Heatmap_Wnt_signaling_pathway_pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04150",
  y_all_object = y_all,
  output_filename = "Heatmap_m-TOR_signaling_pathway_pathway.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04152",
  y_all_object = y_all,
  output_filename = "Heatmap_AMPK_signaling_pathway_pathway.pdf"
)

# Cancer Pathway Heatmaps

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa05200",
  y_all_object = y_all,
  output_filename = "Heatmap_Pathways_in_cancer.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa05202",
  y_all_object = y_all,
  output_filename = "Heatmap_Transcriptional_misregulation_in_cancer.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa05204",
  y_all_object = y_all,
  output_filename = "Heatmap_Chemical_carcinogenesis-DNA_adducts.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa05208",
  y_all_object = y_all,
  output_filename = "Heatmap_Chemical_carcinogenesis-reactive_oxygen_species.pdf"
)

cytokine_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa05222",
  y_all_object = y_all,
  output_filename = "Heatmap_Small_cell_lung_cancer.pdf"
)


# -----------------------------------------------------------------------------
# KEGG TABLE SECTION (FIGURES 5-9)
# -----------------------------------------------------------------------------

# Example usage in your main pipeline:
source(paste0(SOURCE_DIR, "KEGG_TABLE.R"))

 # Define pathways of interest

pathway_list <- c("hsa04060", "hsa04620", "hsa04622", "hsa04668", "hsa04630", 
                  "hsa04064", "hsa04217", "hsa04310", "hsa04152", "hsa04010", "hsa04512", "hsa04150", "hsa04151")
# Cancer pathways
pathway_list <- c("hsa05200", "hsa05202", "hsa05204", "hsa05208", "hsa05222")

# I had to increase this one to logFC 4 to get it under 50 genes
pathway_list <- c("hsa05200")
# Run batch analysis
pathway_results <- batch_pathway_analysis_keggrest_only(
  deg_24hr = results_24hr_vs_mock,
  deg_48hr = results_48hr_vs_mock,
  pathway_list = pathway_list,
  logFC_threshold = 3,
  output_dir = OUTPUT_DIR
)
########################################################################
#
# KEGG BARPLOT SECTION
#

# Example usage in your main pipeline:
source(paste0(SOURCE_DIR, "KEGG_BARPLOT.R"))

# Single pathway - faceted style
result <- create_pathway_barplot(
  deg_24hr = results_24hr_vs_mock,
  deg_48hr = results_48hr_vs_mock,
  kegg_pathway_id = "hsa04060",
  logFC_threshold = 3,
  output_dir = OUTPUT_DIR,
  plot_style = "faceted"
)

# Single pathway - side-by-side style
result <- create_pathway_barplot(
  deg_24hr = results_24hr_vs_mock,
  deg_48hr = results_48hr_vs_mock,
  kegg_pathway_id = "hsa04060",
  logFC_threshold = 3,
  output_dir = OUTPUT_DIR,
  plot_style = "sidebyside"
)

# Batch processing for multiple pathways
pathway_list <- c("hsa04060", "hsa04620", "hsa04622", "hsa04668", "hsa04630", 
                  "hsa04064", "hsa04217", "hsa04310", "hsa04152", "hsa04010", 
                  "hsa04512", "hsa04150", "hsa04151")
results <- batch_pathway_barplots(
  deg_24hr = results_24hr_vs_mock,
  deg_48hr = results_48hr_vs_mock,
  pathway_list = pathway_list,
  logFC_threshold = 3,
  output_dir = OUTPUT_DIR,
  plot_style = "sidebyside"
)

# -----------------------------------------------------------------------------
# COSMIC ONCOGENE SECTION
# -----------------------------------------------------------------------------

# Source the script
source(paste0(SOURCE_DIR, "COSMIC_ONCOGENE.R"))

# Run the analysis
cosmic_results <- cosmic_oncogene_analysis(
  dgelist = y_all,                        
  results_24hr = results_24hr_vs_mock,    
  results_48hr = results_48hr_vs_mock,    
  rpkm_matrix = RPKM,                     
  cosmic_dir = "COSMIC_data",            
  output_file = "COSMIC_oncogene_analysis.tsv",
  fdr_threshold = 0.05
)

# -----------------------------------------------------------------------------
# EXIT
# -----------------------------------------------------------------------------
print('Analysis complete!')
q()
