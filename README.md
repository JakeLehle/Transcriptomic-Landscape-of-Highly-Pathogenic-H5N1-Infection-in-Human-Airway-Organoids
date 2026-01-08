# Transcriptomic Landscape of Highly Pathogenic H5N1 Infection in Human Airway Organoids

This repository contains the bioinformatics pipeline and analysis scripts for bulk RNA-seq analysis of H5N1 influenza virus infection in human airway organoids (HAOs). The workflow covers read alignment through differential expression analysis, pathway enrichment, and publication-ready visualizations.

## Repository Structure

```
Transcriptomic-Landscape-of-Highly-Pathogenic-H5N1-Infection-in-Human-Airway-Organoids/
├── CONDA_ENV/
│   └── environment_R_bulk_RNA-seq_NovoGene.yaml   # Conda environment specification
├── SCRIPTS/
│   ├── HISAT2_ALIGNMENT.sh      # Read alignment pipeline (SLURM)
│   ├── RNA-SEQ_WORK_FLOW.R      # Master analysis workflow
│   ├── VOLCANO_PLOT.R           # Volcano plot generation
│   ├── ENRICHMENT.R             # Pathway enrichment analysis
│   ├── HEATMAP.R                # KEGG pathway heatmaps
│   └── KEGG_TABLE.R             # Pathway comparison tables
├── LICENSE
└── README.md
```

## Requirements

### Setting Up the Conda Environment

```bash
# Create the environment from the YAML file
conda env create -f CONDA_ENV/environment_R_bulk_RNA-seq_NovoGene.yaml

# Activate the environment
conda activate RNA-seq_NovoGene
```

### Additional R Packages

The `RNA-SEQ_WORK_FLOW.R` script automatically installs missing R packages on first run. Key dependencies include:

**CRAN packages:** rlang, devtools, tibble, dplyr, ggplot2, pheatmap, RColorBrewer, ggrepel, gridExtra, tidyr

**Bioconductor packages:** edgeR, DESeq2, Rsubread, org.Hs.eg.db, clusterProfiler, DOSE, ReactomePA, pathview, enrichplot, KEGGREST, biomaRt

## Pipeline Overview

The analysis follows this workflow:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  1. ALIGNMENT (HISAT2_ALIGNMENT.sh)                                         │
│     Raw FASTQ → Trimmed reads → HISAT2 alignment → Sorted BAM files         │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  2. QUANTIFICATION & NORMALIZATION (RNA-SEQ_WORK_FLOW.R)                    │
│     BAM files → featureCounts → DGEList → TMM normalization → CPM/RPKM      │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  3. DIFFERENTIAL EXPRESSION (RNA-SEQ_WORK_FLOW.R)                           │
│     edgeR GLM-QL framework → Three contrasts:                               │
│       • H5N1_24hr vs Mock                                                   │
│       • H5N1_48hr vs Mock                                                   │
│       • H5N1_48hr vs H5N1_24hr                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  4. VISUALIZATION & ENRICHMENT                                              │
│     • Volcano plots (VOLCANO_PLOT.R)                                        │
│     • GO/KEGG/Reactome/Disease enrichment (ENRICHMENT.R)                    │
│     • Pathway heatmaps (HEATMAP.R)                                          │
│     • Pathway comparison tables (KEGG_TABLE.R)                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Script Descriptions

### 1. HISAT2_ALIGNMENT.sh

**Purpose:** SLURM batch script for aligning paired-end RNA-seq reads to a reference genome using HISAT2.

**Key Steps:**
- Builds HISAT2 genome index (if not present)
- Aligns trimmed reads in parallel
- Sorts and indexes BAM files with samtools
- Generates alignment statistics (flagstat)
- Creates MultiQC summary report

**Usage:**
```bash
# Edit paths in script, then submit to SLURM
sbatch SCRIPTS/HISAT2_ALIGNMENT.sh
```

**Inputs:** Trimmed FASTQ files (`*_1_val_1.fq.gz`, `*_2_val_2.fq.gz`), reference genome (FASTA + GTF)

**Outputs:** Sorted BAM files, BAM indices (`.bai`), alignment logs, MultiQC report

---

### 2. RNA-SEQ_WORK_FLOW.R

**Purpose:** Master orchestration script that runs the complete downstream RNA-seq analysis pipeline.

**Key Steps:**
- Installs/loads all required packages
- Runs featureCounts for read quantification
- Creates DGEList with comprehensive gene annotation (ENSEMBL → gene symbol mapping via biomaRt)
- Performs TMM normalization
- Conducts differential expression analysis using edgeR's GLM quasi-likelihood framework
- Exports normalized expression matrices and DEG results
- Calls visualization and enrichment scripts

**Usage:**
```r
# Edit OUTPUT_DIR and SOURCE_DIR paths, then run in R
source("SCRIPTS/RNA-SEQ_WORK_FLOW.R")
```

**Inputs:** BAM files in working directory, GTF annotation file

**Outputs:**
- `logCPM_with_gene_symbols.csv` - Log2 CPM normalized expression matrix
- `logRPKM_with_gene_symbols.csv` - Log2 RPKM normalized expression matrix
- `H5N1_24hr_vs_Mock_DEGs.csv` - Differential expression results
- `H5N1_48hr_vs_Mock_DEGs.csv` - Differential expression results
- `H5N1_48hr_vs_24hr_DEGs.csv` - Differential expression results

---

### 3. VOLCANO_PLOT.R

**Purpose:** Creates publication-ready volcano plots for each differential expression comparison.

**Features:**
- Color-coded significance (red = upregulated, blue = downregulated, grey = not significant)
- Customizable thresholds (default: |logFC| > 1.5, FDR < 0.05)
- Annotated gene counts for up/down regulated genes
- Custom legend with large text formatting

**Usage:** Called automatically by `RNA-SEQ_WORK_FLOW.R`, or:
```r
source("SCRIPTS/VOLCANO_PLOT.R")
volcano_plot <- create_enhanced_volcano(
  deg_results = results_24hr_vs_mock,
  comparison_name = "H5N1_24hr_vs_Mock"
)
```

**Outputs:** PNG volcano plots (e.g., `Volcano_Plot_H5N1_24hr_vs_Mock.png`)

---

### 4. ENRICHMENT.R

**Purpose:** Comprehensive functional enrichment analysis of differentially expressed genes.

**Analyses Performed:**
- Gene Ontology (GO) enrichment (Biological Process, Cellular Component, Molecular Function)
- Disease Ontology (DO) enrichment
- KEGG pathway enrichment
- Reactome pathway enrichment
- Gene Set Enrichment Analysis (GSEA) for KEGG

**Usage:** Called automatically by `RNA-SEQ_WORK_FLOW.R`, or:
```r
source("SCRIPTS/ENRICHMENT.R")
enrichment_results <- run_comprehensive_enrichment(
  deg_results = results_24hr_vs_mock,
  comparison_name = "H5N1_24hr_vs_Mock",
  y_all_object = y_all
)
```

**Outputs:** For each comparison, creates a directory containing:
- CSV files for each enrichment type
- Dotplot PDFs for significant pathways
- Enrichment map visualizations
- Pathview pathway diagrams
- Combined RData object

---

### 5. HEATMAP.R

**Purpose:** Creates expression heatmaps for genes within specific KEGG pathways.

**Features:**
- Fetches pathway genes via KEGGREST API
- Row-scaled (z-score) expression values
- Hierarchical clustering of genes and samples
- Sample group color annotations
- Displays normalized expression values in cells

**Usage:** Called automatically by `RNA-SEQ_WORK_FLOW.R`, or:
```r
source("SCRIPTS/HEATMAP.R")
pathway_genes <- create_kegg_pathway_heatmap(
  pathway_id = "hsa04060",  # Cytokine-cytokine receptor interaction
  y_all_object = y_all,
  output_filename = "Heatmap_Cytokine_Pathway.pdf"
)
```

**Outputs:** PDF heatmaps for specified KEGG pathways

---

### 6. KEGG_TABLE.R

**Purpose:** Generates colored comparison tables showing logFC values across timepoints for pathway-specific genes.

**Features:**
- Retrieves pathway genes via KEGGREST API with retry logic
- Color-gradient visualization (blue = downregulated, white = unchanged, red = upregulated)
- Filters genes by logFC threshold
- Batch processing for multiple pathways
- Summary statistics report

**Usage:** Called automatically by `RNA-SEQ_WORK_FLOW.R`, or:
```r
source("SCRIPTS/KEGG_TABLE.R")

# Single pathway
results <- create_pathway_comparison_table_keggrest_only(
  deg_24hr = results_24hr_vs_mock,
  deg_48hr = results_48hr_vs_mock,
  kegg_pathway_id = "hsa04060",
  logFC_threshold = 1.5,
  output_dir = getwd()
)

# Batch processing
pathway_list <- c("hsa04060", "hsa04620", "hsa04622")
batch_results <- batch_pathway_analysis_keggrest_only(
  deg_24hr = results_24hr_vs_mock,
  deg_48hr = results_48hr_vs_mock,
  pathway_list = pathway_list,
  logFC_threshold = 1.5,
  output_dir = getwd()
)
```

**Outputs:**
- CSV files with pathway gene expression data
- Color-coded PDF tables
- Summary report (`Pathway_Analysis_Summary_KEGGREST.txt`)

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/[username]/Transcriptomic-Landscape-of-Highly-Pathogenic-H5N1-Infection-in-Human-Airway-Organoids.git
cd Transcriptomic-Landscape-of-Highly-Pathogenic-H5N1-Infection-in-Human-Airway-Organoids

# 2. Set up conda environment
conda env create -f CONDA_ENV/environment_R_bulk_RNA-seq_NovoGene.yaml
conda activate RNA-seq_NovoGene

# 3. Prepare your reference genome and place trimmed FASTQs in input directory

# 4. Edit paths in HISAT2_ALIGNMENT.sh and submit alignment job
sbatch SCRIPTS/HISAT2_ALIGNMENT.sh

# 5. After alignment completes, edit paths in RNA-SEQ_WORK_FLOW.R
#    Set OUTPUT_DIR to your BAM file directory
#    Set SOURCE_DIR to the SCRIPTS directory

# 6. Run the analysis in R
Rscript SCRIPTS/RNA-SEQ_WORK_FLOW.R
# Or interactively in RStudio
```

## Expected Output Files

| Category | Files | Description |
|----------|-------|-------------|
| Normalized Expression | `logCPM_with_gene_symbols.csv`, `logRPKM_with_gene_symbols.csv` | Gene expression matrices |
| Differential Expression | `H5N1_*_DEGs.csv` | DEG results with logFC, FDR, significance |
| Volcano Plots | `Volcano_Plot_*.png` | Publication-ready volcano plots |
| Enrichment Results | `Enrichment_Results_*/` | Directory per comparison with all enrichment outputs |
| Pathway Heatmaps | `Heatmap_*_Pathway.pdf` | Expression heatmaps for specific pathways |
| Pathway Tables | `Pathway_Table_Colored_*.pdf`, `Pathway_Data_*.csv` | Comparative pathway gene tables |

## Citation

*[Citation information will be added upon publication]*

## License

See [LICENSE](LICENSE) for details.

## Contact

For questions regarding the analysis pipeline or data, please open an issue in this repository.
