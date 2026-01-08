#!/bin/bash

#SBATCH -J H5N1_RNA_ALIGN                         			 # Job name
#SBATCH -o ~/WORKING/LOGS/H5N1_RNA_ALIGN.o.log                # Name of the stdout output file
#SBATCH -e ~/WORKING/LOGS/H5N1_RNA_ALIGN.e.log                # Name of the stderr error file
#SBATCH -t 7-00:00:01                       			 # Time of job
#SBATCH -p compute2                       			 # Queue (partition) name
#SBATCH -N 1                             			 # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                              			 # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 80							 # Total number of cores 80 max
# Load one of these

module load anaconda3 && conda activate RNA-seq_NovoGene

# Set variables
INPUT_DIR=~/WORKING/RNA-seq/H5N1/fastq
OUTPUT_DIR=~/WORKING/RNA-seq/H5N1/aligned
REF_DIR=~/WORKING/RNA-seq/H5N1_BACKUP/Ref
THREADS=80

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

# Step 1: Prepare reference genome
echo "Preparing reference genome..."
cd "$REF_DIR"

# Decompress reference files if needed
if [ -f "genome.fa.gz" ]; then
    echo "Decompressing genome.fa.gz..."
    gunzip -c genome.fa.gz > genome.fa
fi

if [ -f "genome.gtf.gz" ]; then
    echo "Decompressing genome.gtf.gz..."
    gunzip -c genome.gtf.gz > genome.gtf
fi

# Generate genome index for HISAT2 if it doesn't exist
HISAT2_INDEX="$REF_DIR/genome_index"
if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
    echo "Building HISAT2 index..."
    hisat2-build -p $THREADS genome.fa "$HISAT2_INDEX"
else
    echo "HISAT2 index already exists, skipping index building..."
fi

# Generate genome.fai index if it doesn't exist
if [ ! -f "genome.fa.fai" ]; then
    echo "Generating genome.fai index..."
    samtools faidx genome.fa
fi

GENOME_FAI="$REF_DIR/genome.fa.fai"

# Step 2: Function to run HISAT2 alignment
run_hisat2() {
    local input_dir="$1"
    local output_dir="$2"
    local hisat_index="$3"
    local genome_fai="$4"
    local base="$5"
    
    echo "Aligning sample: $base"
    
    # Define file paths
    local trimmed_r1="${input_dir}/${base}_1_val_1.fq.gz"
    local trimmed_r2="${input_dir}/${base}_2_val_2.fq.gz"
    local output_bam="${output_dir}/${base}.bam"
    local log_file="${output_dir}/logs/${base}_hisat2.log"
    
    # Check if input files exist
    if [ ! -f "$trimmed_r1" ] || [ ! -f "$trimmed_r2" ]; then
        echo "Error: Input files not found for sample $base"
        echo "R1: $trimmed_r1"
        echo "R2: $trimmed_r2"
        return 1
    fi
    
    # Run HISAT2 alignment and sort BAM
    echo "Running HISAT2 for $base..."
    hisat2 -p 2 \
        -x "$hisat_index" \
        -1 "$trimmed_r1" \
        -2 "$trimmed_r2" \
        2> "$log_file" | \
    samtools sort -@ 2 -o "$output_bam" - && \
    
    # Index the BAM file
    echo "Indexing BAM for $base..."
    samtools index -@ 2 "$output_bam" && \
    
    echo "Completed alignment for: $base"
}

export -f run_hisat2

# Step 3: Run HISAT2 in parallel for all trimmed samples
echo "Starting alignment of trimmed samples..."
find "$INPUT_DIR" -maxdepth 1 -type f -name "*_1_val_1.fq.gz" | xargs -P $THREADS -I {} sh -c '
    base=$(basename "{}" _1_val_1.fq.gz)
    run_hisat2 "$0" "$1" "$2" "$3" "$base"
' "$INPUT_DIR" "$OUTPUT_DIR" "$HISAT2_INDEX" "$GENOME_FAI"

# Wait for all alignment processes to complete
wait

# Step 4: Generate alignment statistics
echo "Generating alignment statistics..."
find "$OUTPUT_DIR" -name "*.bam" | xargs -P $THREADS -I {} sh -c '
    base=$(basename "{}" .bam)
    samtools flagstat "{}" > "$0/logs/${base}_flagstat.txt"
    echo "Generated flagstat for $base"
' "$OUTPUT_DIR"

# Step 5: Run MultiQC on alignment results
echo "Running MultiQC on alignment results..."
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR" --filename "multiqc_alignment_report"

echo "HISAT2 alignment pipeline completed successfully!"
echo "Results available in: $OUTPUT_DIR"
echo "MultiQC report: $OUTPUT_DIR/multiqc_alignment_report.html"

exit
