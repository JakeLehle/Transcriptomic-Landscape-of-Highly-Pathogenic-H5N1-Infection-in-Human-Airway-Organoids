#!/bin/bash

#SBATCH -J H5N1_RNA_TRIM                         			 # Job name
#SBATCH -o /work/sdz852/WORKING/LOGS/H5N1_RNA_TRIM.o.log                # Name of the stdout output file
#SBATCH -e /work/sdz852/WORKING/LOGS/H5N1_RNA_TRIM.e.log                # Name of the stderr error file
#SBATCH --mail-user=jake.lehle@utsa.edu   			 # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00:01                       			 # Time of job
#SBATCH -p compute2                       			 # Queue (partition) name
#SBATCH -N 1                             			 # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                              			 # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 80							 # Total number of cores 80 max
# Load one of these

module load anaconda3 && conda activate RNA-seq_NovoGene

# Set variables
INPUT_DIR=/work/sdz852/WORKING/RNA-seq/H5N1_BACKUP/01.RawData
OUTPUT_DIR=/work/sdz852/WORKING/RNA-seq/H5N1/JAKE/fastq
THREADS=4

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run trim_galore in parallel for gzipped files
find "$INPUT_DIR" -maxdepth 1 -type f -name "*_1.fq.gz" | xargs -P $THREADS -I {} sh -c '
    base=$(basename "{}" _1.fq.gz)
    dir=$(dirname "{}")
    echo "Processing sample: $base"
    trim_galore --paired --cores 2 --output_dir "$0" \
        "$dir/${base}_1.fq.gz" \
        "$dir/${base}_2.fq.gz"
' "$OUTPUT_DIR"

# Wait for all trim_galore processes to complete
wait

# Run MultiQC on the output directory
echo "Running MultiQC on trimmed files..."
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR" --filename "multiqc_report"


exit
