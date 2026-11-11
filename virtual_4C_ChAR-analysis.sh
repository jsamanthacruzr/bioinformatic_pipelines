#!/usr/bin/env bash
# ==============================================================================
# Virtual 4C Analysis & Genomic Interactions Pipeline (ChAR-seq / RIC-seq / Hi-C)
# Description: Processes SRA datasets, runs BLAST searches for lncRNAs, 
#              filters unique reads, formats BEDPE interaction files, 
#              and extracts virtual 4C profiles using Juicer Tools.
# Author: Bioinformatics Workflow
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. Download SRA Data & Database Setup
# ------------------------------------------------------------------------------
echo "=== Step 1: Downloading SRA datasets and building BLAST DB ==="
SRA_ID="SRR33765893"

# Download paired-end FASTQ files
fasterq-dump --split-files "${SRA_ID}"

# Make BLAST database for lncRNAs
makeblastdb -in lncRNAs.fasta -dbtype nucl -out lncRNAs_db

# ------------------------------------------------------------------------------
# 2. Convert & Split FASTA Files for Parallel Processing
# ------------------------------------------------------------------------------
echo "=== Step 2: Converting FASTQ to FASTA and splitting into 20 parts ==="
seqkit fq2fa "${SRA_ID}_1.fastq" -o "${SRA_ID}_R1.fasta"
seqkit fq2fa "${SRA_ID}_2.fastq" -o "${SRA_ID}_R2.fasta"

seqkit split -p 20 "${SRA_ID}_R1.fasta"
seqkit split -p 20 "${SRA_ID}_R2.fasta"

# ------------------------------------------------------------------------------
# 3. Parallel BLAST Search & Hit Concatenation
# ------------------------------------------------------------------------------
echo "=== Step 3: Running BLASTn on split files ==="
# Loop over split FASTA files for R1
for file in "${SRA_ID}_R1.fasta.split"/*.fasta; do
    blastn -db lncRNAs_db -query "${file}" -outfmt 6 -out "${file}.hits.tsv"
done

# Loop over split FASTA files for R2
for file in "${SRA_ID}_R2.fasta.split"/*.fasta; do
    blastn -db lncRNAs_db -query "${file}" -outfmt 6 -out "${file}.hits.tsv"
done

# Concatenate BLAST hits for R1 and R2 separately
cat "${SRA_ID}_R1.fasta.split"/*.hits.tsv > "${SRA_ID}_blast-hits_R1.tsv"
cat "${SRA_ID}_R2.fasta.split"/*.hits.tsv > "${SRA_ID}_blast-hits_R2.tsv"

# ------------------------------------------------------------------------------
# 4. Filter Alignment Length (> 20 nt) & Extract Unique Read IDs
# ------------------------------------------------------------------------------
echo "=== Step 4: Filtering hits by alignment length (> 20 nt) ==="
awk -v minlen=20 '$4 > minlen {print}' "${SRA_ID}_blast-hits_R1.tsv" > "${SRA_ID}_blast-hits_filtered_R1.tsv"
awk -v minlen=20 '$4 > minlen {print}' "${SRA_ID}_blast-hits_R2.tsv" > "${SRA_ID}_blast-hits_filtered_R2.tsv"

# Log line counts for QA
wc -l "${SRA_ID}_blast-hits_filtered_R1.tsv"
wc -l "${SRA_ID}_blast-hits_filtered_R2.tsv"

# Extract unique read IDs
cut -f1 "${SRA_ID}_blast-hits_filtered_R1.tsv" | sort -u > "${SRA_ID}_IDs_R1.txt"
cut -f1 "${SRA_ID}_blast-hits_filtered_R2.tsv" | sort -u > "${SRA_ID}_IDs_R2.txt"

# Combine IDs from both reads
cat "${SRA_ID}_IDs_R1.txt" "${SRA_ID}_IDs_R2.txt" | sort -u > ids_combined.txt

# Extract matching reads from original FASTA files
seqkit grep -f ids_combined.txt "${SRA_ID}_R1.fasta" -o "${SRA_ID}_filtered_R1.fasta"
seqkit grep -f ids_combined.txt "${SRA_ID}_R2.fasta" -o "${SRA_ID}_filtered_R2.fasta"

# ------------------------------------------------------------------------------
# 5. Experimental Notes: ChAR-seq Bridge Sequence
# ------------------------------------------------------------------------------
# 25 µM annealed ChAR-seq bridge design:
# Top strand:    5'- /5rApp/ AANNNAAACCGGCGTCCAAGGATCTTTAATTAAGTCGCAG /3SpC3/ -3'
# Bottom strand: 5'- /5Phos/ GATCTGCGACTTAATTAAAGATCCTTGGACGCCGG /iBiodT/ T -3'

# ------------------------------------------------------------------------------
# 6. Join RNA & DNA Coordinates and Generate Weighted BEDPE
# ------------------------------------------------------------------------------
echo "=== Step 6: Joining RNA-DNA datasets and generating weighted BEDPE ==="

# Ensure input BED files are sorted by Read ID (Column 4) before joining
sort -k4,4 RNA.clean.bed > RNA.clean.sorted.bed
sort -k4,4 DNA.clean.bed > DNA.clean.sorted.bed

# Join RNA and DNA interactions by Read ID
join -t $'\t' -1 4 -2 4 RNA.clean.sorted.bed DNA.clean.sorted.bed > joined_lncs.highq.tmp

# Convert joined file to official 10-column BEDPE format with normalized score
awk 'BEGIN{OFS="\t"}{
    id = $1;
    chrR = $2; startR = $3; endR = $4; nhR = $5; strandR = $6;
    chrD = $7; startD = $8; endD = $9; nhD = $10; strandD = $11;

    # Protect against NH = 0
    if(nhR < 1) nhR = 1;
    if(nhD < 1) nhD = 1;

    # Calculate average weight based on multi-mapper hits (NH)
    weight = (1/nhR + 1/nhD) / 2;

    # Scale score for IGV visualization
    score = int(weight * 1000);

    # Output standard 10-column BEDPE format
    print chrR, startR, endR, chrD, startD, endD, id, score, strandR, strandD;
}' joined_lncs.highq.tmp > lncs_interactions.weighted.bedpe

# ------------------------------------------------------------------------------
# 7. Virtual 4C Extraction & Conversion (Juicer Tools)
# ------------------------------------------------------------------------------
echo "=== Step 7: Generating Virtual 4C profiles from Hi-C & RIC-seq ==="

# Extract observed contacts around target locus (K562 Hi-C at 10 kb resolution)
java -Xmx4g -jar juicer_tools_1.22.01.jar dump observed NONE \
    K562_combined.hic chr17:41380229:41471132 chr17 BP 10000 \
    K562_virtual4C_10k.txt

# Convert Hi-C matrix output to BEDPE format (10 kb resolution)
awk -v binsize=10000 '$3 >= 1 {
    print "chr17\t"$1"\t"($1+binsize)"\tchr17\t"$2"\t"($2+binsize)"\t"$3
}' K562_virtual4C_10k.txt > K562_virtual4C_10k.bedpe

# Extract virtual 4C at promoter target locus (TILR1 - 5 kb resolution)
java -Xmx4g -jar juicer_tools_1.22.01.jar dump observed NONE \
    K562_combined.hic chr17:41379500:41384500 chr17 BP 5000 \
    K562_virtual4C_TILR1_5k.txt

# Extract virtual 4C contacts from RIC-seq (HeLa dataset at 5 kb resolution)
java -Xmx4g -jar juicer_tools_1.22.01.jar dump observed NONE \
    GSM3629921_HeLa_Total_rep1.hic chr17:41376457:41469114 chr17 BP 5000 \
    HeLa_RICseq_5k.txt

# Convert RIC-seq dump to BEDPE format (5 kb resolution)
awk -v binsize=5000 '$3 >= 1 {
    print "chr17\t"$1"\t"($1+binsize)"\tchr17\t"$2"\t"($2+binsize)"\t"$3
}' HeLa_RICseq_5k.txt > HeLa_RICseq_5k.bedpe

echo "=== Pipeline finished successfully ==="
