#!/usr/bin/env bash
# These are the commands used for processing the whole metagenomic sequencing data from the halite longitudinal study.
# The pipeline assumes that metaWRAP is installed with default settings, and is available in the PATH.
# It is assumed that the fastq raw read files are in teh "RAW_READS" directory.

# 1. Preprocessing the reads, trimming, human read contaminaiton removal:
tar -xvf RAW_READS/*tar.gz
mkdir READ_QC
for i in RAW_READS/*_1.fastq; do
	metawrap read_qc -1 ${i%_*}_1.fastq -2 ${i%_*}_2.fastq -t 24 -o READ_QC/${i%_*}
done

mkdir CLEAN_READS
for i in READ_QC/*; do
	b=${i#*/}
	mv ${i}/final_pure_reads_1.fastq CLEAN_READS/${b}_1.fastq
	mv ${i}/final_pure_reads_2.fastq CLEAN_READS/${b}_2.fastq
done

# At this point make sure that the files in "CLEAN_READS" are named as something meaningfull. For this example, I will rename the samples with the naming convention SITE-YEAR-MONTH-REPLICATE. For example: SG1-N-2017-02-1_1.fastq or SG1-2016-02-4_2.fastq


# 2. Metagenomic co-assembly of all the samples with metaSPAdes:
cat CLEAN_READS/*_1.fastq >> CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/*_2.fastq >> CLEAN_READS/ALL_READS_2.fastq
metawrap assembly --metaspades -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -m 1000 -t 48 -o ASSEMBLY
# With 48 threads, this took ~72 hours. Signifficant memory allocaiton (600BG+) was also necessary. If this is not possible, use MegaHIT instead (--megahit option), which requires signifficantly less resources and produces a comparable assembly.


# 3. Running taxonomic composition prediciton with KRAKEN:
metawrap kraken -o KRAKEN -t 48 -s 1000000 CLEAN_READS/*fastq ASSEMBLY/final_assembly.fasta
# Subsetting the reads with the -s option is recommended in the interest of time.


# 4. Binning of scaffolds to produce metagenome assembled genomes:
# initial binning with metabat2, concoct, and maxbin2:
metawrap binning -o INITIAL_BINNING -t 48 -a ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct CLEAN_READS/*fastq
# inerative refinement with metaWRAP bin_refinement module (-c defines minimum desided completion, and -x maximum contamination):
metawrap bin_refinement -o BIN_REFINEMENT -t 48 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 70 -x 5
# MAG reassembly with the Reassemble_bins module was skipped for this analysis due to heterogeneity concerns, but is optional for thi analysis:
#metawrap reassemble_bins -o BIN_REASSEMBLY -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -t 48 -m 800 -c 70 -x 10 -b BIN_REFINEMENT/metaWRAP_bins


# 5. MAG and contig quantification and annotation.
# Note: that from here onwards we only use Site 1 (SG1) reads.
# MAG and contig abundance estimation. This produces both the MAG abundance matrix (QUANT_BINS/abundance_table.tab) and the contig abundance matrix (QUANT_BINS/contig_abundance.tab).
metawrap quant_bins -b BIN_REFINEMENT/metaWRAP_bins -o QUANT_BINS -a ASSEMBLY/final_assembly.fa CLEAN_READS/SG1-201*fastq
# Community visualization with Blobology:
metawrap blobology -a ASSEMBLY/final_assembly.fasta -t 48 -o BLOBOLOGY --bins BIN_REFINEMENT/metaWRAP_bins CLEAN_READS/SG1-201*fastq
# Functional annotations of the entire metagenomic assembly and the MAGs were obtained through the Integrated Microbial Genomes and Microbiomes (IMG) annotation service (through JGI). See manuscript for data availability and hosting details.
# An in-house script (process_img_annotation/rename_contigs.py) was used to rename the contigs in the IMG product_names file to the original contig names.
python rename_contigs.py *names_map *product_names > img_annotation.products
# An in-house script (process_img_annotation/make_master_table.py) was used to link the KO of each predicted gene to the possible Brite KEGG pathways, to produce the img_annotation.master annotation file used to produce the functional analysis and the figures.
python make_master_table.py brite2function.tab ko2brite.tab img_annotation.products > img_annotation.master

# 6. MAG classification:
# First pass MAG taxonomic classification
metawrap classify_bins -b BIN_REFINEMENT/metaWRAP_bins -o CLASSIFICATION -t 48
# More accurate classifications of the MAGs was obtained by manually aligning annotated sequences from each bin to the RefSeq database. Sequences used were: 16S, 23S (nucleotide sequences), and ribosomal proteins L1, L2, L3, L4, S1, S2, S3, S4 (amino acid sequences).



