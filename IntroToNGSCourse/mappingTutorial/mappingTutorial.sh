#!/bin/bash

##### Mapping Tutorial 
##### Available at https://wikis.utexas.edu/display/bioiteam/Mapping+tutorial

### First, download sample reference genome and paired end read data files
# E.coli reference genome on https://www.ncbi.nlm.nih.gov/nuccore/NC_012967.1
# Send -> to File -> Format: GenBank (full) -> Show GI -> Create File

# Paired-end read files (2 total): http://www.ebi.ac.uk/ena/data/view/SRR030257&display=html
# Download File 1 and File 2  

### Upload to leeuwenhoek
# (directory /scratch/atmaGenomes/2017/mappingTutorial/)
# Should have 3 files: NC_012967.1.gbk, SRR030257_1.fastq, SRR030257_2.fastq

### Go to the right directory
cd /scratch/atmaGenomes/2017/mappingTutorial/

### Count the number of sequences in file SRR030257_1.fastq
grep -c ^@SRR030257 SRR030257_1.fastq 

### Count how many bases long the reads in SRR030257_1.fastq are
sed -n 2p SRR030257_1.fastq | awk -F"[ATCGatcg]" '{print NF-1}'

### Convert between file formats
# bp_seqconvert.pl is an in-built bioperl script that converts between diff seq formats
# e.g. to convert the reference genome from genbank to fasta format:
bp_seqconvert.pl --from genbank --to fasta < NC_012967.1.gbk > NC_012967.1.fasta

### Make a new directory for bowtie output
mkdir bowtie

### Index the reference genome file and put the output into the bowtie dir
# Note: indexing can take a long time on large genomes
# Usually have it as a separate step/script
# because you only need to do it once, then you can keep using the same index for future mappings
bowtie-build NC_012967.1.fasta bowtie/NC_012967.1

### Map the reads to the reference to generate SAM output
# Make sure to submit this as a job to the queue rather than running on head
# -t keeps track of the time taken to perform each step
bowtie -t --sam bowtie/NC_012967.1 -1 SRR030257_1.fastq -2 SRR030257_2.fastq bowtie/SRR030257.sam
# Time loading reference: 00:00:00
# Time loading forward index: 00:00:00
# Time loading mirror index: 00:00:00
# Seeded quality full-index search: 00:05:29
# reads processed: 3800180
# reads with at least one reported alignment: 3500240 (92.11%)
# reads that failed to align: 299940 (7.89%)
# Reported 3500240 paired-end alignments to 1 output stream(s)
# Time searching: 00:05:29
# Overall time: 00:05:29

### Do a sanity check of the SAM output
# Note: Don't load the entire file into memory - too big! Always head or grep SAM file
head bowtie/SRR030257.sam
