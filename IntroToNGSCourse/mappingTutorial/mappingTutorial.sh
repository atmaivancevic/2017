#!/bin/bash

##### MAPPING TUTORIAL
##### Available at https://wikis.utexas.edu/display/bioiteam/Mapping+tutorial

##### Learning Objectives:
# Learn how to index a reference genome, map reads, convert output to SAM format for downstream analysis
# Compare bowtie, bwa and bowtie2 on a publicly available E. coli Illumina data set.

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

##### USING BOWTIE
### Make a new directory for bowtie output
mkdir bowtie

### Check which version of bowtie is installed
bowtie --version
#bowtie version 1.0.0
#64-bit

### Index the reference genome file and put the output into the bowtie dir
# Note: indexing can take a long time on large genomes
# Usually have it as a separate step/script
# because you only need to do it once, then you can keep using the same index for future mappings
bowtie-build NC_012967.1.fasta bowtie/NC_012967.1

### Map the reads to the reference to generate SAM output
# Make sure to submit this as a job to the queue rather than running on head
# -t keeps track of the time taken to perform each step
# -p enables multi-threading (set to 8 cores)
bowtie -p 8 -t --sam bowtie/NC_012967.1 -1 SRR030257_1.fastq -2 SRR030257_2.fastq bowtie/SRR030257.sam
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

##### USING BWA
### Make a new dir for bwa output
mkdir bwa

### Check which version of bwa is installed
bwa 
# Version: 0.7.4-r396-beta

### As before, index ref genome and map reads to genome
# Note: bwa doesn't give you a choice of where to put output files
# It creates them in the same dir as the FASTA file that you input
# To get around this, first copy the ref genome to bwa directory
cp NC_012967.1.fasta bwa 

# Then index the genome
bwa index bwa/NC_012967.1.fasta

# Then run the mapping command (aln)
# Note: bwa separates the initial mapping, so need to map each set of reads in the pair separately
# Another note: last time, used 8 cores. Now, using 4 cores per process. 
# If running the commands one by one, can use 8 cores.
# If you were to run these together in a script, 
# and submit it to the queue on a system which only has 8 cores per node,
# then you would need to split it as follows with 4 cores for each command.
bwa aln -t 4 -f bwa/SRR030257_1.sai bwa/NC_012967.1.fasta SRR030257_1.fastq
bwa aln -t 4 -f bwa/SRR030257_2.sai bwa/NC_012967.1.fasta SRR030257_2.fastq
# *.sai output file contains the 'alignment seeds' in a format specific to bwa.
# It's basically an intermediate file
# Still need one more step to convert this to SAM output format
bwa sampe -f bwa/SRR030257.sam bwa/NC_012967.1.fasta bwa/SRR030257_1.sai bwa/SRR030257_2.sai SRR030257_1.fastq SRR030257_2.fastq

##### USING BOWTIE2
### Make a fresh dir, index genome and map reads
# bowtie2 is currently considered the latest and greatest mapper
# Best compromise between configurability, sensistivity and speed
mkdir bowtie2
bowtie2 --version
# /usr/local/bin/bowtie2-align-s version 2.2.3
# 64-bit
bowtie2-build NC_012967.1.fasta bowtie2/NC_012967.1
bowtie2 -p 8 -t -x bowtie2/NC_012967.1 -1 SRR030257_1.fastq -2 SRR030257_2.fastq -S bowtie2/SRR030257.sam
