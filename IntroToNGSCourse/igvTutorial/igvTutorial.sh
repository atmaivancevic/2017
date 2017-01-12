#!/bin/bash

##### INTEGRATIVE GENOMICS VIEWER TUTORIAL
##### Available at https://wikis.utexas.edu/display/bioiteam/Integrative+Genomics+Viewer+(IGV)+tutorial

##### Learning objectives
# Create a custom genome database (e.g. for microbial genomes)
# Load a pre-existing genome assembly (e.g. for eukaryotes and other higher organisms)
# Load output from mapping reads to a reference genome
# Load output from calling genetic variants
# Navigate the view of the genome and interpret the display of this data

### Viewing E. coli data in IGV
# First, copy the genome and read files from Mapping Tutorial to this dir
cd /scratch/atmaGenomes/2017
mkdir igvTutorial
cd igvTutorial
cp ../mappingTutorial/NC_012967.1.gbk .
cp ../mappingTutorial/SRR030257_*.fastq .

# Next, IGV prefers the ref genome in gff format
# bp_seqconvert.pl does not have a gff option
# can convert it with a personal script, or use a different tool (e.g. Readseq)
# Download Readseq
wget http://iubio.bio.indiana.edu/soft/molbio/readseq/java/readseq.jar

# Use a general java help command
java -jar readseq.jar
java -cp readseq.jar run

# Convert the ref genome from genbank to gff format
java -cp readseq.jar run NC_012967.1.gbk -f GFF -o NC_012967.1.gbk.gff






