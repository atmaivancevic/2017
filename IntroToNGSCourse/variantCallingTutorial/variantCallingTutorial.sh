#!/bin/bash

##### VARIANT CALLING TUTORIAL
##### Available at https://wikis.utexas.edu/display/bioiteam/Variant+calling+tutorial

##### Learning objectives


### Calling variants in reads mapped by bowtie
# First, copy the output files from Mapping Tutorial to this dir
# (working on leeuwenhoek)
cd /scratch/atmaGenomes/2017
mkdir variantCallingTutorial
cd variantCallingTutorial
cp -r ../mappingTutorial/bowtie .
cp -r ../mappingTutorial/bwa .
cp -r ../mappingTutorial/bowtie2 .

# Check samtools version and where it is installed
samtools
# Version: 0.1.18 (r982:295)
which samtools
# /usr/local/bin/samtools

# Create new output directory and copy in necessary files
mkdir samtoolsBowtie
cp bowtie/SRR030257.sam samtoolsBowtie/
cp ../mappingTutorial/NC_012967.1.fasta samtoolsBowtie/

# Use samtools to index the reference genome file
samtools faidx samtoolsBowtie/NC_012967.1.fasta

# Convert from SAM to BAM format
samtools view -b -S -o samtoolsBowtie/SRR030257.bam samtoolsBowtie/SRR030257.sam

# Sort and index the BAM file
# Note: when sorting, samtools appends an extra .bam to the end of the output
samtools sort samtoolsBowtie/SRR030257.bam samtoolsBowtie/SRR030257.sorted
samtools index samtoolsBowtie/SRR030257.sorted.bam
# Another note: don't bother ever gzipping BAM files - they're already the compressed version of SAM files

# Now we want to call genome variants
# Use the samtools mpileup command to compile info about the bases mapped to each ref position
# Output is bcf format, which is the binary form of the text Variant Call Format (vcf)
samtools mpileup -u -f samtoolsBowtie/NC_012967.1.fasta samtoolsBowtie/SRR030257.sorted.bam > samtoolsBowtie/SRR030257.bcf

# Convert bcf to human-readable vcf
bcftools view -v -c -g samtoolsBowtie/SRR030257.bcf > samtoolsBowtie/SRR030257.vcf

### Do the same with bwa and bowtie2
# bwa
mkdir samtoolsBwa
cp bwa/SRR030257.sam samtoolsBwa/
cp ../mappingTutorial/NC_012967.1.fasta samtoolsBwa/
samtools view -b -S -o samtoolsBwa/SRR030257.bam samtoolsBwa/SRR030257.sam 
# this step fails...
# Parse error at line 7: missing colon in auxiliary data
# Aborted (core dumped)
# > sed -n 7p samtoolsBwa/SRR030257.sam # to see line 7
#SRR030257.3	83	NC_012967	228261	0	36M	=	228139	-158	TGAGGTCGGTGGTTCAAGTCCACTCAGGCCTACCAA	;<><>-???<?AA>>AAAAA>AAAAAAAAAAAAAAA	XT:A:R	NM:i:0	SM:i:0	AM:i:0	X0:i:3	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:36	XA:Z:
# It doesn't look like it finished....bwa must have done something wrong when generating the SAM output

# bowtie2
mkdir samtoolsBowtie2
cp bowtie2/SRR030257.sam samtoolsBowtie2/
cp ../mappingTutorial/NC_012967.1.fasta samtoolsBowtie2/
samtools view -b -S -o samtoolsBowtie2/SRR030257.bam samtoolsBowtie2/SRR030257.sam 
samtools sort samtoolsBowtie2/SRR030257.bam samtoolsBowtie2/SRR030257.sorted
samtools index samtoolsBowtie2/SRR030257.sorted.bam 
samtools mpileup -u -f samtoolsBowtie2/NC_012967.1.fasta samtoolsBowtie2/SRR030257.sorted.bam > samtoolsBowtie2/SRR030257.bcf
bcftools view -v -c -g samtoolsBowtie2/SRR030257.bcf > samtoolsBowtie2/SRR030257.vcf
# done, worked perfectly!

### Filtering VCF files with grep
# vcf includes alternative Allele Frequency tags, denoted AF1
# Find what values of AF1 are in the vcf files
grep AF1 samtoolsBowtie2/SRR030257.vcf
grep AF1 samtoolsBowtie/SRR030257.vcf
# Note: the E.coli genome is haploid - there aren't any geterozygotes
# So predictions with an allele freq not equal to 1 are not really applicable
# We should really remove these lines (where AF1=0) from the file
# e.g. cat input.vcf | grep -v AF1=0 > output.vcf

### Comparing the results of different mappers using BEDTools
# Make a new output dir and copy all vcf files to it, renaming them informatively
mkdir comparison
cp samtoolsBowtie/SRR030257.vcf comparison/bowtie.vcf
cp samtoolsBowtie2/SRR030257.vcf comparison/bowtie2.vcf
cd comparison

# See which bedtools version is available
bedtools --version
# bedtools v2.22.0

# Find common mutations between the files
bedtools intersect -a bowtie.vcf -b bowtie2.vcf > commonBowtieBowtie2.vcf
# 34 total

# Find mutations that are unique for each mapper
bedtools subtract -a bowtie.vcf -b commonBowtieBowtie2.vcf > uniqueBowtie.vcf
# 186 total
bedtools subtract -a bowtie2.vcf -b commonBowtieBowtie2.vcf > uniqueBowtie2.vcf
# 94 total

# So Bowtie found more variants than Bowtie2. Huh.




