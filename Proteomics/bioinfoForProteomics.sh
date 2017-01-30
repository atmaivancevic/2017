#!/bin/bash

### Bioinformatics for Proteomics
### Available at https://compomics.com/bioinformatics-for-proteomics/

### 1 - Peptide and Protein Identification
# The process of searching mass spectral data can roughly be split into 6 steps:
# 1. Convert the raw, typically binary, output from the MS output into open formats
# 2. Process the MS/MS spectra into peak lists
# 3. Retrieve the desired sequence database and adapt it to the identification strategy
# 4. Search the peak lists against a sequence database using one or more search engines
# 5. Identify the peptides and infer the proteins
# 6. Validate the detected peptides and proteins

# In this chapter, we'll learn the following things:
# - Database Generation
# - Peak List Generation
# - Peptide to Spectrum Matching
# - Browsing Identification Results
# - Peptide and Protein Validation
# - PTM Analysis
# - De Novo Peptide Identification



##### 1.1 Data Generation #####
# Go to the UniProt website (http://www.uniprot.org/) and select "Taxonomy" under "Supporting data"
# Fill in "homo sapiens" (note the quotes!) in the Taxonomy filed at the top and hit enter.
# For the "Homo sapiens (Human)" result, click the link "UniProtKB"
# This shows all the human proteins on UniProt

# Select the "Reviewed" (Swiss-Prot) option in the upper left corner to only list the reviewed proteins
# Then click "Download" button above the table - make sure that "Download all" is selected, the format is "FASTA (canonical)" and "Uncompressed" is selected. Click "Go" to start downloading the protein sequences. 
# Save the file as uniprot_human_reviewed_30jan2017.fasta (in this directory)
### Handy hint: always document the database type/version and date in the file name!

# So, we have a FASTA file with all reviewed human protein sequences.
# But our sample may contain contamination from other species
# Considering sample contamination is esp. important when searching non-human data, because tiny amounts of human keratin, from hair or skin, often end up in samples.
# These need to be filtered out
# A list of common contaminants can be found at the Global Proteome Machine (http://www.thegpm.org/crap/)

# In our sample, there is one non-human protein that needs to be filtered out - the enzyme we used to digest the proteins into peptides (in this case, Trypsin).
# We will therefore add the protein seq for Trypsin to our FASTA file.
# Go back to the main UniProt website (link above) and search for Trypsin from pig: "P00761"
# Click "Sequence" in the left menu and click the FASTA download option
# Add this sequence to the bottom of uniprot_human_reviewed_30jan2017.fasta and save as uniprot_human_reviewed_trypsin_30jan2017.fasta

# Final note: if you want to create your own custom database, go to the Database Help page (http://compomics.github.io/searchgui/wiki/databasehelp.html).



#### 1.2 Peak List Generation #####
# MGF is the preferred format for spectrum identification - contains only MS/MS peak lists with some basic information about the precursor
# Most companies automatically convert the raw data files into MGF
# So I'm going to skip this step (already have the mgf files)

# Final note: Two graphical interfaces are currently available that allow you to look at your data (TOPPview) and draw pipelines (TOPPAS).



##### 1.3 Peptide to Spectrum Matching #####
# We are going to search the mgf file from section 1.2 against the human database from section 1.1
# Both files are in the 1.3-resources subdirectory
# searchgui is a graphical user interface we're going to use (download from http://compomics.github.io/projects/searchgui.html)

# Move the folder to Applications
# Start SearchGUI by double-clicking the file SearchGUI-3.2.5.jar

# To perform the search, need to provide the spectra (mgf file), database (fasta file) and experiment dependent search settings
# Load the mgf file. On the question about adding missing precursor charges choose "No"

# Next add the Search Settings. 
# - Title/Name: Tutorial
# - Click the Spectrum Matching button. 
# 	- Upload the fasta database
#	- Select "yes" when it offers to add decoy sequences
# 	- click "OK" to close the dialog

# Next, specify the modifications to consider. 
# - Fixed modification: Carbamidomethylation of C
# - Variable modification: Oxidation of M
# - Digestion: Enzyme
# - Enzyme: Trypsin
# - Specificity: Specific
# - Max missed cleavages: 2
# - Precurzor m/z Tolerance: 10 ppm
# - Fragment m/z Tolerance: 0.02 Da
# All other values default
# Click "OK" and then "Save"

# Set the output folder to the desired location. 
# This is where the search results will be stored - always use an empty folder, for simplicity
# Set to /Users/atma/2017/Proteomics/1.3-resources/searchOutput/

# Select which search engines you want
# To save time, we'll just select X! Tandem and OMSSA
# Leave the PeptideShake post-processing option unchecked for now

# Press "Start the Search!"



##### 1.4 Browsing Identification Results #####
# The search conducted above generated two files containing the peptides matched by OMSSA and X!Tandem for each spectrum, so-called Peptide to Spectrum Matches (PSMs)
# From those, want to find the identified peptides and proteins.
# This is what PeptideShaker does 
# Download from http://compomics.github.io/projects/peptide-shaker.html
(downloaded and moved to Applications)
.........this is where im up to..........

Plan of action:
All of these tools are available both as graphical interfaces and as command line options.
First, go through each tutorial using the interface, to get a feel of the steps.
Then, try to automate the pipeline using the command line.
Then, try it out with real data (Ramans mouse data).

Another good resource: http://www.proteinspector.com/pdf/oveland2015.pdf
This is a review titled "Viewing the proteome: How to visualise proteomics data?"
To generate a volcano plot, it suggests the tools GProX and Perseus







