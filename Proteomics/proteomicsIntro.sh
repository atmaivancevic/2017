#!/bin/bash

### Proteomics introduction course
### Available at https://www.ebi.ac.uk/training/online/course/proteomics-introduction-ebi-resources/what-proteomics

# Proteomics = large-scale study of proteomes
# Proteome = a set of proteins produced in an organism, system or biological context.
# The proteome is not constant - differs from cell to cell and changes over time
# So, proteomics is used to investigate:
# 	- when and where proteins are expressed
# 	- rates of protein production, degradation, and steady-state abundance
# 	- how proteins are modified (e.g. phosphorylation, other post-translational modifications)
# 	- the movement of proteins between subcellular compartment
# 	- the involvement of proteins in metabolic pathways
# 	- how proteins interact with one another

# Problems it might help us answer:
# 	- Which proteins interact with a particular protein of interest (e.g. p53 tumour suppressor protein)?
#	- Which proteins are localised to a subcellular compartment (e.g. the mitochondrion?)
#	- Which proteins are involved in a biological process (e.g. circadian rhythm)?

# The EBI (https://www.ebi.ac.uk/) hosts up-to-date databases to enable rapid searching and retreival of these data. 
# Four major databases are:
#	- UniProtKB (http://www.uniprot.org/)
#		- contains protein seqs and info avout the known biological functions
#	- IntAct (https://www.ebi.ac.uk/intact/)
#		- contains info about protein interactions
#	- Reactome (http://www.reactome.org/)
#		- contains info about which roles proteins play in human biological pathways, and which processes they contribute to
#	- PRIDE (https://www.ebi.ac.uk/pride/archive/)
#		- contains experimental evidence of published protein and peptide identifications


# Proteomics experiments typically collect data on three properties of proteins in a sample:
#	- location
#	- abundance/turnover
#	- post-translation modifications

# E.g. a typical worflow might be: 
# 	- experimental evidence analysed with UniProtKB to generate new protein identifications
# 	- the raw data is uploaded to PRIDE in support of published results
# 	- curated results demonstrating protein-protein interactions or evidence of associations with a bioligcal pathway are added to IntAct and Reactome respectively
#	- all four databases use UniProt accession numbers as a standard method of referencing proteins





