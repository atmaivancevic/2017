#!/bin/bash

# $1 = Genus (e.g. Homo)
# $2 = Species (e.g. sapiens)
# $3 = ncbi taxid (e.g. 9606)

#echo "Genus="$1
#echo "Species="$2
#echo "NCBItxid="$3

awk 'BEGIN { print "'$1' '$2',Species information," " <h3><i>'$1' '$2'</i></h3>" " <h4>NCBI TaxID: <a target=\047_blank\047 href=\047http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id='$3'\047>'$3'</a></h4>" " <h4>Genome assembly: <a target=\047_blank\047 href=\047https://www.ncbi.nlm.nih.gov/assembly/?term='$1'+'$2'\047>available</a></h4>" " <h5>More information: </h5><ul><li><a target=\047_blank\047 href=\047https://www.google.com.au/search?q='$1'+'$2'\047>Google search</a></li><li><a target=\047_blank\047 href=\047https://en.wikipedia.org/wiki/'$1'_'$2'\047>Wikipedia search</a></li></ul>" " <img src=\047INSERTIMAGELINK\047 width=200>"}'
