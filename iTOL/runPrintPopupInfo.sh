#!/bin/bash

# Use this file to print popup HTML code for all genomes
# Takes as input a text file containing species names and NCBI txid number
# E.g. input txt file looks like:
#Macropus eugenii 9393
#Homo sapiens 9606
#...


while IFS= read -r line; do
  echo "reading input: $line"
./printPopupInfo.sh $line >> htmlForPopup.txt
echo "###" >> htmlForPopup.txt
done < genomesAndTxidNos.txt

