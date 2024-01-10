#!/bin/sh
# 01_Getting_COIs.sh
# By Luis Rodrigo Arce Valdés (05/06/22)
# Using this script we will download COI sequences from the BOLD database hybridizing pair species
# http://www.boldsystems.org/
# https://github.com/CNuge/BOLD-CLI

# 01. Creating output diretories
mkdir -p ../data/02_BOLD

# 02.- First I am doing a list of all hybridizing species
# Cut to select second column (First species).
# Tail -n +2 to remove header
cut -f2 ../data/01_Insectarium_hybrids/Hybrids.txt | tail -n +2 > tmp1

# Second species
cut -f3 ../data/01_Insectarium_hybrids/Hybrids.txt | tail -n +2 > tmp2

# Concatenating both columns, sorting, removing duplicates and spaces at the end of strings
cat tmp1 tmp2 | sort | uniq | sed 's/ *$//' > ../data/02_BOLD/Species.txt
rm tmp1 tmp2

# 03.- Searching COI sequences within the BOLD systems database
./bold-cli -output ../data/02_BOLD/COIs.fasta -query sequence -marker COI-5P -taxon ../data/02_BOLD/Species.txt
