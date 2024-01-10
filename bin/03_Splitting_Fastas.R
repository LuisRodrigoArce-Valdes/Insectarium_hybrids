#!/usr/bin/Rscript
# 03_Splitting_Fastas.R
# Made by Luis Rodrigo Arce Valdés, to filter fasta files, choose COI sequences only and split per species
# Calling up libraries
library(phylotools)
library(dplyr)
library(tidyr)
library(stringr)
rm(list = ls())

# List of expected species per group
species <- read.table("../data/02_BOLD/Species.txt", sep = ",", header = F, col.names = "Species")
# Reading fasta files files per order
fastas <- read.fasta("../data/02_BOLD/COIs.fasta", clean_name = F)
# Filtering for COI
fastas <- fastas[grepl("COI-5P", fastas$seq.name),]
# Splitting column names
fastas <- separate(data = fastas, col = seq.name, into = c("BOLD", "Species", "Gene", "GenBank"), sep = "\\|", fill = "right")
# Removing "-" from the beginning and the end of the sequences
fastas$seq.text <- gsub("^-*","", fastas$seq.text)
fastas$seq.text <- gsub("-*$","", fastas$seq.text)
# Sorting by species names
fastas <- fastas[order(fastas$Species, decreasing = F),]
row.names(fastas) <- 1:nrow(fastas)
# Looking for missing species
missings <- setdiff(species$Species,fastas$Species)
# Writing missing species
write(missings, paste0("../data/02_BOLD/missings.txt"), ncolumns = 1)

# Splitting per species
dir.create(paste0("../data/03_Species_Fastas"), showWarnings = F)

# Filtering for species within our database only
spp <- unique(fastas[fastas$Species %in% species$Species,"Species"])
for (n in spp) {
  species <- fastas[fastas$Species==n,]
  # Transforming table to phylotools data
  species <- data.frame(seq.name=paste(species$BOLD, species$Species, species$Gene, species$GenBank, sep = "|"), seq.text=species$seq.text)
  # Writing output fastas
  n <- gsub(" ","_",n)
  phylotools::dat2fasta(species, outfile = paste0("../data/03_Species_Fastas/",n,".fasta"))
}
