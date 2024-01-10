#!/usr/bin/Rscript
# 05_Paired_Fastas.R
# Made by Luis Rodrigo Arce Valdés, to merge into a single fasta file hybridizing species COIs
# Calling up libraries
rm(list = ls())
library(phylotools)
library(tidyr)
library(dplyr)

# List of hybridizing species pairs
species <- read.delim("../data/01_Insectarium_hybrids/Hybrids.txt", header = T)
# Keeping only the species columns
species <- species[,c("Sp1","Sp2")]
# Replacing " " by "_"
species <- as.data.frame(apply(species, 2, function(x) gsub(" ", "_", x)))
# Adding number of pair
species$pair <- 1:nrow(species)
# Tidying
species <- gather(species, "Number", "Species", 1:2)
# Removing "_" from the end of species names:
species$Species <- gsub("_$","",species$Species)

# Reading fasta files files
fastas <- list()
# Listing all .fasta files per group
temp = list.files(path = "../data/05_Consensus/", pattern = "*.fasta")
# Editing names
temp <- gsub(".fasta","",temp)

# Reading all .fasta files and assigning them to their name
for (i in temp) {
  print(i)
  fastas[[i]] <- read.fasta(paste0("../data/05_Consensus/",i,".fasta"), clean_name = F)
}
fastas <- bind_rows(fastas, .id = "Species")
fastas <- fastas[,-2]

# Now we will add each species sequence
seqs <- vector()
for (i in 1:nrow(species)) {
  sp <- species[i,"Species"]
  if (sp %in% fastas$Species) {
    fa <- fastas[fastas$Species==sp,"seq.text"]
    seqs <- append(seqs, fa)
  }
  else {
    seqs <- append(seqs, NA)
  }
}
species$COI <- seqs

# Removing species without COIs
species <- species[complete.cases(species),]

# Writing output fasta file per pair of species
species <- species[,-2]
fastas <- list()
dir.create("../data/06_Paired_Fastas", showWarnings = F)
species$pair <- as.character(species$pair)

# Creating indiviual data.frames per pair of species
for (i in unique(species$pair)) {
  fastas[[i]] <- species[species$pair==i,-1]
  colnames(fastas[[i]]) <- c("seq.name","seq.text")
  if (nrow(fastas[[i]])==2) {
    print(nrow(fastas[[i]]))
    phylotools::dat2fasta(fastas[[i]], outfile = paste0("../data/06_Paired_Fastas/",fastas[[i]]$seq.name[1], "_X_", fastas[[i]]$seq.name[2],".fasta"))
  }
}
