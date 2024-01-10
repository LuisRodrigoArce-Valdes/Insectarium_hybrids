#!/usr/bin/Rscript
rm(list = ls())
# 01_SpeciesPairs.R
# This script reads the butterflies species pairs in the Montreal insectarium 
# and makes pairs between species of the same genera to estimate COI genetic distances between them.
library(dplyr)
library(tidyr)

# Reading input file
butterflies <- read.csv("../data/00_raw/Butterflies.csv")

# Filtering for genera that has more than two species
butters <- data.frame()
for(i in unique(butterflies$Genus)){
  if(length(butterflies[butterflies$Genus==i,"Genus"])>1){
    butters <- rbind(butters, butterflies[butterflies$Genus==i,])
  }
}
rm(butterflies)

# Removing unidentified species
butters %>% 
  filter(species != "spp.") %>% 
  unique() %>%
  mutate(Tax = paste(Genus, species)) -> butters

# Pairing all possible combination of species
butterflies <- data.frame()
for(i in unique(butters$Genus)) {
  butters %>% 
    filter(Genus == i) -> tmp
  combn(tmp$Tax, 2) %>%
    t() %>%
    as.data.frame() %>% 
    rename(Sp1 = 1, Sp2 = 2) %>%
    mutate(Genus = i) -> tmp
  butterflies <- rbind(butterflies, tmp)
}
rm(butters, tmp)

# Reading hybrids database
hybrids <- read.csv("../data/00_raw/Lepidoptera_Hybrids.csv")

# Adding species pair string columns
hybrids %>% 
  mutate(Code1 = paste0(Sp1,"_",Sp2)) %>% 
  mutate(Code2 = paste0(Sp2,"_",Sp1)) %>%
  select(Code1, Code2, Condition, Authors, Year, Article) -> hybrids

# Searching if any pair of the insectarium butterflies have a record of hybridization
butterflies %>%
  mutate(Code1 = paste0(Sp1,"_",Sp2)) %>%
  mutate(Code2 = paste0(Sp2,"_",Sp1)) %>% 
  mutate(C1V1 = Code1 %in% hybrids$Code1) %>%
  mutate(C1V2 = Code1 %in% hybrids$Code2) %>%
  mutate(C2V1 = Code2 %in% hybrids$Code1) %>%
  mutate(C2V2 = Code2 %in% hybrids$Code2) %>%
  mutate(hybrid = ifelse(C1V1 == T, T, ifelse(C1V2 == T, T, ifelse(C2V1 == T, T, ifelse(C2V2 == T, T, F))))) %>% 
  # merging dataframes using Code 1 since all three hybridazing species pairs were true on the C1V1 comparison
  left_join(hybrids, by = join_by(Code1)) %>%
  select(Genus, Sp1, Sp2, hybrid, Condition, Authors, Year, Article) %>%
  arrange(desc(Year)) %>%
  distinct(Sp1, Sp2, .keep_all = T) -> butterflies

# Saving data frame
dir.create("../data/01_Insectarium_hybrids", showWarnings = F)
write.table(butterflies, "../data/01_Insectarium_hybrids/Hybrids.txt", sep = "\t", quote = F, row.names = F)
