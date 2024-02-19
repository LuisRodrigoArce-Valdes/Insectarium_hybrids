#!/usr/bin/Rscript
# 07_Genetic_Distancess.R
# Made by Luis Rodrigo Arce Valdés, to estimate genetic distances between COIs of hybridising species
rm(list = ls())

# Calling up libraries
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)

# Genetic Distance Model (help: dist.dna):
model <- "raw"

# Reading FASTAs
COIs <- list()
consensus <- list()
# Listing all .fasta files per group
temp <- list.files(path = "../data/07_Paired_Alignments/", pattern = "*.fasta")
# Editing names
temp <- gsub(".fasta","",temp)
# Reading all .fasta files and assigning them to their name
for (i in temp) {
  print(i)
  COIs[[i]] <- read.FASTA(paste0("../data/07_Paired_Alignments/",i,".fasta"), type = "DNA")
  consensus[[i]] <- dist.dna(COIs[[i]], model = model, variance = F, as.matrix = T)
  consensus[[i]] <- consensus[[i]][1,2]
}
consensus <- as.data.frame(t.data.frame(as.data.frame(consensus)))
consensus$Cross <- row.names(consensus)
row.names(consensus) <- 1:nrow(consensus)
consensus <- consensus[,c(2,1)]
colnames(consensus)[2] <- "Distance"

# Reading insectarium database
insectarium <- read.delim("../data/01_Insectarium_hybrids/Hybrids.txt")
insectarium %>% 
  mutate(Cross = paste0(Sp1," X ",Sp2)) -> insectarium

# Merging insectarium data to consensus distances
consensus %>% 
  mutate(Cross = gsub("_"," ",Cross)) %>% 
  left_join(insectarium, by = join_by(Cross)) %>% 
  select(Genus, Sp1, Sp2, Distance, hybrid, Condition, Authors, Year, Article) -> insectarium
rm(consensus, COIs)

# Reading reproductive barriers database
barriers <- read.csv("../data/00_raw/Lepidoptera_Barriers.csv", header = T)
barriers %>% 
  mutate(code1 = paste0(Sp1," ",Sp2)) %>%
  mutate(code2 = paste0(Sp2," ",Sp1)) -> barriers

insectarium %>%
  mutate(code1 = paste0(Sp1," ",Sp2)) %>%
  mutate(code2 = paste0(Sp2," ",Sp1)) %>%
  mutate(C1V1 = code1 %in% barriers$code1) %>%
  mutate(C2V1 = code2 %in% barriers$code1) %>%
  mutate(C1V2 = code1 %in% barriers$code2) %>%
  mutate(C2V2 = code1 %in% barriers$code2) %>%
  mutate(barriers = ifelse(C1V1 == T, T, ifelse(C1V2 == T, T, ifelse(C2V1 == T, T, ifelse(C2V2 == T, T, F))))) %>% 
  # merging dataframes using Code 1 since all three hybridazing species pairs were true on the C1V1 comparison
  left_join(barriers[,c("code1","Authors","Year","Condition","Choice.Mating","Hybrid.Egg","Hybrid.Adult","One.Sex.Fertile.Hybrids","Fertile.Hybrids")], by = join_by(code1)) %>% 
  select(Genus, Sp1, Sp2, Distance, hybrid, Condition.x, Authors.x, Year.x, Article,
         barriers, Authors.y, Year.y, Condition.y, Choice.Mating, Hybrid.Egg, Hybrid.Adult,
         One.Sex.Fertile.Hybrids, Fertile.Hybrids)-> insectarium

# Changing colnames
colnames(insectarium) <- gsub("\\.x",".hybridization", colnames(insectarium))
colnames(insectarium) <- gsub("\\.y",".barriers", colnames(insectarium))
rm(barriers)

# Replacing logical values
insectarium[insectarium==TRUE] <- "Yes"
insectarium[insectarium=="FALSE"] <- "No"

# Writing results
write.table(insectarium, "../results/01_Genetic_Distances.tsv", sep = "\t", quote = F, row.names = F)

# Plotting
ggplot(insectarium) +
  geom_violin(aes(y=Distance, x=Genus)) +
  geom_point(aes(y=Distance, x=Genus, color = hybrid, shape = barriers), size = 3, alpha = 0.75) +
  scale_color_manual(values = c("black","red")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(y="Genetic distance (% COI substitutions)", color = "Record of hybridization?", shape = "Reproductive barriers measurements") +
  theme_classic() +
  theme(legend.position = "bottom",
        text = element_text(family = "serif", size = 14),
        axis.text.x = element_text(face = "italic"),
        axis.title.x = element_blank()) -> p

pdf("../results/02_MontrealInsectarium.pdf", width = 10, height = 7)
p
dev.off()

png("../results/02_MontrealInsectarium.png", width = 10, height = 7, units = "in", res = 300)
p
dev.off()
