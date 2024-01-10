#!/bin/sh
# 04_Alignments_and_Consensus.sh
# By Luis Rodrigo Arce Valdés (28/03/22)
# With this script we will align all COI sequences for each species and then create a Consensus sequence for each one.
# To use for de novo analyses, without saved alignemts and consensus sequences
mkdir -p ../data/04_Alignments
mkdir -p ../data/05_Consensus
for n in $(ls ../data/03_Species_Fastas/)
do
	# Counting number of sequences per fasta file
	echo "########"
	echo "$n"
	grep ">" ../data/03_Species_Fastas/$n | wc -l
	
	# Aligning with muscle's "super5" algorithm
	muscle -super5 ../data/03_Species_Fastas/$n -output temp.fasta
	
	# Converting ambiguous nucleotides into " " (to avoid considering them in consensus sequence creation)
	sed 's/n/ /g' temp.fasta > ../data/04_Alignments/$n
	rm temp.fasta
	
	# Creating consesus sequence [cons only works if we have more than one sequencue]
	# Using -plurality = 0, so if there is least one sequence per position consensus sequence will have a called nucleotide
	# Using -setcase = 0 so every nucleotide will be uppercase
	if [ $(grep ">" ../data/04_Alignments/$n | wc -l) -gt 1 ]
	then
		echo "More than 1 sequence, creating a consensus"
		cons -sequence ../data/04_Alignments/$n -outseq ../data/05_Consensus/$n -name $n -plurality 0 -setcase 0 -verbose
		
		# Replacing - in some files that we couln't make consensus sequences better
		sed 's/-/N/g' ../data/05_Consensus/$n | sponge ../data/05_Consensus/$n
	else
		cp ../data/04_Alignments/$n ../data/05_Consensus/$n
		sp=">$n"
		sed -i "1s/.*/$sp/" ../data/05_Consensus/$n
	fi
done
