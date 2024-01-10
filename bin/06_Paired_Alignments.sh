#!/bin/sh
# 06_Paired_Alignments.sh
# By Luis Rodrigo Arce Valdés (28/03/22)
# With this script we will align all consensus COI sequences per pair of species

mkdir -p ../data/07_Paired_Alignments
for i in $(ls ../data/06_Paired_Fastas/)
	do
	echo "########"
	echo "$i"
	# Aligning with muscle
	muscle -align ../data/06_Paired_Fastas/$i -output ../data/07_Paired_Alignments/$i
done
