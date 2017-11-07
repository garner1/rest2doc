#!/usr/bin/env bash

filein=$1
dim=$2
alpha=$3			# alpha estimator in [0,1] 

# THE CANONICAL REPRESENTATION OF JELLYFISH NEEDS TO INTRODUCE THE REVERSE COMPLEMENT OF EACH COUNTED KMER!!!

awk -v dim="$dim" '{OFS="\t";print $2,substr($1,1,dim-1),substr($1,dim,1)}' "$filein" |
LC_ALL=C sort -k2,2 > "$filein"_count-kmer-kmer.tsv # produce the count per k-1 to 1 transitions

awk -v dim="$dim" '{OFS="\t";print $2,substr($1,1,dim-1)}' "$filein" |
LC_ALL=C datamash -s -g 2 sum 1 > "$filein"_kmer-sumCount.tsv # produce the normalization for each k-1 mer

LC_ALL=C join -1 2 -2 1 -o 1.2,1.3,1.1,2.2 "$filein"_count-kmer-kmer.tsv "$filein"_kmer-sumCount.tsv |
awk -v alpha=$alpha '{OFS="\t"; print $1$2,-log(($3+alpha)/($4+alpha*4))}' > "$filein".table.tsv # use the alpha-estimator
