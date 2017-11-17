#!/usr/bin/env bash

bedfile=$1
datadir=$2
fastafile=$3

# CONSTRUCT THE WORD EMBEDDING FROM REFERENCE GENOME
# PREPARE A LIST OF 200bp REGIONS CENTERED AT THE CUTSITE, IN THE + STRAND 

echo "Prepare extended feature bedfile ..."
awk '{OFS="\t";print $1,$2-100,$2+100,"forward","1","+"}' $bedfile | #enlarge the boundary
awk '$2 > 0' | #remove negative locations
LC_ALL=C sort -k1,1 -k2,2n > $datadir/forward.bed
bedtools merge -i $datadir/forward.bed -s | awk '{OFS="\t";print $1,$2,$3,"forward","1","+"}' > $datadir/forward.merged.bed # merge overlapping intervals
echo "Prepare extended reference document ..."
bedtools getfasta -fi $fastafile -bed $datadir/forward.merged.bed -bedOut -s | LC_ALL=C grep -v N > $datadir/refExtendedDoc.bed

echo 'Prepare kmer table'
g++ -std=c++11 ./kmers.cpp -o ./kmers # compile kmers                                                                                                                                                   
mkdir -p "$datadir"/aux_plus
wd=$PWD                                                                                                                                                     
cd "$datadir/aux_plus"
awk '{print $7 >> $1; close($1)}' $datadir/refExtendedDoc.bed # split by chromosome
parallel "sed -i 's/$/NNNNNN/' {}" ::: chr* #introduce the Ns to count the kmers
cd $wd

echo "Count kmers ..."   
mkdir -p "$datadir"/6mer_plus           
rm -f "$datadir"/6mer_plus/*                                                                                                                                                                          
wd=$PWD                                                                                                                                                     
cd "$datadir/aux_plus"
time parallel "cat {} | tr -d '\n' | $wd/kmers 6 6 |LC_ALL=C grep -v 'N' | awk 'NF == 2' > ../6mer_plus/{}.tsv" ::: chr*                                                                            
cd $wd

echo "Prepare information tables ..."   
wd=$PWD                                                                                                                                                     
cd "$datadir/6mer_plus"
parallel "bash $wd/parse_kmer_table.sh {} 6 1.0" ::: chr{?,??}.tsv
rm *kmer*
rename 's/tsv.table.tsv/table.tsv/' *tsv 
cd $wd

echo "Tokenize the documents"
mkdir -p "$datadir"/docs_plus
# SEGMENT THE GENOME
g++ -std=c++11 ./tokenizer_withMean.cpp -o ./mean # compile the tokenizer
wd=$PWD
cd "$datadir"/aux_plus
parallel "sed -i 's/NNNNNN//' {}" ::: chr*
time parallel "$wd/mean {} ../6mer_plus/{.}.table.tsv | cut -d' ' -f2- > ../docs_plus/{.}.txt" ::: chr* # it takes approx 6 min
cd $wd

echo "Run word2vec on the reference corpus ..."
time python word2vector.py $datadir/docs_plus $datadir/modelPlus 
