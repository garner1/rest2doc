#!/usr/bin/env bash

fastafile=~/igv/genomes/hg19/hg19.fa	# reference genome fasta file. You might need to delete and re-generate the index files in case of problems, getfasta will generate the index if not found
bedfile=~/Work/pipelines/data/CATG.bed			# the bed file containing the list of cutsite intervals 
datadir=/home/garner1/Work/dataset/rest2vec

# PREPARE A LIST OF 75bp REGIONS STARTING AT THE CUTSITE, IN THE + AND - STRAND 
echo "Prepare feature bedfile ..."
time awk '{OFS="\t";print $1,$2-1,$2+74,"forward","1","+"}' $bedfile | awk '$2 > 0' > $datadir/forward.bed
time awk '{OFS="\t";print $1,$2-72,$2+3,"reverse","1","-"}' $bedfile | awk '$2 > 0' > $datadir/reverse.bed
time cat $datadir/forward.bed $datadir/reverse.bed > $datadir/refloc.bed
rm -f $datadir/forward.bed $datadir/reverse.bed

# PREPARE THE LIST OF SEQUENCES IN THE CUTSITE REGIONS
# echo "Prepare reference document ..."
# time bedtools getfasta -fi $fastafile -bed $datadir/refloc.bed -bedOut -s | LC_ALL=C grep -v N > $datadir/refdocs.bed

# PREPARE A LIST OF 200bp REGIONS CENTERED AT THE CUTSITE, IN THE + AND - STRAND 
# echo "Prepare extended feature bedfile ..."
# time awk '{OFS="\t";print $1,$2-100,$2+100,"forward","1","+"}' $bedfile | awk '$2 > 0' | LC_ALL=C sort -k1,1 -k2,2n > $datadir/forward.bed
# time awk '{OFS="\t";print $1,$2-100,$2+100,"reverse","1","-"}' $bedfile | awk '$2 > 0' | LC_ALL=C sort -k1,1 -k2,2n > $datadir/reverse.bed
# time bedtools merge -i $datadir/forward.bed -s | awk '{OFS="\t";print $1,$2,$3,"forward","1","+"}' > $datadir/forward.merged.bed # merge overlapping intervals
# time bedtools merge -i $datadir/reverse.bed -s | awk '{OFS="\t";print $1,$2,$3,"reverse","1","-"}' > $datadir/reverse.merged.bed # merge overlapping intervals
# cat $datadir/forward.merged.bed $datadir/reverse.merged.bed | LC_ALL=C sort -k1,1 -k2,2n > $datadir/refExtendedLocations.merged.bed
# rm -f $datadir/forward* $datadir/reverse*

# PREPARE THE LIST OF SEQUENCES CENTER AT THE CUTSITE REGIONS
# echo "Prepare extended reference document ..."
# time bedtools getfasta -fi $fastafile -bed $datadir/refExtendedLocations.merged.bed -bedOut -s | LC_ALL=C grep -v N > $datadir/refExtendedDoc.bed

#####                                                                                                                                                                                                      
# echo 'Prepare kmer table'
# g++ -std=c++11 ./kmers.cpp -o ./kmers # compile kmers                                                                                                                                                   
# mkdir -p "$datadir"/aux_plus "$datadir"/aux_minus
# wd=$PWD                                                                                                                                                     
# cd "$datadir/aux_plus"
# awk '$6 == "+"' ../refExtendedDoc.bed | awk '{print $7 >> $1; close($1)}' - # SPLIT BY CHROMOSOME
# parallel "sed -i 's/$/NNNNNN/' {}" ::: chr*
# cd "$datadir/aux_minus"
# awk '$6 == "-"' ../refExtendedDoc.bed | awk '{print $7 >> $1; close($1)}' -
# parallel "sed -i 's/$/NNNNNN/' {}" ::: chr*
# cd $wd

# echo "Count kmers ..."   
# mkdir -p "$datadir"/6mer_{plus,minus}           
# rm -f "$datadir"/6mer_*/*                                                                                                                                                                          
# wd=$PWD                                                                                                                                                     
# cd "$datadir/aux_plus"
# parallel "cat {} | tr -d '\n' | $wd/kmers 6 6 |LC_ALL=C grep -v 'N' | awk 'NF == 2' > ../6mer_plus/{}.tsv" ::: chr*                                                                            
# cd "$datadir/aux_minus"
# parallel "cat {} | tr -d '\n' | $wd/kmers 6 6 |LC_ALL=C grep -v 'N' | awk 'NF == 2' > ../6mer_minus/{}.tsv" ::: chr*                                                                            
# cd $wd

############################
# NOW YOU NEED TO TOKENIZE THE DOCS AND THEN RUN WORD2VEC
# ####    

# echo "Convert bam to bed"
# time parallel "bedtools bamtobed -i {} | LC_ALL=C sort -k4,4 > {}2bed" ::: /home/garner1/Work/dataset/restseq/xz37/outdata/*.bam
# echo "Convert bam to sam"
# time parallel "samtools view {} | LC_ALL=C sort -k1,1 > {.}.sam" ::: /home/garner1/Work/dataset/restseq/xz37/outdata/*.bam

# echo "Prepare the tsv files"
# time parallel "LC_ALL=C join -o 1.1,1.2,1.3,1.4,1.5,1.6,2.10 -1 4 -2 1 {.}.bam2bed {} | tr ' ' '\t' > {.}.tsv" ::: /home/garner1/Work/dataset/restseq/xz37/outdata/*.sam
# mv /home/garner1/Work/dataset/restseq/xz37/outdata/*.{sam,tsv,bam2bed} $datadir/docs

