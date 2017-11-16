#!/usr/bin/env bash

fastafile=~/igv/genomes/hg19/hg19.fa	# reference genome fasta file. You might need to delete and re-generate the index files in case of problems, getfasta will generate the index if not found
bedfile=~/Work/pipelines/data/AAGCTT.bed			# the bed file containing the list of cutsite intervals 
datadir=/home/garner1/Work/dataset/rest2vec

# ####################
# # GIVEN THE OUTPUT FROM RESTSEQ PIPELINE GENERATE A chr-start-strand-umi-quality-tag .TSV FILE
# ####################
# echo "Remove pcr duplicates [to be optimized] ..."
# wd=$PWD
# cd /home/garner1/Work/dataset/restseq/xz37/outdata
# # FILTER READ WITH LESS THAN 5 AS QUALITY
# parallel "cut -f-8 {}|LC_ALL=C grep chr -|cut -f1,2,5,6,7,8|awk '\$4>5'| 
# LC_ALL=C sort -k1,1 -k2,2n -k3,3 -k4,4nr -k5,5|awk '{print \$1,\$2,\$3,\$5,\$4,\$6}'|
# rev|LC_ALL=C uniq -f 2 | rev | tr ' ' '\t' | LC_ALL=C sort -k6,6 > {.}.tsv" ::: cutsite_dist_strand_qScore_UMI_ID_start__*_q1.bed
# cd $wd


# ####################
# # FILTER FROM THE BAM FILE ONLY THE chr-start-strand-umi-quality-tag .TSV FILES TO GET IN THE END A tag-chr-start-end-quality-strand-seq .TSV FILE
# ####################
# echo "Convert bam to bed"
# wd=$PWD
# cd /home/garner1/Work/dataset/restseq/xz37/outdata
# time parallel "bedtools bamtobed -i {} | LC_ALL=C sort -k4,4 | LC_ALL=C join -1 4 -2 6 - cutsite_dist_strand_qScore_UMI_ID_start__{.}_q1.tsv|tr ' ' '\t' > {}2bed" ::: *.bam
# echo "Convert bam to sam"
# time parallel "samtools view {} | LC_ALL=C sort -k1,1 > {.}.sam" ::: *.bam # this takes 30 min
# echo "Prepare the tsv files"
# time parallel "LC_ALL=C join -o 1.1,1.2,1.3,1.4,1.5,1.6,2.10 -1 1 -2 1 {.}.bam2bed {} | tr ' ' '\t' > {.}.tsv" ::: *.sam # takes 7min
# cd $wd

# ####################
# # PRODUCE A BED FILE WITH EXTENDED EXPERIMENTAL AND REFERENCE SEQUENCES 
# ####################
# wd=$PWD
# cd /home/garner1/Work/dataset/restseq/xz37/outdata
# parallel "cut -f2,3,4,7 {} > {.}.bed" ::: ????????.tsv # get a bed file format from experiment data
# parallel "bedtools getfasta -fi $fastafile -bed {} -bedOut > {.}.exp_ref.bed " ::: ????????.bed # get location, experimental-seq and reference-seq
# parallel "awk '{print \$1,\$2-25,\$3+25}' {} | awk '\$2 > 0'|tr ' ' '\t' > {.}.extended" ::: ????????.exp_ref.bed # enlarge the location 
# parallel "bedtools getfasta -fi $fastafile -bed {} -bedOut | awk '{print \$1,\$2,\$3,substr(\$4,1,25),substr(\$4,length(\$4)-24,length(\$4))}'|tr ' ' '\t'> {}.seq " ::: ????????.exp_ref.extended # enlarge the reference sequence
# parallel "paste {.}.extended.seq {} | awk '{print \$1,\$2,\$3,\$4\$9\$5,\$4\$10\$5}'| tr ' ' '\t' > {}.doc " ::: ????????.exp_ref.bed # enlarged locations,experimental sequence and reference sequence
# cd $wd

####################################################
####################################################
####################################################
# CONSTRUCT THE WORD EMBEDDING FROM REFERENCE GENOME
####################################################
# # PREPARE A LIST OF 200bp REGIONS CENTERED AT THE CUTSITE, IN THE + STRAND 
# echo "Prepare extended feature bedfile ..."
# awk '{OFS="\t";print $1,$2-100,$2+100,"forward","1","+"}' $bedfile | awk '$2 > 0' | #remove negative locations
# LC_ALL=C sort -k1,1 -k2,2n > $datadir/forward.bed
# bedtools merge -i $datadir/forward.bed -s | awk '{OFS="\t";print $1,$2,$3,"forward","1","+"}' > $datadir/forward.merged.bed # merge overlapping intervals
# echo "Prepare extended reference document ..."
# bedtools getfasta -fi $fastafile -bed $datadir/forward.merged.bed -bedOut -s | LC_ALL=C grep -v N > $datadir/refExtendedDoc.bed

# echo 'Prepare kmer table'
# g++ -std=c++11 ./kmers.cpp -o ./kmers # compile kmers                                                                                                                                                   
# mkdir -p "$datadir"/aux_plus
# wd=$PWD                                                                                                                                                     
# cd "$datadir/aux_plus"
# awk '{print $7 >> $1; close($1)}' $datadir/refExtendedDoc.bed # SPLIT BY CHROMOSOME
# parallel "sed -i 's/$/NNNNNN/' {}" ::: chr* #introduce the Ns to count the kmers
# cd $wd

# echo "Count kmers ..."   
# mkdir -p "$datadir"/6mer_plus           
# rm -f "$datadir"/6mer_plus/*                                                                                                                                                                          
# wd=$PWD                                                                                                                                                     
# cd "$datadir/aux_plus"
# time parallel "cat {} | tr -d '\n' | $wd/kmers 6 6 |LC_ALL=C grep -v 'N' | awk 'NF == 2' > ../6mer_plus/{}.tsv" ::: chr*                                                                            
# cd $wd

# echo "Prepare information tables ..."   
# wd=$PWD                                                                                                                                                     
# cd "$datadir/6mer_plus"
# parallel "bash $wd/parse_kmer_table.sh {} 6 1.0" ::: chr{?,??}.tsv
# rm *kmer*
# rename 's/tsv.table.tsv/table.tsv/' *tsv 
# cd $wd

# echo "Tokenize the documents"
# mkdir -p "$datadir"/docs_plus
# # SEGMENT THE GENOME
# g++ -std=c++11 ./tokenizer_withMean.cpp -o ./mean # compile the tokenizer
# wd=$PWD
# cd "$datadir"/aux_plus
# parallel "sed -i 's/NNNNNN//' {}" ::: chr*
# time parallel "$wd/mean {} ../6mer_plus/{.}.table.tsv | cut -d' ' -f2- > ../docs_plus/{.}.txt" ::: chr* # it takes approx 6 min
# cd $wd

# echo "Run word2vec on the reference corpus ..."
# time python word2vector.py $datadir/docs_plus $datadir/modelPlus 
################################
################################
################################
################################
# TOKENIZE DATA

# wd=$PWD
# cd /home/garner1/Work/dataset/restseq/xz37/outdata
# mkdir -p $datadir/restseq_plus
# for file in $( ls ????????.exp_ref.bed.doc ); do
#     barcode=`echo $file | cut -d'.' -f1`
#     echo $barcode
#     mkdir -p $barcode
#     cd $barcode
#     rm -f chr*
#     awk '{print $4 >> $1".experimental"; print $5 >> $1".reference"; print $1,$2,$3 >> $1".locations"}' /home/garner1/Work/dataset/restseq/xz37/outdata/"$barcode".exp_ref.bed.doc
#     cd ..
#     mv $barcode $datadir/restseq_plus
# done
# cd $wd

# wd=$PWD
# for dir in $( ls $datadir/restseq_plus );do
#     echo $dir
#     cd $datadir/restseq_plus/$dir
#     rm -f *.txt
#     parallel "$wd/mean {} $datadir/6mer_plus/{.}.table.tsv | cut -d' ' -f2- > {}.txt" ::: *.{experimental,reference}
#     parallel "paste -d '|' {} {.}.experimental.txt {.}.reference.txt > {.}.loc-docs" ::: chr*.locations
# done
# cd $wd

# WRITE PYTHON CODE TO PROCESS *.loc-docs files
################################################
