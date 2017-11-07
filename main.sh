#!/usr/bin/env bash

fastafile=~/igv/genomes/hg19/hg19.fa	# reference genome fasta file. You might need to delete and re-generate the index files in case of problems, getfasta will generate the index if not found
bedfile=~/Work/pipelines/data/AAGCTT.bed			# the bed file containing the list of cutsite intervals 
datadir=/home/garner1/Work/dataset/rest2vec

# # PREPARE A LIST OF 75bp REGIONS STARTING AT THE CUTSITE, IN THE + AND - STRAND 
# echo "Prepare feature bedfile ..."
# time awk '{OFS="\t";print $1,$2-1,$2+74,"forward","1","+"}' $bedfile | awk '$2 > 0' > $datadir/forward.bed
# time awk '{OFS="\t";print $1,$2-72,$2+3,"reverse","1","-"}' $bedfile | awk '$2 > 0' > $datadir/reverse.bed
# time cat $datadir/forward.bed $datadir/reverse.bed > $datadir/refloc.bed
# rm -f $datadir/forward.bed $datadir/reverse.bed

# # PREPARE THE LIST OF SEQUENCES IN THE CUTSITE REGIONS
# echo "Prepare reference document ..."
# time bedtools getfasta -fi $fastafile -bed $datadir/refloc.bed -bedOut -s | LC_ALL=C grep -v N > $datadir/refdocs.bed
####################################################
# # PREPARE A LIST OF 200bp REGIONS CENTERED AT THE CUTSITE, IN THE + AND - STRAND 
# echo "Prepare extended feature bedfile ..."
# time awk '{OFS="\t";print $1,$2-100,$2+100,"forward","1","+"}' $bedfile | awk '$2 > 0' | LC_ALL=C sort -k1,1 -k2,2n > $datadir/forward.bed
# time awk '{OFS="\t";print $1,$2-100,$2+100,"reverse","1","-"}' $bedfile | awk '$2 > 0' | LC_ALL=C sort -k1,1 -k2,2n > $datadir/reverse.bed
# time bedtools merge -i $datadir/forward.bed -s | awk '{OFS="\t";print $1,$2,$3,"forward","1","+"}' > $datadir/forward.merged.bed # merge overlapping intervals
# time bedtools merge -i $datadir/reverse.bed -s | awk '{OFS="\t";print $1,$2,$3,"reverse","1","-"}' > $datadir/reverse.merged.bed # merge overlapping intervals
# cat $datadir/forward.merged.bed $datadir/reverse.merged.bed | LC_ALL=C sort -k1,1 -k2,2n > $datadir/refExtendedLocations.merged.bed
# rm -f $datadir/forward* $datadir/reverse*

# # PREPARE THE LIST OF SEQUENCES CENTER AT THE CUTSITE REGIONS
# echo "Prepare extended reference document ..."
# time bedtools getfasta -fi $fastafile -bed $datadir/refExtendedLocations.merged.bed -bedOut -s | LC_ALL=C grep -v N > $datadir/refExtendedDoc.bed

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

# echo "Prepare information tables ..."   
# wd=$PWD                                                                                                                                                     
# cd "$datadir/6mer_plus"
# parallel "bash $wd/parse_kmer_table.sh {} 6 1.0" ::: chr{?,??}.tsv
# rm *kmer*
# rename 's/tsv.table.tsv/table.tsv/' *tsv 
# cd "$datadir/6mer_minus"
# parallel "bash $wd/parse_kmer_table.sh {} 6 1.0" ::: chr{?,??}.tsv
# rm *kmer*
# rename 's/tsv.table.tsv/table.tsv/' *tsv 
# cd $wd

# echo "Tokenize the documents"
# mkdir -p "$datadir"/docs_plus
# mkdir -p "$datadir"/docs_minus
# # SEGMENT THE GENOME
# g++ -std=c++11 ./tokenizer_withMean.cpp -o ./mean # compile the tokenizer
# wd=$PWD
# cd "$datadir"/aux_plus
# parallel "sed -i 's/NNNNNN//' {}" ::: chr*
# time parallel "$wd/mean {} ../6mer_plus/{.}.table.tsv | cut -d' ' -f2- > ../docs_plus/{.}.txt" ::: chr*
# cd "$datadir"/aux_minus
# parallel "sed -i 's/NNNNNN//' {}" ::: chr*
# time parallel "$wd/mean {} ../6mer_minus/{.}.table.tsv | cut -d' ' -f2- > ../docs_minus/{.}.txt" ::: chr*
# cd $wd

# echo "Run word2vec on the reference corpus ..."
# python word2vector.py $datadir/docs_plus $datadir/modelPlus 
# python word2vector.py $datadir/docs_minus $datadir/modelMinus
# ################################

# echo "Remove pcr duplicates [to be optimized] ..."
# wd=$PWD
# cd /home/garner1/Work/dataset/restseq/xz37/outdata
# # FILTER READ WITH LESS THAN 5 AS QUALITY
# parallel "cut -f-8 {}|LC_ALL=C grep chr -|cut -f1,2,5,6,7,8|awk '\$4>5'| 
# LC_ALL=C sort -k1,1 -k2,2n -k3,3 -k4,4nr -k5,5|awk '{print \$1,\$2,\$3,\$5,\$4,\$6}'|
# rev|LC_ALL=C uniq -f 2 | rev | tr ' ' '\t' | LC_ALL=C sort -k6,6 > {.}.tsv" ::: cutsite_dist_strand_qScore_UMI_ID_start__*_q1.bed
# cd $wd

# echo "Convert bam to bed"
# wd=$PWD
# cd /home/garner1/Work/dataset/restseq/xz37/outdata
# time parallel "bedtools bamtobed -i {} | LC_ALL=C sort -k4,4 | LC_ALL=C join -1 4 -2 6 - cutsite_dist_strand_qScore_UMI_ID_start__{.}_q1.tsv|tr ' ' '\t' > {}2bed" ::: *.bam
# echo "Convert bam to sam"
# time parallel "samtools view {} | LC_ALL=C sort -k1,1 > {.}.sam" ::: *.bam # this takes 30 min
# echo "Prepare the tsv files"
# time parallel "LC_ALL=C join -o 1.1,1.2,1.3,1.4,1.5,1.6,2.10 -1 1 -2 1 {.}.bam2bed {} | tr ' ' '\t' > {.}.tsv" ::: *.sam # takes 7min
# cd $wd

# mkdir -p $datadir/restseq_plus
# mkdir -p $datadir/restseq_minus
# rm -fr $datadir/restseq_{plus,minus}/*
# wd=$PWD
# cd /home/garner1/Work/dataset/restseq/xz37/outdata
# parallel "mkdir -p $datadir/restseq_plus/{.} && cd $datadir/restseq_plus/{.} && awk '\$6 == \"+\"' /home/garner1/Work/dataset/restseq/xz37/outdata/{} | awk '{print \$0 >> \$2}'" ::: ????????.tsv
# parallel "mkdir -p $datadir/restseq_minus/{.} && cd $datadir/restseq_minus/{.} && awk '\$6 == \"-\"' /home/garner1/Work/dataset/restseq/xz37/outdata/{} | awk '{print \$0 >> \$2}'" ::: ????????.tsv

# MINUS STRANDS NEEDS TO HAVE READ REV-COMPLEMENTED
# cd $datadir/restseq_minus
# parallel "cut -f7 {} | rev | tr 'ATGC' 'TACG' > {}.rc" ::: */chr{?,??}
# parallel "paste {} {}.rc | cut -f-6,8 > {}.pasted" ::: */chr{?,??}
# parallel "mv {}.pasted {} && rm {}.rc " ::: */chr{?,??}
# cd $wd

# echo "Bin the genome ..."
# window=1000000
# genome=hg19
# sliding=500000
# bash ~/Work/pipelines/aux.scripts/fetchChromSizes.sh $genome | grep -v "_" > sizes
# bedtools makewindows -g sizes -w "$window" -s "$sliding" -i winnum > $datadir/"$window"_"$genome"
# rm -f sizes

# echo "Create the cutsite documents on the reference genome ..."
# bedtools intersect -a $datadir/1000000_hg19 -b $datadir/refdocs.bed -wa -wb > $datadir/binned_refdocs.tsv
# mkdir -p $datadir/reference_docs
# wd=$PWD
# cd $datadir/reference_docs
# rm -f $datadir/reference_docs/*
# awk '{print $11 >> $1"_"$4"_"$10}' $datadir/binned_refdocs.tsv
# cd $wd

# echo "Create the cutsite documents on the restseq dataset ..."
# parallel "cat {} |awk '{OFS=\"\t\"; print \$2,\$3,\$4,\$5,\$6,\$7,\$1}'|bedtools intersect -a $datadir/1000000_hg19 -b - -wa -wb|cut -f1,4,10 > {}.binned" ::: $datadir/restseq_{plus,minus}/*/chr{?,??}

# for dir in $( ls $datadir/restseq_plus );do
#     echo $dir
#     cd $datadir/restseq_plus/$dir
#     parallel "awk '{print \$3 >> \$1\"_\"\$2}' {}" ::: *.binned
# done
# for dir in $( ls $datadir/restseq_minus );do
#     echo $dir
#     cd $datadir/restseq_minus/$dir
#     parallel "awk '{print \$3 >> \$1\"_\"\$2}' {}" ::: *.binned
# done

#NOW YOU NEED TO TOKENIZE THE BINNED DOCUMENTS AND PREPARE THE WORD COUNTS
