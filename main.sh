#!/usr/bin/env bash

fastafile=~/igv/genomes/hg19/hg19.fa	# reference genome fasta file. You might need to delete and re-generate the index files in case of problems, getfasta will generate the index if not found
bedfile=~/Work/pipelines/data/AAGCTT.bed			# the bed file containing the list of cutsite intervals 
datadir=/home/garner1/Work/dataset/rest2vec
restdata=/home/garner1/Work/dataset/restseq/xz37/outdata

####################
# GIVEN THE OUTPUT FROM RESTSEQ PIPELINE GENERATE A chr-start-strand-umi-quality-tag .TSV FILE
####################
echo "Remove pcr duplicates [to be optimized] ..."
wd=$PWD
cd $restdata
# FILTER READ WITH LESS THAN 5 AS QUALITY
parallel "cut -f-8 {}|LC_ALL=C grep chr -|cut -f1,2,5,6,7,8|awk '\$4>5'| 
LC_ALL=C sort -k1,1 -k2,2n -k3,3 -k4,4nr -k5,5|awk '{print \$1,\$2,\$3,\$5,\$4,\$6}'|
rev|LC_ALL=C uniq -f 2 | rev | tr ' ' '\t' | LC_ALL=C sort -k6,6 > {.}.tsv" ::: cutsite_dist_strand_qScore_UMI_ID_start__*_q1.bed
cd $wd


####################
# FILTER FROM THE BAM FILE ONLY THE chr-start-strand-umi-quality-tag .TSV FILES TO GET IN THE END A tag-chr-start-end-quality-strand-seq .TSV FILE
####################
echo "Convert bam to bed"
wd=$PWD
cd $restdata
time parallel "bedtools bamtobed -i {} | LC_ALL=C sort -k4,4 | LC_ALL=C join -1 4 -2 6 - cutsite_dist_strand_qScore_UMI_ID_start__{.}_q1.tsv|tr ' ' '\t' > {}2bed" ::: *.bam
echo "Convert bam to sam"
time parallel "samtools view {} | LC_ALL=C sort -k1,1 > {.}.sam" ::: *.bam # this takes 30 min
echo "Prepare the tsv files"
time parallel "LC_ALL=C join -o 1.1,1.2,1.3,1.4,1.5,1.6,2.10 -1 1 -2 1 {.}.bam2bed {} | tr ' ' '\t' > {.}.tsv" ::: *.sam # takes 7min to get location and sequence
cd $wd

####################
# PRODUCE A BED FILE WITH EXTENDED EXPERIMENTAL AND REFERENCE SEQUENCES 
####################
wd=$PWD
cd $restdata
parallel "cut -f2,3,4,7 {} > {.}.bed" ::: ????????.tsv # get a bed file format from experiment data
parallel "bedtools getfasta -fi $fastafile -bed {} -bedOut > {.}.exp_ref.bed " ::: ????????.bed # get location, experimental-seq and reference-seq
parallel "awk '{print \$1,\$2-25,\$3+25}' {} | awk '\$2 > 0'|tr ' ' '\t' > {.}.extended" ::: ????????.exp_ref.bed # enlarge the location 
parallel "bedtools getfasta -fi $fastafile -bed {} -bedOut | awk '{print \$1,\$2,\$3,substr(\$4,1,25),substr(\$4,length(\$4)-24,length(\$4))}'|tr ' ' '\t'> {}.seq " ::: ????????.exp_ref.extended # enlarge the reference sequence
parallel "paste {.}.extended.seq {} | awk '{print \$1,\$2,\$3,\$4\$9\$5,\$4\$10\$5}'| tr ' ' '\t' > {}.doc " ::: ????????.exp_ref.bed # enlarged locations,experimental sequence and reference sequence
cd $wd

###################################
# CONSTRUCT THE WORD EMBEDDING FROM REFERENCE GENOME
bash process_ref_genome.sh $bedfile $datadir $fastafile
#####################################

# TOKENIZE DATA

wd=$PWD
cd $restdata
mkdir -p $datadir/restseq_plus
for file in $( ls ????????.exp_ref.bed.doc ); do
    barcode=`echo $file | cut -d'.' -f1`
    echo $barcode
    mkdir -p $barcode
    cd $barcode
    rm -f chr*
    awk '{print $4 >> $1".experimental"; print $5 >> $1".reference"; print $1,$2,$3 >> $1".locations"}' $restdata/"$barcode".exp_ref.bed.doc
    cd ..
    mv $barcode $datadir/restseq_plus
done
cd $wd

wd=$PWD
for dir in $( ls $datadir/restseq_plus );do
    echo $dir
    cd $datadir/restseq_plus/$dir
    rm -f *.txt
    parallel "$wd/mean {} $datadir/6mer_plus/{.}.table.tsv | cut -d' ' -f2- > {}.txt" ::: *.{experimental,reference}
    parallel "paste -d '|' {} {.}.experimental.txt {.}.reference.txt | LC_ALL=C sort -k2,2n > {.}.loc-docs" ::: chr*.locations
done
cd $wd

wd=$PWD
for dir in $( ls $datadir/restseq_plus );do
    echo $dir
    cd $datadir/restseq_plus/$dir
    parallel "python $wd/compare-experiment2reference.py {}" ::: chr{?,??}.loc-docs
    parallel "paste {} {.}.loc-docs.signal | tr ' ' '\t' > {.}.signal.bed" ::: chr{?,??}.locations
done
cd $wd

echo "Bin the genome ..."
window=10000000
genome=hg19
sliding=1000000
bash ~/Work/pipelines/aux.scripts/fetchChromSizes.sh $genome | grep -v "_" > sizes
bedtools makewindows -g sizes -w "$window" -s "$sliding" -i winnum > $datadir/"$window"_"$genome"
rm -f sizes

echo "Create the cutsite documents on the reference genome ..."
wd=$PWD
for dir in $( ls $datadir/restseq_plus );do
    echo $dir
    cd $datadir/restseq_plus/$dir
    parallel "grep -v 'nan' {} | sponge {}" ::: chr*.signal.bed
    parallel "bedtools intersect -a $datadir/"$window"_"$genome" -b {} -wa -wb | tr '.' ','|datamash -s -g 1,2,3 sum 8 | tr ',' '.' > {}.binned" ::: chr*.signal.bed
done
cd $wd
