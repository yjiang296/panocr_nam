#!/bin/bash
sample_name=$1
PATH_to_sort.rmdup.shift.bed=$2
PATH_to_initial_narrowpeak=$3
PATH_to_genome_fai=$4
fdr=$5
##提取reads的Tn5 site
awk 'BEGIN {OFS = "\t"} ; {print $1,($2+$3)/2,($2+$3)/2+1}' $2 > $2.Tn5
##计算每个ocr的tn5 site densify（tn5 count/ ocr length)
bedtools coverage -a $3 -b $2.Tn5 | awk '{print $1,$2,$3,$4,$11}' > input_$1_peaks_and_tn5count.bed
##bedtools shuffle specifically excluding ACRs from the randomized selection space
bedtools shuffle -i $3 -g $4 |awk '{print $1,$2,$3}' | tr ' ' '\t' > $1.shuffle_randomized.bed ##fai先把scaff去掉
##计算randomized区域的tn5 density
bedtools coverage -a $1.shuffle_randomized.bed -b $2.Tn5 | awk '{print $1,$2,$3,$4}' > input_$1_randomized_and_tn5count.bed #四列 chr start end tn5count

Rscript ./filter_eFDR.R input_$1_peaks_and_tn5count.bed input_$1_randomized_and_tn5count.bed $5 $1_peaks_filted_by_tn5density.fdr$5.bed
##usage example: ./filter_ocr_by_tn5_density.sh b73 ./b73.sort.rmdup.q30.shift.bed ./b73_peaks.q005.narrowPeak ./Zm-B73-REFERENCE-NAM-5.0.fa.fai.woscaf 0.01
##会有个报错，但是不影响使用
