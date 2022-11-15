#!/usr/bin/python
import pandas as pd
import sys
query = sys.argv[1]
file = pd.read_csv(query+'q10.minim.bam.detail',sep='\t')  ###can be modified

def cal_coverage_and_identity(df):
    identity = 1 - df['NM'].sum()/df['MID_length'].sum()
    unionSet = [0 for i in range(10000)]
    for index,row in df.iterrows():
        start = row['Query_aligned_start']
        end = row['Query_aligned_end']
        unionSet[start:end] = [1 for i in range(end-start)]
    cover = sum(unionSet)
    coverage = cover / df['Query_length'].mean()

    return coverage, identity, coverage < 0.9 or identity < 0.9


groups = file.groupby('Query_name')
res = groups.apply(cal_coverage_and_identity)
res.to_csv('./calculate_out/'+query+'.cov_or_ide.lessthan90.csv')

##### usage: python calculate_bam_detail_coverage_identity.py *.bam.detail  #####
##### author: yjiang last_update: 2022-09-07######
'''
This program is used to calculate a comprehensive coverage and identity of query seqs that have secondary alignments.
*bam.detail is the output of get_the_detailed_mapping_information_from_minimap2_bam_jymodified.pl.

'''