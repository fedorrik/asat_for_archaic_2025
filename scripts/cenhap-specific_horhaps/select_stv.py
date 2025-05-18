# script takes stv.bed annotation (path to file as argument) and for each cenhap (last number in chrom after "_")
# selects samples of number_to_sample of all stvs with frequncy >= frequency_threshold
from sys import argv
import pandas as pd
from collections import Counter

frequency_threshold = 0.05
number_to_sample = int(argv[2])
stv_annotation_path = argv[1]

# parse input bed
bed_header = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'rgb']
stv = pd.read_csv(argv[1], sep='\t', names=bed_header)

# remove stv numbering
stv['name'] = stv['name'].apply(lambda x: x.split(':')[1])

# '/' to h in hybrid names
stv['name'] = stv['name'].apply(lambda x: x.replace('/', 'h'))

# add cenhap assignment
stv['cenhap'] = stv['chrom'].apply(lambda x: x.split('_')[-1])

# iterate through cenhaps
for cenhap, stvs in stv.groupby('cenhap')['name']:
    # choose frequent stvs
    frequent_stvs = {stv: cnt for stv, cnt in Counter(stvs).most_common() if cnt / len(stvs) >= frequency_threshold}
    # get dataFrame of cenhap
    cenhap_df = stv[stv['cenhap'] == cenhap]
    # iterate through frequent stvs
    for frequent_stv in frequent_stvs:
        # select frequent stvs in cenhap dataFrame
        stv_df = cenhap_df[cenhap_df['name'] == frequent_stv]
        # make a sample if there are a lot of them (otherwise keep all of them)
        if len(stv_df) > number_to_sample:
            stv_df = stv_df.sample(number_to_sample, replace=True)
        # drop column with cenhap assignment
        stv_df = stv_df.drop(['cenhap'], axis=1)
        # write to file in output directory
        stv_df.to_csv('bed/ch{}.{}.bed'.format(cenhap, frequent_stv.split('.')[1]), sep='\t', index=False, header=None)

