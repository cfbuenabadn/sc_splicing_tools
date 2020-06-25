#!/usr/bin/env python

'''
Executable for preprocessing alignment outputs from RSEM and STAR. These are intended for 
single cell RNA-seq datasets, but I see no reasons why it wouldn't work for bulk RNA either.

The input for this script is a directory with the following format:
MAIN_DIRECTORY/<SAMPLE>/
    star_output/ (output directory of STAR alignment run)
    rsem_output/ (output directory of RSEM alignment and quantification run)
    
The output is a directory with the following files:
    effective_length.tab (effective length column from RSEM gene alignment)
    rsem_gene_counts.tab (expected count column from RSEM gene alignment)
    rsem_gene_tpm.tab    (TPM column from RSEM gene alignment)
    star_counts.tab      (read counts from STAR alignment)
    
This code should be easy to modify to include the isoform alignments of RSEM.
    
Additionally, it has the option to collect the SJ.out.tab tables from the STAR run and 
store them in a SJ_tables directory in a format compatible for process_SJ.py

Latest update:
Carlos F Buen Abad Najar, June 25, 2020
cfbuenabadn@berkeley.edu
'''

import numpy as np
import pandas as pd
import subprocess as sp   
import os
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='Process RNASeq alignments to get counts and TPM tables.')

parser.add_argument('-ad', '--alignment_dir', type=str, nargs=1, required=True,
                   help='directory with <sample>/<aligner>_output/rsem_output files')

parser.add_argument('-g', '--gene_symbols', type=str, nargs=1, required=True,
                   help='table file with encode-to-gene symbol translations')

parser.add_argument('-rn', '--rsem_name', type=str, nargs=1, required=True, default='',
                   help='prefix for rsem results tables. E.g.: "rsem_output."')

parser.add_argument('-s', '--strand', type=str, nargs=1, required=True, default='no_specific',
                   help='"no_specific", "first_strand" or "reverse_strand", depending on the dataset')

parser.add_argument('-sj', '--move_SJ',  action='store_true',
                   help='collect SJ.out.tab files and put them in a directory for process_SJ.py')

parser.add_argument('-o', '--out_dir', type=str, nargs=1, required=True,
                   help='output directory')

if __name__ == '__main__':
    
    args = parser.parse_args()
    alignment_dir = args.alignment_dir[0]
    out_dir = args.out_dir[0]
    strand = args.strand[0]
    
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    
    if args.move_SJ:
        SJ_dir = out_dir + '/SJ_tables/'
        if not os.path.isdir(SJ_dir):
            os.mkdir(SJ_dir)
    
    samples = sorted(os.listdir(alignment_dir))

    table_el = pd.DataFrame()
    table_counts = pd.DataFrame()
    table_tpm = pd.DataFrame()
    star_counts = pd.DataFrame()
    
    for sample in tqdm(samples):
        tabla = pd.read_csv(os.path.join(alignment_dir,sample, 'rsem_output/' + args.rsem_name[0] + 'genes.results'), 
                            sep = '\t', index_col=0)
        table_el[sample] = tabla['effective_length']
        table_counts[sample] = tabla['expected_count']
        table_tpm[sample] = tabla['TPM']
        
        star_tab = pd.read_csv(star_dir + sample + '/star_output/ReadsPerGene.out.tab', 
                           names = ['no_specific', 'first_strand', 'reverse_strand'],
                           sep='\t', index_col=0, skiprows=4)
        
        star_counts[sample] = star_tab[strand]
        
        if args.move_SJ:
            SJ_tab = star_dir + sample + '/star_output/SJ.out.tab'
            SJ_tab_out = SJ_dir + sample + '.SJ.out.tab'
            sp.run('cp {SJ} {out}'.format(SJ=SJ_tab, out=SJ_tab_out), shell=True)
    
        
    table_el.index = [x.split('.')[0] for x in table_el.index]
    table_counts.index = [x.split('.')[0] for x in table_counts.index]
    table_tpm.index = [x.split('.')[0] for x in table_tpm.index]
    star_counts.index = [x.split('.')[0] for x in star_counts_tab.index]
    
    mart = pd.read_csv(args.gene_symbols[0], sep='\t', index_col=0)
    mart_clean = mart.drop_duplicates()
    mart_clean = clean.groupby(clean.index).first()

    good_genes = [x for x in table_tpm.index if x in mart_clean.index]
    
    table_filtered_el = table_el.loc[good_genes]
    table_filtered_el.index = mart_clean.loc[good_genes].mgi_symbol
    table_filtered_el.columns = [x.split('_')[0] for x in table_filtered_el.columns]
    table_filtered_el.to_csv(out_dir + '/effective_length.tab', sep='\t', header=True, index=True)
    
    table_filtered_counts = table_counts.loc[good_genes]
    table_filtered_counts.index = mart_clean.loc[good_genes].mgi_symbol
    table_filtered_counts.columns = [x.split('_')[0] for x in table_filtered_counts.columns]
    table_filtered_counts.to_csv(out_dir + '/rsem_gene_counts.tab', sep='\t', header=True, index=True)
    
    table_filtered_tpm = table_tpm.loc[good_genes]
    table_filtered_tpm.index = mart_clean.loc[good_genes].mgi_symbol
    table_filtered_tpm.columns = [x.split('_')[0] for x in table_filtered_tpm.columns]
    table_filtered_tpm.to_csv(out_dir + '/rsem_gene_tpm.tab', sep='\t', header=True, index=True)
    
    star_filtered = star_counts.loc[good_genes]
    star_filtered.index = mart_clean.loc[good_genes].mgi_symbol
    star_filtered.columns = [x.split('_')[0] for x in star_filtered.columns]
    star_filtered.to_csv(out_dir + '/star_counts.tab', sep='\t', header=True, index=True)
    
    

    
    
    
    
    
    