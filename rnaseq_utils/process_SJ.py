#!/usr/bin/env python

'''
Executable to process SJ.tab.out tables from STAR alignment and obtain a table with
counts for splice junction counts of skipped exon events, and of constitutive introns.

Input should be a directory of the format:
MAIN_DIRECTORY/<sample>.SJ.out.tab

Output is a directory with the following files:
    SE_counts.tab            (Table with skipped exon splice junction counts; should be relatively easy to adapt for other events)
    constitutive_introns.tab (Table with constitutive intron counts)
    
Additional input requirements:

    --skipped_exons skipped exons bed file with the following format:
        chr1    4293013 4311269 Rp1_1_I1        Rp1_1   -
        chr1    4293013 4351909 Rp1_1_SE        Rp1_1   -
        chr1    4311434 4351909 Rp1_1_I2        Rp1_1   -
        ...
        
        Each skipped exon is represented by three consecutive rows. I1 and I2 are the flanking introns when the exon is included; 
        SE is the intron when the exon is excluded.
        
    --constitutive_introns constitutive introns bed file with the following format:
    
        chr1    3216969 3421701 Xkr4    Xkr4_1  -
        chr1    3421902 3670551 Xkr4    Xkr4_2  -
        chr1    4352082 4352201 Rp1     Rp1_1   -
        ...
        
        Each constitutive intron corresponds to one row. These should be introns that are always removed in mature mRNA molecules of 
        a given gene (e.g., introns that appear in all isoforms of a given gene in a genome annotation file).
        
Latest update:
Carlos F Buen Abad Najar, June 25, 2020
cfbuenabadn@berkeley.edu
'''

import numpy as np
import pandas as pd
import subprocess as sp
import os
from collections import Counter
from tqdm import tqdm

import argparse

parser = argparse.ArgumentParser(description='Process SJ tables.')
parser.add_argument('-sj', '--sj_dir', type=str, nargs=1, required=True,
                   help='directory with <cell>.SJ.out.tab files')

parser.add_argument('-se', '--skipped_exons', type=str, nargs=1, required=True,
                   help='skipped exons bed file')

parser.add_argument('-ci', '--constitutive_introns', type=str, nargs=1, required=True,
                   help='constitutive introns bed file')

parser.add_argument('-o', '--out_dir', type=str, nargs=1, required=True,
                   help='output directory')

if __name__ == '__main__':
    
    # Parsing arguments
    
    args = parser.parse_args()
    SJ_dir = args.sj_dir[0]
    
    if not os.path.isdir(SJ_dir):
        os.mkdir(SJ_dir)
        
    out_dir = args.out_dir[0]
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    
    # Get cells and SJ files
    cells = [x.split('.')[0] for x in os.listdir(SJ_dir)]
    
    SJ_counts = {}
    SJ_converter = {}

    # Create directory with splice junctions in skipped exons
    with open(args.skipped_exons[0], 'r') as ase_file:
        for row in ase_file:
            sj = row.rstrip().split('\t')
            key = ':'.join([sj[0], sj[1], sj[2], sj[5]])
            if key in SJ_converter.keys():
                SJ_converter[key].append(sj[3])
            else:
                SJ_converter.update({key:[sj[3]]})
            SJ_counts.update({sj[3]:[0]*len(cells)})
            
            
    # Parse through each cell, obtain alternative splice junctions
    print('Creating skipped exons table')
    counter = 0
    for cell in tqdm(cells):
        SJ_table = SJ_dir + cell + '.SJ.out.tab'
        with open(SJ_table, 'r') as SJ_file:
            for row in SJ_file:
                sj = row.rstrip().split('\t')


                if sj[3] == '1':
                    strand = '+'
                elif sj[3] == '2':
                    strand = '-'
                else:
                    strand = '0'

                row_key = sj[0] + ':' + sj[1] + ':' + sj[2] + ':' + strand

                if row_key in SJ_converter.keys():
                    ase_names = SJ_converter[row_key]
                    for ase in ase_names:
                        SJ_counts[ase][counter] = int(sj[6]) #+ int(sj[7])

                row_key = sj[0] + ':' + sj[2] + ':' + sj[1]
                if row_key in SJ_converter.keys():
                    ase_names = SJ_converter[row_key]
                    for ase in ase_names:
                        SJ_counts[ase][counter] = int(sj[6]) #+ int(sj[7])
        counter += 1


    # Create table and output
    SJ_counts_table = pd.DataFrame.from_dict(SJ_counts)
    SJ_counts_table.index = cells
    SJ_counts_table.T.to_csv(out_dir + '/SE_counts.tab', sep='\t', index=True, header=True)


    # Create splice junction directories for constitutive introns
    SJ_counts = {}
    SJ_converter = {}

    with open(args.constitutive_introns[0], 'r') as ase_file:
        for row in ase_file:
            sj = row.rstrip().split('\t')
            key = ':'.join([sj[0], sj[1], sj[2], sj[5]])
            if key in SJ_converter.keys():
                SJ_converter[key].append(sj[4])
            else:
                SJ_converter.update({key:[sj[4]]})
            SJ_counts.update({sj[4]:[0]*len(cells)})


    # Parse through each cell, obtain constitutive splice junctions
    print('Creating constitutive introns table')
    counter = 0
    for cell in tqdm(cells):
        SJ_table = SJ_dir + cell + '.SJ.out.tab'
        with open(SJ_table, 'r') as SJ_file:
            for row in SJ_file:
                sj = row.rstrip().split('\t')


                if sj[3] == '1':
                    strand = '+'
                elif sj[3] == '2':
                    strand = '-'
                else:
                    strand = '0'

                row_key = sj[0] + ':' + sj[1] + ':' + sj[2] + ':' + strand

                if row_key in SJ_converter.keys():
                    ase_names = SJ_converter[row_key]
                    for ase in ase_names:
                        SJ_counts[ase][counter] = int(sj[6]) 

                row_key = sj[0] + ':' + sj[2] + ':' + sj[1]
                if row_key in SJ_converter.keys():
                    ase_names = SJ_converter[row_key]
                    for ase in ase_names:
                        SJ_counts[ase][counter] = int(sj[6]) 
        counter += 1
                        
    SJ_counts_table = pd.DataFrame.from_dict(SJ_counts)
    SJ_counts_table.index = cells
    
    SJ_counts_table.T.to_csv(out_dir + '/constitutive_introns.tab', sep='\t', index=True, header=True)
        