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

parser.add_argument('-el', '--effective_length',  action='store_true',
                   help='collect effective length tables')

parser.add_argument('-gz', '--gzipped_files', action='store_true',
                   help='input files are zipped')

parser.add_argument('-se', '--skipped_exons', type=str, nargs=1, required=False,
                   help='skipped exons bed file')

parser.add_argument('-ci', '--constitutive_introns', type=str, nargs=1, required=False,
                   help='constitutive introns bed file')

parser.add_argument('-sp', '--split_cells', type=int, nargs=1, required=False, default = 1,
                   help='constitutive introns bed file')

parser.add_argument('-o', '--out_dir', type=str, nargs=1, required=True,
                   help='output directory')

def split_samples(cells, sp):
    step = int(len(cells)/sp)
    for i in range(0, len(cells), step):
        yield cells[i:i + step]
        
        
if __name__ == '__main__':
    
    args = parser.parse_args()
    alignment_dir = args.alignment_dir[0]
    out_dir = args.out_dir[0]
    strand = args.strand[0]
    effective_length = args.effective_length
    if args.gzipped_files:
        gz = '.gz'
    else:
        gz = ''
    
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
        
    
    samples_all = sorted(os.listdir(alignment_dir))    
    cells = [x.split('_')[0] for x in samples_all]
    
    cell_chunks = split_samples(samples_all, args.split_cells) 
    
    #print('Working on ' + str(args.split_cells) + ' chunks')
    
    counter_chunk = 1
    
    for samples in cell_chunks:
    
        if args.move_SJ:
            SJ_dir = out_dir + '/SJ_tables/'
            if not os.path.isdir(SJ_dir):
                os.mkdir(SJ_dir)

        else:
            skipped_exons = args.skipped_exons[0]
            constitutive_introns = args.constitutive_introns[0]

            print('Preparing SJ dictionaries')

            # Get cells and SJ files
            samples

            SJ_counts_se = {}
            SJ_converter_se = {}

            # Create directory with splice junctions in skipped exons
            with open(skipped_exons, 'r') as ase_file:
                for row in ase_file:
                    sj = row.rstrip().split('\t')
                    key = ':'.join([sj[0], sj[1], sj[2], sj[5]])
                    if key in SJ_converter_se.keys():
                        SJ_converter_se[key].append(sj[3])
                    else:
                        SJ_converter_se.update({key:[sj[3]]})
                    SJ_counts_se.update({sj[3]:[0]*len(samples)})

            # Create splice junction directories for constitutive introns
            SJ_counts_ci = {}
            SJ_converter_ci = {}

            with open(constitutive_introns, 'r') as ase_file:
                for row in ase_file:
                    sj = row.rstrip().split('\t')
                    key = ':'.join([sj[0], sj[1], sj[2], sj[5]])
                    if key in SJ_converter_ci.keys():
                        SJ_converter_ci[key].append(sj[4])
                    else:
                        SJ_converter_ci.update({key:[sj[4]]})
                    SJ_counts_ci.update({sj[4]:[0]*len(samples)})


#         table_el = pd.DataFrame()
        table_counts = pd.DataFrame()
        table_tpm = pd.DataFrame()
        ###
#         iso_table_el = pd.DataFrame()
        iso_table_counts = pd.DataFrame()
        iso_table_tpm = pd.DataFrame()
        iso_table_psi = pd.DataFrame()
        ###
        star_counts = pd.DataFrame()

        meta_table = pd.DataFrame()

        counter_se = 0
        counter_ci = 0



        for sample in tqdm(samples):
            tabla = pd.read_csv(os.path.join(alignment_dir,sample, 'rsem_output/' + args.rsem_name[0] + 'genes.results'+gz), 
                                sep = '\t', index_col=0)

            sample_name = sample.split('_')[0]

#             table_el[sample_name] = tabla['effective_length']
            table_counts[sample_name] = tabla['expected_count']
            table_tpm[sample_name] = tabla['TPM']


            ###
            iso_tabla = pd.read_csv(os.path.join(alignment_dir,sample, 'rsem_output/' + args.rsem_name[0] + 'isoforms.results'+gz), 
                                sep = '\t', index_col=0)

            sample_name = sample.split('_')[0]

#             iso_table_el[sample_name] = iso_tabla['effective_length']
            iso_table_counts[sample_name] = iso_tabla['expected_count']
            iso_table_tpm[sample_name] = iso_tabla['TPM']
            iso_table_psi[sample_name] = iso_tabla['IsoPct']
            ###

            star_tab = pd.read_csv(alignment_dir + sample + '/star_output/ReadsPerGene.out.tab'+gz, 
                               names = ['no_specific', 'first_strand', 'reverse_strand'],
                               sep='\t', index_col=0, skiprows=4)

            star_counts[sample_name] = star_tab[strand]


            log_df = pd.read_csv(alignment_dir + sample + '/star_output/Log.final.out'+gz, sep='\t', names=['names', 'value'])

            meta_list = [log_df.loc[4, 'value'],
                 log_df.loc[5, 'value'],
                 log_df.loc[7, 'value'],
                 log_df.loc[8, 'value'][:-1],
                 log_df.loc[9, 'value'],
                 log_df.loc[10, 'value'],
                 log_df.loc[11, 'value'],
                 log_df.loc[12, 'value'],
                 log_df.loc[13, 'value'],
                 log_df.loc[14, 'value'],
                 log_df.loc[15, 'value'],
                 log_df.loc[22, 'value'],
                 log_df.loc[23, 'value'][:-1],
                 log_df.loc[24, 'value'],
                 log_df.loc[25, 'value'][:-1]]

            meta_table[sample_name] = meta_list

            if args.move_SJ:
                SJ_tab = alignment_dir + sample + '/star_output/SJ.out.tab'+gz
                SJ_tab_out = SJ_dir + sample_name + '.SJ.out.tab'
                sp.run('cp {SJ} {out}'.format(SJ=SJ_tab, out=SJ_tab_out), shell=True)

            else:

                SJ_table = alignment_dir + sample + '/star_output/SJ.out.tab'+gz
                with open(SJ_table, 'r') as SJ_file:
                    for row in SJ_file:
                        sj = row.rstrip().split('\t')


                        if sj[3] == '1':
                            strand_se = '+'
                        elif sj[3] == '2':
                            strand_se = '-'
                        else:
                            strand_se = '0'

                        row_key = sj[0] + ':' + sj[1] + ':' + sj[2] + ':' + strand_se

                        if row_key in SJ_converter_se.keys():
                            ase_names = SJ_converter_se[row_key]
                            for ase in ase_names:
                                SJ_counts_se[ase][counter_se] = int(sj[6]) #+ int(sj[7])

                        row_key = sj[0] + ':' + sj[2] + ':' + sj[1]
                        if row_key in SJ_converter_se.keys():
                            ase_names = SJ_converter_se[row_key]
                            for ase in ase_names:
                                SJ_counts_se[ase][counter_se] = int(sj[6]) #+ int(sj[7])
                counter_se += 1


                with open(SJ_table, 'r') as SJ_file:
                    for row in SJ_file:
                        sj = row.rstrip().split('\t')


                        if sj[3] == '1':
                            strand_ci = '+'
                        elif sj[3] == '2':
                            strand_ci = '-'
                        else:
                            strand_ci = '0'

                        row_key = sj[0] + ':' + sj[1] + ':' + sj[2] + ':' + strand_ci

                        if row_key in SJ_converter_ci.keys():
                            ase_names = SJ_converter_ci[row_key]
                            for ase in ase_names:
                                SJ_counts_ci[ase][counter_ci] = int(sj[6]) 

                        row_key = sj[0] + ':' + sj[2] + ':' + sj[1]
                        if row_key in SJ_converter_ci.keys():
                            ase_names = SJ_converter_ci[row_key]
                            for ase in ase_names:
                                SJ_counts_ci[ase][counter_ci] = int(sj[6]) 
                counter_ci += 1


        if not args.move_SJ:


            SJ_counts_table_se = pd.DataFrame.from_dict(SJ_counts_se)
            SJ_counts_table_se.index = [x.split('_')[0] for x in samples]
            SJ_counts_table_se.T.to_csv(out_dir + '/SE_counts_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
                                        sep='\t', index=True, header=True, compression='gzip')


            SJ_counts_table_ci = pd.DataFrame.from_dict(SJ_counts_ci)
            SJ_counts_table_ci.index = [x.split('_')[0] for x in samples]

            SJ_counts_table_ci.T.to_csv(out_dir + '/constitutive_introns_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
                                        sep='\t', index=True, header=True, compression='gzip')



        meta_table.index = meta_idx = ['input_reads', 'read_len', 'uniquely_mapped', 'unique_percent', 'map_len', 'splice_junctions',
                'annotated_SJ', 'GT/AG', 'GC/AG', 'AT/AC', 'non_canonical', 'multimap_reads', 'multimap_percent',
                'unmapped_reads', 'unmapped_percent']

        meta_table.to_csv(out_dir + '/star_meta_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
                          sep='\t', header=True, index=True, compression='gzip')

#         table_el.index = [x.split('.')[0] for x in table_el.index]
        table_counts.index = [x.split('.')[0] for x in table_counts.index]
        table_tpm.index = [x.split('.')[0] for x in table_tpm.index]
        star_counts.index = [x.split('.')[0] for x in star_counts.index]

        mart = pd.read_csv(args.gene_symbols[0], sep='\t', index_col=0)
        mart_clean = mart.drop_duplicates()
        mart_clean = mart_clean.groupby(mart_clean.index).first()

        good_genes = [x for x in table_tpm.index if x in mart_clean.index]

        symbol_name = mart_clean.columns[0]

#         table_filtered_el = table_el.loc[good_genes]
#         table_filtered_el.index = mart_clean.loc[good_genes, symbol_name]
#         table_filtered_el.columns = [x.split('_')[0] for x in table_filtered_el.columns]
#         table_filtered_el.to_csv(out_dir + '/effective_length_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
#                                  sep='\t', header=True, index=True, compression='gzip')

        table_filtered_counts = table_counts.loc[good_genes]
        table_filtered_counts.index = mart_clean.loc[good_genes, symbol_name]
        table_filtered_counts.columns = [x.split('_')[0] for x in table_filtered_counts.columns]
        table_filtered_counts.to_csv(out_dir + '/rsem_gene_counts_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
                                     sep='\t', header=True, index=True, compression='gzip')

        table_filtered_tpm = table_tpm.loc[good_genes]
        table_filtered_tpm.index = mart_clean.loc[good_genes, symbol_name]
        table_filtered_tpm.columns = [x.split('_')[0] for x in table_filtered_tpm.columns]
        table_filtered_tpm.to_csv(out_dir + '/rsem_gene_tpm_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
                                  sep='\t', header=True, index=True, compression='gzip')

        star_filtered = star_counts.loc[good_genes]
        star_filtered.index = mart_clean.loc[good_genes, symbol_name]
        star_filtered.columns = [x.split('_')[0] for x in star_filtered.columns]
        star_filtered.to_csv(out_dir + '/star_counts_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
                             sep='\t', header=True, index=True, compression='gzip')

#         iso_table_el.columns = [x.split('_')[0] for x in iso_table_el.columns]
#         iso_table_el.to_csv(out_dir + '/isoform_effective_length_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
#                             sep='\t', header=True, index=True, compression='gzip')

        iso_table_counts.columns = [x.split('_')[0] for x in iso_table_counts.columns]
        iso_table_counts.to_csv(out_dir + '/rsem_isoform_counts_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
                                sep='\t', header=True, index=True, compression='gzip')

        iso_table_tpm.columns = [x.split('_')[0] for x in iso_table_tpm.columns]
        iso_table_tpm.to_csv(out_dir + '/rsem_isoform_tpm_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
                             sep='\t', header=True, index=True, compression='gzip')

        iso_table_psi.columns = [x.split('_')[0] for x in iso_table_psi.columns]
        iso_table_psi.to_csv(out_dir + '/rsem_isoform_psi_{counter_chunk}.tab.gz'.format(counter_chunk=str(counter_chunk)), 
                             sep='\t', header=True, index=True, compression='gzip')
        
        counter_chunk += 1

    if args.split_cells == 1:
        
        sp.run('mv ' + out_dir + '/rsem_isoform_psi_1.tab.gz ' + out_dir + '/rsem_isoform_psi.tab.gz', shell=True)
        
        sp.run('mv ' + out_dir + '/rsem_isoform_tpm_1.tab.gz ' + out_dir + '/rsem_isoform_tpm.tab.gz', shell=True)
        
        sp.run('mv ' + out_dir + '/rsem_isoform_counts_1.tab.gz ' + out_dir + '/rsem_isoform_counts.tab.gz', shell=True)
        
        #sp.run('mv ' + out_dir + '/rsem_isoform_effective_length_1.tab ' + out_dir + '/rsem_isoform_effective_length.tab', shell=True)
        
        sp.run('mv ' + out_dir + '/star_counts_1.tab.gz ' + out_dir + '/star_counts.tab.gz', shell=True)
        
        sp.run('mv ' + out_dir + '/rsem_gene_tpm_1.tab.gz ' + out_dir + '/rsem_gene_tpm.tab.gz', shell=True)
        
        sp.run('mv ' + out_dir + '/rsem_gene_counts_1.tab.gz ' + out_dir + '/rsem_gene_counts.tab.gz', shell=True)
        
        #sp.run('mv ' + out_dir + '/effective_length_1.tab ' + out_dir + '/effective_length.tab', shell=True)
        
        sp.run('mv ' + out_dir + '/star_meta_1.tab.gz ' + out_dir + '/star_meta.tab.gz', shell=True)
        
        sp.run('mv ' + out_dir + '/SE_counts_1.tab.gz ' + out_dir + '/SE_counts.tab.gz', shell=True)
        
        sp.run('mv ' + out_dir + '/constitutive_introns_1.tab.gz ' + out_dir + '/constitutive_introns.tab.gz', shell=True)
        
#         ######
        
#     else:
#         del SJ_counts_table_se
#         del SJ_counts_table_ci
#         del meta_table
#         del table_el
#         del table_counts
#         del table_tpm
#         del star_counts
#         del table_filtered_el
#         del table_filtered_counts
#         del table_filtered_tpm 
#         del star_filtered 
#         del iso_table_el
#         del iso_table_counts
#         del iso_table_tpm
#         del iso_table_psi
        
#         for counter_chunk in range(1, args.split_cells[0]+1):
            
#         sp.run('mv ' + out_dir + '/rsem_isoform_psi_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/rsem_isoform_psi.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/rsem_isoform_tpm_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/rsem_isoform_tpm.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/rsem_isoform_counts_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/rsem_isoform_counts.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/rsem_isoform_effective_length_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/rsem_isoform_effective_length.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/star_counts_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/star_counts.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/rsem_gene_tpm_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/rsem_gene_tpm.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/rsem_gene_counts_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/rsem_gene_counts.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/effective_length_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/effective_length.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/star_meta_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/star_meta.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/SE_counts_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/SE_counts.tab', shell=True)
        
#         sp.run('mv ' + out_dir + '/constitutive_introns_{counter_chunk}.tab'.format(counter_chunk=str(counter_chunk)) + ' ' + out_dir + '/constitutive_introns.tab', shell=True)
        
               
    
    
    
