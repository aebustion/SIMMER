# metacyc_microbiome_homologs.py
# Annamarie Bustion
# 2022_10
# 
#example query
#$ python3 metacyc_microbiome_homologs.py -i <path_to>/SIMMER_files -o <path_to_output_dir> -q RXN-11469
######################################

import pickle as pkl
import pandas as pd
import numpy as np
import sys
import argparse
import os
import shutil
import glob as glob


def pull_tsv_results(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t', header=None)
    
    df.columns= ['hit_id', 'hit_descript', 'Lineage', 'Phylum',
             'full_seq_evalue','full_seq_bitscore', 'total_len',
             'hit_num_domains', 'hit_domain_indices', 'best_domain_evalue',
             'Genome', 'Genome_type', 'Original_name', 'GC_content', 'genome','uhgg_name','genome2','?','prev',
                  'abund']
    
    return df


def return_query_results(input_dir, output_dir, prot_dict, query):

    read_files=[]

    #find MetaCyc rxn homologs
    files = prot_dict[query]
    for file in files:
        f = glob.glob(input_dir + "/prot_data/" + file + '*tsv')
        read_files.append(f[0])
        
    #tsvs
    output_df=pd.DataFrame()
    for tsv in read_files:
        try:
            output_df = pd.concat([output_df, pull_tsv_results(tsv)])
        except:
            continue
    try:
        output_df[['hit_id', 'hit_descript', 'Lineage', 'Phylum',
             'full_seq_evalue','full_seq_bitscore', 'total_len',
             'hit_num_domains', 'hit_domain_indices', 'best_domain_evalue',
             'Genome', 'Genome_type', 'Original_name', 'GC_content', 'uhgg_name','prev',
                  'abund']].to_csv(output_dir + '/' + query + '_enzyme_predictions.tsv', sep='\t', index=None)
    except:
        output_df.to_csv(output_dir + '/' + query + '_enzyme_predictions.tsv', sep='\t', index=None)
    
    #fastas
    read_files=[]
    for file in files:
        f = glob.glob(input_dir + "/prot_data/" + file + '*fasta')
        read_files.append(f[0])
    with open(output_dir + '/' + query + '_enzyme_predictions.fasta','wb') as wfd:
        for fasta in read_files:
            with open(fasta,'rb') as fd:
                shutil.copyfileobj(fd, wfd)

    print("All output files now in: " + output_dir)

    
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_dir',
                    help='input directory')
    parser.add_argument('-o', action='store', dest='output_dir',
                    help='output directory')
    parser.add_argument('-q', action='store', dest='query',
                    help='MetaCyc rxn for which the user wants microbiome homologs')
    return parser.parse_args()


def main():
    arguments = get_arguments()
    input_dir = arguments.input_dir
    output_dir = arguments.output_dir
    query = arguments.query

    prot_dict = pkl.load(open(input_dir + "/prot_data/prot_dict.p", 'rb'))
    print('\nPrecomputed data loaded')
    
    #load and run query
    if query != None:
        return_query_results(input_dir, output_dir, prot_dict, query)
    else:
        print("No query provided")

if __name__ == '__main__':
    main()
