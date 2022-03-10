import os
import re
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import SearchIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle as pkl
import sys
import argparse


# functions
def pull_hms_results(input_dir, output_dir, dm, tax_df):
    try:
        df_existing = pd.read_pickle(output_dir + dm + "_UHGG_hms_tsv.pkl")
    except FileNotFoundError:
        df_existing = pd.DataFrame()
    
    dict = {}
    i = 0
    my_recs = []
    query_len=0
    
    for f in glob.glob(input_dir + "/MGYG*/*" + dm + "_*.hms"):
        if query_len == 0:
            with open(f, 'r') as fub:
                lines = fub.readlines()
                for line in lines:
                    if re.search(r'Query:', line):
                        query_len = int(line.split('[')[-1].split('=')[-1][:-2])

        for qresult in SearchIO.read(f, "hmmer3-text"):
            rec = qresult[0].hit
            best_domain_evalue = qresult[0].evalue
            seq_id = qresult.id
            full_seq_evalue = qresult.evalue  
            bitscore = qresult.bitscore
            bias = qresult.bias
            description = qresult.description
            hmm_from = qresult[0].query_start
            hmm_to = qresult[0].query_end
            aln_span = qresult[0].aln_span
            hmm_span = hmm_to - hmm_from
            my_recs.append(rec)
            
            dict[i] = {'hit_id': seq_id, 
                       'hit_descript': description, 
                       'full_seq_evalue': full_seq_evalue,
                       'full_seq_bitscore': bitscore,
                       'full_seq_bias': bias,
                       'best_domain_evalue': best_domain_evalue, 
                       'hmm_from': hmm_from, 
                       'hmm_to': hmm_to, 
                       'hmm_span': hmm_span,
                       'aln_span': aln_span}
            i = i + 1
        
    df = pd.DataFrame.from_dict(dict, "index")
    genome = ['_'.join(i.split('_')[0:2]) for i in df['hit_id']]
    df['Genome']=genome
    df = pd.merge(df, tax_df, on='Genome')
    df = df[['hit_id', 'hit_descript', 'Lineage','full_seq_evalue','full_seq_bitscore', 'aln_span', 'Genome', 'Genome_type', 'Original_name', 'GC_content']]
    
    df = df.loc[df['aln_span']>=query_len*.5]
    df_new = df_existing.append(df).drop_duplicates()
    df_new.to_pickle(output_dir + dm + "_UHGG_hms_tsv.pkl", protocol=pkl.HIGHEST_PROTOCOL)
    
    try:
        recs = list(SeqIO.parse(output_dir + dm + "_UHGG_hits.fasta", "fasta"))
    except FileNotFoundError:
        recs = []
    for seq_record in my_recs:
        if seq_record.id in df_new['hit_id'].values:
            seq_record.seq = seq_record.seq.ungap('-')
            recs.append(seq_record)
    SeqIO.write(recs, output_dir + dm + "_UHGG_hits.fasta", "fasta")


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_dir',
                    help='input directory')
    parser.add_argument('-o', action='store', dest='output_dir',
                    help='output directory')
    parser.add_argument('-q', action='store', dest='query',
                    help='query tsv file location. leave empty to input directly')
    return parser.parse_args()


def main():
    arguments = get_arguments()
    input_dir = arguments.input_dir
    output_dir = arguments.output_dir
    query = arguments.query

    #load data
    tax_df = pd.read_csv("/pollard/data/wynton/consortia/UHGG_v1/genomes-all_metadata.tsv", sep='\t')
    
    #load and run query
    pull_hms_results(input_dir, output_dir, query, tax_df=tax_df)


if __name__ == '__main__':
    main()