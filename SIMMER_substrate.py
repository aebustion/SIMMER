# SIMMER.py
# Annamarie Bustion
# 2022_02
# 
#example query
#$ python3 SIMMER.py -i <path_to>/SIMMER_files -o <path_to_output_dir>
######################################

#import mkl
import pickle as pkl
import pandas as pd
import numpy as np
from numpy import linalg as LA
import scipy.stats as stats
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.AtomPairs import Pairs
import sys
import time
import argparse
import os
import shutil
import glob as glob


def run_substrate(row, df, output_dir):
#    print(df.iloc[row,3])
    sms_l = df.iloc[row,3].split('.')[0]
#    print(sms_l)
    m_l = Chem.MolFromSmiles(sms_l, sanitize=True)
    
    return m_l

def find_closest_rxns(query_id, X_mol, id_to_index):
    dists = {}
    vec = id_to_index[query_id]
    for i in range(len(X_mol)):
        dist = np.linalg.norm(X_mol[vec] - X_mol[i])
        dists[i] = dist
    s_dists = sorted(dists.items(),key=lambda x:x[1])
    results = []
    for ind, dist in s_dists:
        results.append([list(id_to_index.keys())[list(id_to_index.values()).index(ind)], str(dist)]) 
    [results.remove(result) for result in results if query_id in result]
    
    return [x for x in results if not x[0].startswith('DM')]
        

def fp_queries(dms_df, fps, output_dir):
    print("\nDescribing input reactions...")
    t0 = time.time()
    for i in range(len(dms_df)):
        try:
            sub = run_substrate(i, dms_df, output_dir)
            fp = Pairs.GetAtomPairFingerprintAsBitVect(sub)
            fps.append(fp)
        except:
            print("Improperly formatted SMILES; please refer to your compound's pubchem entry for appropriate SMILES")
            sys.exit()
    t1 = time.time()
    print("\nFinished creating fingerprints of queries in " + "{:.2f}".format(t1-t0) + " seconds")
    return fps


def add_queries_to_tanimoto(fps, tan_ar):
    t0 = time.time()
    metacyc_fps = fps[:8914]
    dms_fps = fps[8914:]
    tot_fps = metacyc_fps+dms_fps    
    tan_ar_dm = []
    for i in range(len(dms_fps)):
        tan = DataStructs.cDataStructs.BulkTanimotoSimilarity(dms_fps[i], tot_fps)
        tan_ar_dm.append(tan)
    tot_array = np.append(np.array(tan_ar),np.array(tan_ar_dm)[:,:8914], axis=0)
    X_mol = np.append(tot_array, np.array(tan_ar_dm).T, axis=1)  
    t1 = time.time()
    print("\nFinished computing and adding queries to tanimoto similarity matrix in " 
          + "{:.2f}".format(t1-t0) + " seconds")
    return X_mol


def take_ES_walk(ec_cat, ec_df, level):
    running_tally = []
    i = 0
    for ec in ec_df[level]:
        if ec==ec_cat:
            i = i+1
            running_tally.append(i)
        elif ec =='NIL':
            i = i+0
            running_tally.append(i)
        elif ec == 'DM':
            i = i+0
            running_tally.append(i)
        else:
            i = i-1
            running_tally.append(i)
    return running_tally


def find_EC_pval(ec_cats, ec_df, level, perm):
    sig_es_score_results=[]
    for ec_cat in ec_cats:
        if ec_cat=='NIL':
            continue
        elif ec_cat=='DM':
            continue
        else: 
            es_tally = take_ES_walk(ec_cat, ec_df, level)
            es_score = max(es_tally)
            where = es_tally.index(es_score)
            es_score_norm = es_score/(ec_df[level].value_counts()[ec_cat])
            perm.append(es_score_norm)
            perm.sort(reverse=True)
            p_val = (perm.index(es_score_norm)+1)/len(perm)
            perm.remove(es_score_norm)
            if p_val<0.05:
                if es_score_norm>0:
                    sig_es_score_results.append([ec_cat, es_score_norm, p_val, where])
                
    return pd.DataFrame(sig_es_score_results, columns=['EC', 'enrich_score', 'p_val', 'where'])


def predict_all_ECs(ranked_list, id_to_index, perm_df, rxn_to_ec):
    ec_array = []
    for result in ranked_list:
        ec_array.append([rxn_to_ec[(result[0])], result[1]])
    ec_df = pd.DataFrame(ec_array)
    ec_df.columns = ['EC4', 'euc']
    for ec in range(1,4):
        ec_df['EC' + str(ec)]=['.'.join(e.split('.')[0:ec]) for e in ec_df['EC4']]
    ec_df = ec_df[['euc', 'EC1', 'EC2', 'EC3', 'EC4']]
    
    # first level
    sig_results1 = find_EC_pval(ec_df['EC1'].unique(), ec_df, 'EC1', perm_df['EC1'].to_list()).sort_values(by='where')
    message = "a signficant EC class prediction"
    
    # second level
    if sig_results1.shape[0]>0:
        ec_sub_df = ec_df.loc[ec_df['EC1'].isin(sig_results1['EC'])]
        sig_results2 = find_EC_pval(ec_sub_df['EC2'].unique(), ec_df, 'EC2', perm_df['EC2'].to_list()).sort_values(by='where')
        sig_results = pd.concat([sig_results1, sig_results2])
        
        # third level
        if sig_results2.shape[0]>0:
            ec_sub_df = ec_df.loc[ec_df['EC2'].isin(sig_results2['EC'])]
            sig_results3 = find_EC_pval(ec_sub_df['EC3'].unique(), ec_df, 'EC3', perm_df['EC3'].to_list()).sort_values(by='where')
            sig_results = pd.concat([sig_results, sig_results3])
            
        else:
            message = message + ", but no significant sub-class"        
    else:
        message = "no significant EC prediction"
        sig_results=sig_results1
    
    return sig_results.reset_index(drop=True), message


def pull_tsv_results(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t', header=None)
    
    df.columns= ['hit_id', 'hit_descript', 'Lineage', 'Phylum',
             'full_seq_evalue','full_seq_bitscore', 'total_len',
             'hit_num_domains', 'hit_domain_indices', 'best_domain_evalue',
             'Genome', 'Genome_type', 'Original_name', 'GC_content', 'genome','uhgg_name','genome2','?','prev',
                  'abund']
    
    return df


def return_query_results(DM, X_mol, id_to_index, rxn_to_ec, final, input_dir, output_dir,
                         perm_df, prot_dict):
    t0=time.time()
    ranked_list = find_closest_rxns(DM, X_mol, id_to_index)
    
    ecdf, message = predict_all_ECs(ranked_list, id_to_index, perm_df, rxn_to_ec)
    t1=time.time()
    print("\nFor " + DM + ":\nFinished predicting EC codes in " + "{:.2f}".format(t1-t0) + ' seconds,')
    print("and there was " + message)

    #find closest rxn with homologs
    for i in range(len(ranked_list)):
        read_files=[]
        num_lines=0
        closest_rxn, euc_dist = ranked_list[i]
        files = prot_dict[closest_rxn.split('_')[0]]
        for file in files:
            f = glob.glob(input_dir + "/prot_data/" + file + '*tsv')
            num_lines = num_lines + sum(1 for line in open(f[0]))
            read_files.append(f[0])
        if num_lines > 1:
            break

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
                  'abund']].to_csv(output_dir + '/' + DM + '_enzyme_predictions.tsv', sep='\t', index=None)
    except:
        output_df.to_csv(output_dir + '/' + DM + '_enzyme_predictions.tsv', sep='\t', index=None)
    
    #fastas
    read_files=[]
    for file in files:
        f = glob.glob(input_dir + "/prot_data/" + file + '*fasta')
        read_files.append(f[0])
    with open(output_dir + '/' + DM + '_enzyme_predictions.fasta','wb') as wfd:
        for fasta in read_files:
            with open(fasta,'rb') as fd:
                shutil.copyfileobj(fd, wfd)

    #EC predictions
    ecdf.to_csv(output_dir + '/' + DM + '_EC_predictions.tsv', sep='\t', index=None)
    
    #Ranked reactions
    pd.DataFrame(ranked_list, 
                 columns=['MetaCyc Rxn', 
                          'Euclidean distance to query']).to_csv(output_dir + '/' + DM + "_distance_ranked_reactions.tsv",
                                                                 sep='\t', index=None)
    
    print("All output files now in: " + output_dir + "\nand the closest MetaCyc reaction with gut homologs is " 
          + closest_rxn.split('_')[0])

    
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_dir',
                    help='input directory')
    parser.add_argument('-o', action='store', dest='output_dir',
                    help='output directory')
    parser.add_argument('-q', action='store', dest='query',
                    help='query tsv file location. leave empty to input directly')
    #parser.add_argument('-t', action='store', dest='num_threads',
    #                help='number of threads to use')
    #parser.add_argument('-f', action='store', dest='fp_type',
    #                help='store type of fingerprint: type options here')
    #parser.add_argument('-p', action='store', dest='pca',
    #                help='whether or not to use dimension reduction, yes/no')
    return parser.parse_args()


def main():
    arguments = get_arguments()
    #mkl.set_num_threads(arguments.num_threads)
    input_dir = arguments.input_dir
    output_dir = arguments.output_dir
    query = arguments.query
    #fp_type = arguments.fp_type
    #pca = arguments.pca

    #load MetaCyc data
    #load data
    fps = pkl.load(open(input_dir + "/chem_data/MC_sub_fps.p", 'rb'))
    tan_ar = pkl.load(open(input_dir + "/chem_data/MC_sub_tanimoto_matrix.p", 'rb'))
    id_to_index = pkl.load(open(input_dir + "/chem_data/MC_rxn_to_matrix_index_dict.p", 'rb'))
    rxn_to_ec = pkl.load(open(input_dir + "/chem_data/MC_rxn_ec_dict.p", 'rb'))
    final = pd.read_pickle(input_dir + "/chem_data/MetaCyc_reactions_tsv.p")
    perm_df = pd.read_csv(input_dir + "/chem_data/ec_permutations.csv")
    prot_dict = pkl.load(open(input_dir + "/prot_data/prot_dict.p", 'rb'))
    print('\nPrecomputed data loaded')
    
    #load and run query
    if query != None:
        dms_df = pd.read_csv(query, sep='\t')
        if dms_df.shape[1] != 5:
            dms_df = pd.read_csv(query, sep=',')
            if dms_df.shape[1] != 5:
                print("Improperly formatted input tsv. Make sure file is tab delimited, and contains 5 columns of information (DM<id>, substrate name, product name, substrate(s) SMILES delimited by . , product(s) SMILES delimited by . )")
                sys.exit()
    else:
        reaction = input("\nEnter a substrate identifier (e.g. DM1): ")
        left_comp = input("\nEnter substrate name (e.g. gemcitabine): ")
        right_comp = 'placeholder' #input("\nEnter product name (e.g. 2′, 2′-difluorodeoxyuridine): ")
        left_smiles = input("\nEnter the SMILES of substrate \nexample OC[C@H]1O[C@@H](N2C=CC(=O)NC2=O)C(F)(F)[C@@H]1O : ") 
        right_smiles = 'placeholder' #input("\nEnter the SMILES of product(s) separated by '.' \nexample OC[C@H]1O[C@@H](N2C=CC(=O)NC2=O)C(F)(F)[C@@H]1O.[NH4+] : ") 

        if len(left_smiles) + len(right_smiles) != 0:
        
            dms_df = pd.DataFrame([[reaction, 
                                left_comp, 
                                right_comp, 
                                left_smiles,
                                right_smiles]], 
                              columns = ['reaction',
                                         'left_comp', 'right_comp', 
                                         'left_smiles','right_smiles'])

        elif len(left_smiles) + len(right_smiles) ==0:
            print('No queries provided')
            sys.exit()
        
    X_mol = add_queries_to_tanimoto(fp_queries(dms_df, fps, output_dir), tan_ar)
    
    ### add in eigendecomposition stuff later (will require the eigen fx, num pcs, and num_threads ###
    
    for i in range(len(dms_df)):
        id_to_index[dms_df.iloc[i,0]] = i+8914
    for i in range(len(dms_df)):
        rxn_to_ec[dms_df.iloc[i,0]] = 'DM'
    
    for DM in dms_df['reaction'].values:
        return_query_results(DM, X_mol, id_to_index, rxn_to_ec, final, input_dir,
                             output_dir, perm_df, prot_dict)


if __name__ == '__main__':
    main()
