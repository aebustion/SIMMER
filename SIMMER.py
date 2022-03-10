# SIMMER.py
# Annamarie Bustion
# 2022_02
# 
# Parameters:
#    example query
#$ python3 SIMMER.py -i /pollard/data/projects/drug_metabolism/SIMMER_files -o ./ 
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
import sys
import time
import argparse
#import seaborn as sns


def run_rxn(row, df):
    sms_l = df.iloc[row,3].split(' // ')
    sms_r = df.iloc[row,4].split(' // ')
    
    smas_l = []
    for sm in sms_l:
        m = Chem.MolFromSmiles(sm.replace('R', '*'))
        sma = Chem.MolToSmarts(m, isomericSmiles=True)
        smas_l.append(sma)
    smas_r = []
    for sm in sms_r:
        m = Chem.MolFromSmiles(sm.replace('R', '*'))
        sma = Chem.MolToSmarts(m, isomericSmiles=True)
        smas_r.append(sma)
        
    left = '.'.join(map(str,smas_l))
    right = '.'.join(map(str,smas_r))
    
    rxn = rdChemReactions.ReactionFromSmarts(left + '>>' + right)
    return rxn


def do_eigen_decomp(similarities,output_dir,n_factors):
    print("\nStarting eigen decomp...")
    D, V = LA.eigh(np.array(similarities))
    for f in n_factors:
        t0 = time.time()
        V_fin = V[:,-int(f):]
        D_fin = D[-int(f):]
        D_half = np.sqrt(D_fin)
        D_half = np.diag(D_half)
        X_mol = np.dot(V_fin,D_half)
        t1 = time.time()
        print("\nFinished eigen decomp for " + str(f) + " factors in " + "{:.2f}".format(t1-t0) + " seconds")
        pkl.dump(X_mol, open(str(output_dir) + str(f) + '_factors.p', "wb"))

        
#necessary for handling the NILS in ec_dict while building list comprehension
def catch(f):
    try:
        return f
    except IndexError:
        return 'NIL'

    
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
        

def calc_odds_ratio(query, X_mol, ec_level, id_to_index, rxn_to_ec):
    df_result = pd.DataFrame([catch([rxn[0],
                                     rxn[1], 
                                     '.'.join(rxn_to_ec[rxn[0]].split('.')[:ec_level])]) for rxn in find_closest_rxns(query,X_mol, id_to_index)])
    
    #df_result = df_result[~df_result[0].str.startswith('DM')]
    df_result = df_result.loc[df_result[0]!=query]
    df_result = df_result.loc[df_result[2]!='DM']
    df_result = df_result.iloc[:25,:]
    df_result = df_result.loc[df_result[2]!='NIL']
    
    best_result = []
    for ec in df_result[2].unique():
        x1 = df_result.iloc[:25,:].groupby(2).count().iloc[df_result.iloc[:25,:].groupby(2).count().index.get_loc(ec)][0]
        x2 = 25-x1
        tot_ec_rxns = df_result.groupby(2).count().iloc[df_result.groupby(2).count().index.get_loc(ec)][0]
        x3 = tot_ec_rxns - x1
        x4 = len(X_mol)-tot_ec_rxns-x2
        table = np.array([[(x1+1), (x2+1)],[(x3+1), (x4+1)]])
        #this pvalue is the same as doing a survival function (1-cdf of hypergeometric)
        oddsratio, pvalue = stats.fisher_exact(table)
        best_result.append([query, 
                            ec, 
                            "{:.2f}".format(float(df_result.loc[df_result[2]==ec][1].min())), 
                            "{:.2f}".format(oddsratio), 
                            float("{:.2E}".format(pvalue))])
    
    best_result_df = pd.DataFrame(best_result)
    best_result_df.columns = ['query', 'predicted_EC', 'min_euc', 'OR', 'p_value']
    
    return best_result_df.loc[best_result_df['p_value']<0.05].sort_values(by='min_euc')


def make_ec_plot(DM):
    sns.set_style("white")
    df1 = pd.DataFrame(find_closest_rxns(DM, X_mol))
    df1['ec'] = [rxn_to_ec[rxn].split('.')[0] for rxn in df1[0]]
    df1.columns = ['rxn', 'distance', 'ec']
    skip = ['DM', 'NIL', 'missing']
    plot_df = df1[~df1['ec'].isin(skip)]
    plot_df['distance'] = plot_df['distance'].astype('float')
    g = sns.FacetGrid(plot_df, 
                      hue="ec", 
                      height=8.27, 
                      aspect=11.7/8.27, 
                      palette='Set2')
    g = g.map(sns.distplot, "distance",  hist=True, rug=True, kde=False)
    plt.legend()
    plt.title((DM), fontsize=20)
    plt.tick_params(labelsize=12)
    plt.xlabel('Euclidean distance', fontsize=20)
    plt.ylabel('Count', fontsize=20)
    plt.savefig(output_dir + "/" + DM + "_ec_histogram.png")


def fp_queries(dms_df, fps):
    print("\nDescribing input reactions...")
    t0 = time.time()
    for i in range(len(dms_df)):
        rxn = run_rxn(i, dms_df)
        fp = Chem.rdChemReactions.CreateDifferenceFingerprintForReaction(rxn)
        fps.append(fp)
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


def return_query_results(DM, X_mol, id_to_index, rxn_to_ec, final, output_dir):
    t0=time.time()
    closest_rxn, euc_dist = find_closest_rxns(DM, X_mol, id_to_index)[0]
    closest_ec = rxn_to_ec[closest_rxn]

    ec1df = calc_odds_ratio(DM, X_mol, 1, id_to_index, rxn_to_ec)
    ec1 = ec1df.iloc[0,1]
    ec2df = calc_odds_ratio(DM, X_mol, 2, id_to_index, rxn_to_ec)
    ec2df = ec2df.loc[ec2df['predicted_EC'].str.startswith(ec1)]
    ec3df = calc_odds_ratio(DM, X_mol, 3, id_to_index, rxn_to_ec)
    ec3df = ec3df.loc[ec3df['predicted_EC'].str.startswith(ec1)]
    
    #img = run_rxn(final.loc[final['reaction']==closest_rxn].index[0],final)
    #img.save(output_dir + '/' + DM + '_top_hit_' + closest_rxn + closest_ec +'.png')

    result_df = pd.concat([ec1df,ec2df,ec3df], 
          keys=['EC_1places','EC_2places', 'EC_3places'],
         names=['Series name'])
    
    t1=time.time()
    print("\nFor " + DM + ":\nFinished predicting EC codes in " + "{:.2f}".format(t1-t0) + ' seconds')
    print("All output files now in: " + output_dir + "\nand the closest MetaCyc reaction is " 
          + closest_rxn.split('_')[0] 
          + "_EC" 
          + closest_ec + '\n')
    
    result_df.to_csv(output_dir + '/' + DM + '_EC_predictions.tsv', sep='\t')
    with open(output_dir + '/' + DM + "_closest_rxn.txt", "w") as file:
        file.write(str([DM, closest_rxn.split('_')[0], euc_dist]))

    
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_dir',
                    help='input directory')
    parser.add_argument('-o', action='store', dest='output_dir',
                    help='output directory')
    parser.add_argument('-q', action='store', dest='query',
                    help='query tsv file location. leave empty to input directly')
    parser.add_argument('-t', action='store', dest='num_threads',
                    help='number of threads to use')
    parser.add_argument('-f', action='store', dest='fp_type',
                    help='store type of fingerprint: type options here')
    parser.add_argument('-p', action='store', dest='pca',
                    help='whether or not to use dimension reduction, yes/no')
    return parser.parse_args()


def main():
    arguments = get_arguments()
    #mkl.set_num_threads(arguments.num_threads)
    input_dir = arguments.input_dir
    output_dir = arguments.output_dir
    query = arguments.query
    fp_type = arguments.fp_type
    pca = arguments.pca

    #load MetaCyc data
    #load data
    fps = pkl.load(open(input_dir + "/MC_rxn_fps.p", 'rb'))
    tan_ar = pkl.load(open(input_dir + "/MC_rxn_tanimoto_matrix.p", 'rb'))
    id_to_index = pkl.load(open(input_dir + "/MC_rxn_to_matrix_index_dict.p", 'rb'))
    rxn_to_ec = pkl.load(open(input_dir + "/MC_rxn_ec_dict.p", 'rb'))
    final = pd.read_pickle(input_dir + "/MetaCyc_reactions_tsv.p")
    
    print('\nPrecomputed data loaded')
    
    #load and run query
    if query != None:
        dms_df = pd.read_csv(query, sep='\t')
    else:
        reaction = input("\nEnter a reaction identifier (e.g. DM1): ")
        left_comp = input("\nEnter substrate name (e.g. 4-asa): ")
        right_comp = input("\nEnter product name (e.g. n-acetyl-4-asa): ")
        left_smiles = input("\nEnter the SMILES of substrate(s) and known cofactors separated by '//' \ne.g. C1=CC(=C(C=C1N)O)C(=O)O // CC(=O)SCCNC(=O)CCNC(=O)[C@@H](C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O : ") 
        right_smiles = input("\nEnter the SMILES of product(s) separated by '//' \ne.g. CC(=O)NC1=CC(=C(C=C1)C(=O)O)O // CC(C)(COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)[C@H](C(=O)NCCC(=O)NCCS)O : ") 
        
        dms_df = pd.DataFrame([[reaction, 
                                left_comp, 
                                right_comp, 
                                left_smiles,
                                right_smiles]], 
                              columns = ['reaction',
                                         'left_comp', 'right_comp', 
                                         'left_smiles','right_smiles'])
        
    #fps = fp_queries(dms_df)
    X_mol = add_queries_to_tanimoto(fp_queries(dms_df, fps), tan_ar)
    
    #add in eigendecomposition stuff later ####
    
    for i in range(len(dms_df)):
        id_to_index[dms_df.iloc[i,0]] = i+8914
    for i in range(len(dms_df)):
        rxn_to_ec[dms_df.iloc[i,0]] = 'DM'
    
    for DM in dms_df['reaction'].values:
        return_query_results(DM, X_mol, id_to_index, rxn_to_ec, final, output_dir)


if __name__ == '__main__':
    main()