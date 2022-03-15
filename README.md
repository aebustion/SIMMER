# SIMMER
Use chemical and protein similarity to identify microbiome-mediated enzymatic reactions

## Conda environment
To run command-line SIMMER, the user will need the appropriate conda environment.

`$ conda env create --file SIMMER.yml`

`$ source activate SIMMER`

## Precomputed Data
This can be found at /pollard/data/projects/drug_metabolism/SIMMER_files or can be downloaded at https://www.dropbox.com/s/1u21oha2y0cxv2q/SIMMER_files.tar.gz?dl=0

These files contain all the precomputed data needed to run any SIMMER queries: 
* precomputed hmmsearch results (tsvs and pngs of phylogenies) of the UHGG database (this is just a very small sample of the full database at the moment 03/2022)
* fingerprints of all MetaCyc reactions, 
* a pairwise tanimoto similarity matrix for all MetaCyc fingerprints, and
* three dictionaries
  * MetaCyc reaction:EC commission numbers, 
  * MetaCyc reaction:PANTHER subfamily, and 
  * MetaCyc reaction:tanimoto matrix index

 
 ## Input format
 If the user has more than a few queries, they will want to input their queries as a tsv file organized as follows:
 
| reaction | left_comp | right_comp | left_smiles | right_smiles |
| -------- | --------- | ---------- | ----------- | ------------ |
| DM1 | brivudine | bromovinyluracil | OC[C@H]1O[C@H](C[C@@H]1O)N1C=C(\C=C\Br)C(=O)NC1=O | Br\C=C\C1=CNC(=O)NC1=O |
| DM2 | gemcitabine | 2′, 2′-difluorodeoxyuridine | NC1=NC(=O)N(C=C1)[C@@H]1O[C@H](CO)[C@@H](O)C1(F)F // O // [H+]  | OC[C@H]1O[C@@H](N2C=CC(=O)NC2=O)C(F)(F)[C@@H]1O // [NH4+] |
| DM3 | foo | bar | C1=CC=C2C(=C1)C(=NO2)CS(=O)(=O)N  | O=C(CS(=O)(N)=O)C1=CC=CC=C1O |
 
 Notes:
 * The user can enter any text into the left_comp and right_comp columns (as seen in the DM3 example), but for the reaction column **must identify a query using 'DM' followed by any alphanumeric character** (SIMMER distinguishes queries from its database of MetaCyc reactions by 'seeing' the DM label).
 * If the query includes multiple substrates and/or products (as seen in the DM2 example), their respective SMILES **must be separated by a ' // '**
 * SIMMER can process either canonical or isomeric SMILES; the MetaCyc database uses isomeric SMILES where available.
 
 
 ## Output
 After running a SIMMER query, the user will obtain five files:
 * DM_distance_ranked_rxns.txt (a file containing all MetaCyc reactions' euclidean distances to the input query, sorted from closest to farthest)
 * DM_EC_predictions.tsv (a tsv file containing SIMMER's predicted EC numbers for the query reaction)
 * DM_enzyme_predictions.fasta (a fasta file containing all predicted bacterial enzymes in the human gut capable of the query reaction)
 * DM_enzyme_predictions.tsv.pkl (a pickled tsv file summarizing all predicted bacterial enzymes in the human gut capable of the query reaction)
 * DM_enzyme_predictions.png (a tree representation of the above predicted bacterial enzymes in the human gut)
 
 ## Run an example query using command line SIMMER tool

to run with multiple queries in tsv format:

`$ python3 SIMMER.py -i <path_to_precomputed_data> -o <path_to_output_director> -q <path_to_query_df>`

`$ python3 SIMMER.py -i /pollard/data/projects/drug_metabolism/SIMMER_files -o ./ -q ./queries_example.tsv`

or to run with single command line entry (follow command line prompts to input the query):

`$ python3 SIMMER.py -i /pollard/data/projects/drug_metabolism/SIMMER_files -o ./`
