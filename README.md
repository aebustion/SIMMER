# SIMMER
Use chemical and protein Similarity to Identify MicrobioMe Enzymatic Reactions.

Bustion AE, Nayak RR, Agrawal A, Turnbaugh PJ, Pollard KS. 2023. SIMMER employs similarity algorithms to accurately identify human gut microbiome species and enzymes capable of known chemical transformations. Elife 12. doi:10.7554/eLife.82401

This tool can be used directly on the web at https://simmer.pollard.gladstone.org/. To query more than ten input reactions at a time, or to consider more than one closest MetaCyc reaction, please read on for use of the command-line tool.

![alt text](https://github.com/aebustion/SIMMER/blob/main/Figure1.png?raw=true)

Raise an issue or email aebustion@gmail.com if you encounter any problems. This code is still getting actively updated.

## Precomputed Data
This can be downloaded from https://simmer.pollard.gladstone.org//SIMMER_files.tar.gz

`$ wget https://simmer.pollard.gladstone.org//SIMMER_files.tar.gz`\
Once unzipped, the relevant files are located in <your_filepath>/usr/share/nginx/projects/simmer/SIMMER-website/SIMMER_files

(This precomputed data takes up ~18G once unzipped. If space is a constraint, please use the webtool version of SIMMER available at https://simmer.pollard.gladstone.org/.

These files contain all the precomputed data needed to run any SIMMER queries: 
* precomputed hmmsearch results (tsvs and fastas) of the UHGG database
* fingerprints and pairwise tanimoto similarity matrix for all MetaCyc fingerprints
* dictionaries needed to link chemical and protein data

## Conda environment
To run command-line SIMMER, the user will need the appropriate conda environment.

`$ conda env create --file SIMMER_local.yml`\
`$ source activate SIMMER_local`

If "ResolvePackageNotFound" error is returned, you can manually create the SIMMER env and install the following packages:

`$ conda install -c conda-forge rdkit`\
`$ conda install -c conda-forge scipy`

If for some reason an error such as ModuleNotFoundError: No module named 'module_name' occurs when running SIMMER, please find and install the relevant package from [anaconda](https://anaconda.org/).

Once ready to run SIMMER, use\
`$ conda activate SIMMER_local`\
`$ python3 SIMMER.py -h`\
example commands are at the end of this README.
 
 ## Input format
 If the user has more than a few queries, they will want to input their queries as a tsv or csv file organized as seen in the multiple_queries_example.csv and as seen here:
 
| reaction | left_comp | right_comp | left_smiles | right_smiles |
| -------- | --------- | ---------- | ----------- | ------------ |
| DM1 | brivudine | bromovinyluracil | OC[C@H]1O[C@H](C[C@@H]1O)N1C=C(\C=C\Br)C(=O)NC1=O | Br\C=C\C1=CNC(=O)NC1=O |
| DM2 | gemcitabine | 2′, 2′-difluorodeoxyuridine | NC1=NC(=O)N(C=C1)[C@@H]1O[C@H](CO)[C@@H](O)C1(F)F.O.[H+]  | OC[C@H]1O[C@@H](N2C=CC(=O)NC2=O)C(F)(F)[C@@H]1O.[NH4+] |
| DM3 | foo | bar | C1=CC=C2C(=C1)C(=NO2)CS(=O)(=O)N  | O=C(CS(=O)(N)=O)C1=CC=CC=C1O |
 
 Notes:
 * The user can enter any text into the left_comp and right_comp columns (as seen in the DM3 example), but for the reaction column **must identify a query using 'DM' followed by any alphanumeric character** (SIMMER distinguishes queries from its database of MetaCyc reactions by 'seeing' the DM label).
 * If the query includes multiple substrates and/or products (as seen in the DM2 example), their respective SMILES **must be separated by a '.'**
 * SIMMER can process either canonical or isomeric SMILES; the MetaCyc database uses isomeric SMILES where available.
 
 
 ## Output
 After running a SIMMER query, the user will obtain five files in about ~10-13 seconds:
 * DM_distance_ranked_rxns.txt (a file containing all MetaCyc reactions' euclidean distances to the input query, sorted from closest to farthest)
 * DM_EC_predictions.tsv (a tsv file containing SIMMER's predicted EC numbers for the query reaction)
 * DM_enzyme_predictions.fasta (a fasta file containing all predicted bacterial enzymes in the human gut capable of the query reaction)
 * DM_enzyme_predictions.tsv (a  tsv file summarizing all predicted bacterial enzymes in the human gut capable of the query reaction)
 * coming soon: DM_enzyme_predictions.png (a tree representation of the above predicted bacterial enzymes in the human gut)
 
 ## Run an example query using command line SIMMER tool
To see all argument options:\
`$ python3 SIMMER.py -h`

to run for a single query at a time (follow command line prompts to input the query):\
`$ python3 SIMMER.py -i <path_to>/SIMMER_files -o <path_to_output_directory>`

or to run with multiple queries in tsv or csv format (try using multiple_queries_example.csv as an example):\
`$ python3 SIMMER.py -i <path_to>/SIMMER_files -o <path_to_output_dir> -q <path_to_query_df>`

Additionally, the user may wish to see which microbiome enzymes may be associated with a MetaCyc reaction of interest (i.e. Reverse SIMMER):\
`$ python3 metacyc_microbiome_homologs.py -i <path_to>/SIMMER_files -o <path_to_output_dir> -q RXN-11469`
