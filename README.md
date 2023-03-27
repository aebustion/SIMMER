# SIMMER
Use chemical and protein Similarity to Identify MicrobioMe Enzymatic Reactions. https://doi.org/10.1101/2022.08.02.502504

This tool can be used directly on the web at https://simmer.pollard.gladstone.org/. To query more than ten input reactions at a time, or to consider more than one closest MetaCyc reaction, please read on for use of the command-line tool.

![alt text](https://github.com/aebustion/SIMMER/blob/main/Figure1.png?raw=true)

Raise an issue or email annamarie.bustion@gladstone.ucsf.edu if you encounter any problems. This code is still getting actively updated.

## Conda environment
To run command-line SIMMER, the user will need the appropriate conda environment.

`$ conda env create --file SIMMER_local.yml`\
`$ source activate SIMMER_local`

If "ResolvePackageNotFound" error is returned, you will manually create the SIMMER env and install the following packages:

`$ conda create --name SIMMER_local python=3`\
`$ source activate SIMMER_local`\
`$ conda install -c conda-forge rdkit`\
`$ conda install -c conda-forge scipy`

If for some reason an error such as ModuleNotFoundError: No module named 'module_name' occurs when running SIMMER, please find and install the relevant package from [anaconda](https://anaconda.org/).

Once ready to run SIMMER, use\
`$ source activate SIMMER_local`

## Precomputed Data (updated files on 2022_07_28 - please redownload to be up-to-date)
This can be downloaded from https://simmer.pollard.gladstone.org//SIMMER_files.tar.gz

`$ wget https://simmer.pollard.gladstone.org//SIMMER_files.tar.gz`\
Once unzipped, the relevant files are located in <your_filepath>/usr/share/nginx/projects/simmer/SIMMER-website/SIMMER_files

(This precomputed data takes up ~18G once unzipped. A website with built-in database to make this lower impact on the user will be available soon. Check back here for the link.)

These files contain all the precomputed data needed to run any SIMMER queries: 
* precomputed hmmsearch results (tsvs and fastas) of the UHGG database
* fingerprints and pairwise tanimoto similarity matrix for all MetaCyc fingerprints
* dictionaries needed to link chemical and protein data
 
 ## Input format
 If the user has more than a few queries, they will want to input their queries as a tsv file organized as follows:
 
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
To see all argument options:
`$ python3 SIMMER.py -h`

to run with multiple queries in tsv format:
`$ python3 SIMMER.py -i <path_to_precomputed_data> -o <path_to_output_director> -q <path_to_query_df>`

or to run with command line entry for a single query at a time (follow command line prompts to input the query):
`$ python3 SIMMER.py -i /pollard/data/projects/drug_metabolism/SIMMER_files -o ./`

Additionally, the user may wish to see which microbiome enzymes may be associated with a MetaCyc reaction of interest (i.e. Reverse SIMMER):
`$ python3 metacyc_microbiome_homologs.py -i <path_to>/SIMMER_files -o <path_to_output_dir> -q RXN-11469`
