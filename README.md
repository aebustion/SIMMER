# SIMMER
Use chemical and protein similarity to identify microbiome-mediated enzymatic reactions

## Conda environment
To run command-line SIMMER, the user will need the appropriate conda environment.

`$ conda env create --file SIMMER.yml`

`$ source activate SIMMER`

## Precomputed Data
Download: https://www.dropbox.com/s/i9dhjnepg0zl555/SIMMER_files.tar.gz?dl=0


`$ wget https://www.dropbox.com/s/i9dhjnepg0zl555/SIMMER_files.tar.gz?dl=0`

These files contain all the precomputed data needed to run any SIMMER queries: 
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
 
 Notes:
 * SIMMER distinguishes queries from its database of MetaCyc reactions by 'seeing' a reaction labeled as DMX (where X can be any alphanumeric identifier) 
 * SIMMER can process either canonical or isomeric SMILES, and the MetaCyc database uses isomeric SMILES where possible.
 * If the query includes multiple substrates and/or products, these smiles must be separated by a ' // ' 
 
 ## Run an example query using command line SIMMER tool
 
`$ python3 SIMMER.py`
