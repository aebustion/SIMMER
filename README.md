# SIMMER
Use chemical and protein similarity to identify microbiome-mediated enzymatic reactions

## Conda environment
To run SIMMER, the user will need the appropriate conda environment.

`$ conda env create --file SIMMER.yml`

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

