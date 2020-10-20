# Halide sites
## Project availability
You can clone whole repo by standart command:
```
git clone https://github.com/rostkick/Halide_sites.git
```
Also if you need to save it **without commit-log** you schould use:
```
git clone —depth=1 https://github.com/rostkick/Halide_sites.git
```
If you want to save **only one** file (for example without saving example data) you should:
  1) copy link of favorite file
  2) go to https://minhaskamal.github.io/DownGit/#/home, paste it, and create Download Link.
## Dependencies
You should install it manually: 
* [Anaconda](https://www.digitalocean.com/community/tutorials/how-to-install-anaconda-on-ubuntu-18-04-quickstart)
* [FreeSASA](https://freesasa.github.io/)  
## Preparation  
Create your conda environment from our smk.yml  
```
conda env create -f smk.yml --name smk
```
Activate your environment (some of this should work)  
```
conda activate smk
```
or
```
source activate smk
```

Make sure everything works on `--dryrun`.
```
snakemake all -n
```
## How to use
```
snakemake all
```
## DAG of jobs  
![alt text](dag.svg)  
## Pipeline input  
Pipeline automatically creates input consisting of all PDB entries with halide sites 

Or you can collect your input by following rule:

Input consist of PDB ID, which separeted by ```\n```:
```
3VRS
2LH5
6RI4
5EVC
...
```
You should put the files in the directory ```data/pdb_ID/```, and give them a name according to a strict template ```pdb_entries_X.txt```, where ```X``` is (```F, Cl, Br, I```)
## Scripts description  
1. rule PDB_par   ser(**PDB_parser.py**) - obtain pdb ID from PDB database
2. rule fetch_structures (**fetch_structures.py**) - obtain structures from PDB;  
3. rule filter (**filter.py**) - applies filters to PDB structures;  
4. rule get_context_data (**get_context.py**) - applies filters to halide binding sites;   
5. rule combine_final_data (**parse_context.py**) - aggregates data from all halides, then return 2 output files:
  * aggregated information about (model_name, distances, angels, fASA, atom_type, residue_type, water distribution, etc);
  * aggregated information about compositions of amino acids.  
6. rule figures_plotting (**figs.R**) - makes graphic report.  
## References  
1. Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web. <https://anaconda.com>.  
2. Mitternacht, S. FreeSASA: An open source C library for solvent accessible surface area calculations. F1000Research 5, 189 (2016).  
