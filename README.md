# TDA-Metaganomas
In this repository we have some tools that we will use to distinguish the presence or absence of species in a metagenome made up of closely related species. In particular, we want to determine the presence of *Clavibacter michiganensis michiganensis* in a sample made up of various *Clavibacter* subspecies.


## Folder examples TDA
In this folder we find some scripts with examples of Topological Data Analysis (TDA) libraries and with the aim of applying them to the study of metagenomic and pangemomic data.
1. Alpha_complex
2. Rips_complex
3. Simplex_tree

## Scripts
1. [/scripts/Sim_Metagenome.ipynb](Sim_Metagenome.ipynb)
With this notebook, we use `PyMetaSeem` to generate random reads from various clavibacter michiganensis subspecies stored in [/genomes/6TP/](genomes/6TP/).
 2. make_csv_from_bowtie.ipynb
 In this case, we can generate CSV files, with a Read to Organism Mapping (ROM) between the simulated metagenomes and the database using `bowtie2`.
 3. TDA_Final.ipynb
 In this script, we used the techniques of TDA was proposed by Aldo and eal. to identify the false positives in a group metagenomics data.
 