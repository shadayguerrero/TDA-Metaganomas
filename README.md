# TDA-Metaganomas
In this repository we have some tools that we will use to distinguish the presence or absence of species in a metagenome made up of closely related species. In particular, we want to determine the presence of *Clavibacter michiganensis michiganensis* in a sample made up of various *Clavibacter* subspecies.


## Folder examples TDA
In this folder we find some scripts with examples of Topological Data Analysis (TDA) libraries and with the aim of applying them to the study of metagenomic and pangemomic data.
1. Alpha_complex
2. Rips_complex
3. Simplex_tree

## Scripts
1. [Sim_Metagenome.ipynb](/scripts/Sim_Metagenome.ipynb)

With this notebook, we use `PyMetaSeem` to generate random reads from various clavibacter michiganensis subspecies stored in [/genomes/6TP/](genomes/6TP/).

 2. [make_csv_from_bowtie.ipynb](/scripts/make_csv_from_bowtie.ipynb)
 
 In this case, we can generate CSV files, with a Read to Organism Mapping (ROM) between the simulated metagenomes and the database using `bowtie2`.
 
 3. [TDA_Final.ipynb](/scripts/TDA_Final.ipynb)
 
 In this script, we used the techniques of TDA was proposed by Aldo and eal. to identify the false positives in a group metagenomics data.

## genomes

1.[Clavibacter](/genomes/6TP/)
~~~
Clavibacter_michiganensis_subsp_capsici_1101.fna
Clavibacter_michiganensis_subsp_insidiosus_ATCC_10253.fna
Clavibacter_michiganensis_subsp_michi_contigs.fna
Clavibacter_michiganensis_subsp_nebraskensis_419B.fna
Clavibacter_michiganensis_subsp_sepedonicus_ATCC33113.fna
Clavibacter_michiganensis_subsp_tessellarius_ATCC_33566.fna
~~~

2. [Salmonella](/genomes/salmonella/)
~~~
CP019409.fasta  CP019417.fasta  CP022135.fasta  CP022494.fasta
CP019116.fasta  CP019410.fasta  CP019418.fasta  CP022138.fasta  CP022497.fasta
CP019403.fasta  CP019411.fasta  CP022015.fasta  CP022139.fasta  CP022500.fasta
CP019404.fasta  CP019412.fasta  CP022019.fasta  CP022142.fasta  CP022503.fasta
CP019405.fasta  CP019413.fasta  CP022034.fasta  CP022467.fasta  CP022504.fasta
CP019406.fasta  CP019414.fasta  CP022116.fasta  CP022489.fasta  
CP019407.fasta  CP019415.fasta  CP022117.fasta  CP022490.fasta
CP019408.fasta  CP019416.fasta  CP022120.fasta  CP022491.fasta

~~~