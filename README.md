## Manuscript title 

The workflow is designed to identification of taxonomically restriced genes (TRGs) in Bacillus. There are sevral python scripts in term of finding TRG, which include downloading all the information (eg:proteome,genome and annotation) from NCBI as well 
finding the TRGs. 


### Requiremnts 
python:3
ETE toolkit

## General workflow

1. The information related to Bacterial genomes available in GenBank and RefSeq is organized into easily parsable JSON files, including information such as updated taxonomical classification (at all taxonomical levels), URLs for the sequence and annotation files. 

```
python 01-ncbi_json_genome_assembly.py --taxid 2 
```
This python script taxonomic identifier as an arguments and provide the output in JSON format. 
```

Before downloading the protein sequences for queried genus/species, information of more than
bacterial genome assemblies available in GenBank and RefSeq is organized into easily parsable JSON files, 
including information such as updated taxonomical classification (at all taxonomical levels), URLs for 
the sequence and annotation files.

```python scripts

01-ncbi_json_genome_assembly.py
02-ncbi_json_summary.py
03-ncbi_json_download.py
```
Step 2: Clustered the protein sequences on species level (generated one protein fasta file for respective
species which include the all protein sequences from all strains under the same species)

```python scripts

04-ncbi_cluster_proteins.py
```

Step 3: Split the protein sequences file into different subfiles (1000 seqs in each file) for performing the
blastp analysis against Uniprot Reference database

```python scripts

05-ncbi_split_proteins.py
```