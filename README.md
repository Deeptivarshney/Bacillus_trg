## General workflow


Step 1: Download the available protein sequences of all species (all different strains) under the queried
genus.

Before downloading the protein sequences for queried genus/species, information of more than
150,000 bacterial genome assemblies (~37,760 species) available in GenBank and RefSeq is
organized into easily parsable JSON files, including information such as updated taxonomical
classification (at all taxonomical levels), URLs for the sequence and annotation files.

```python scripts

ncbi_json_genome_assembly.py
download_from_ensembl.py
download_from_ncbi.py
```
Step 2: Clustered the protein sequences on species level (generated one protein fasta file for respective
species which include the all protein sequences from all strains under the same species)

```python scripts

ncbi_cluster_proteins.py
```

Step 3: Split the protein sequences file into different subfiles (1000 seqs in each file) for performing the
blastp analysis against Uniprot Reference database

```python scripts

ncbi_split_proteins.py
```