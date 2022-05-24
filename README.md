## Manuscript title 

The workflow is designed to identification of taxonomically restriced genes (TRGs) in Bacillus. There are sevral python scripts in term of finding TRG, which include downloading all the information (eg:proteome,genome and annotation) from NCBI as well 
finding the TRGs. 


### Requiremnts 

Python > 3.0

ETE toolkit

## General workflow

1. The information related to Bacterial genomes available in GenBank and RefSeq is organized into easily parsable JSON files, including information such as updated taxonomical classification (at all taxonomical levels), URLs for the sequence and annotation files. 

```
python 01-ncbi_json_genome_assembly.py --taxid 2 
```
This python script takes taxonomic identifier as an arguments and provide the output in JSON format (taxid:2 for whole bacteria)

2. To get the information about assemblies level for all the genomes queried taxonomy ID from the JSON file

```
python 02-ncbi_json_summary.py  --json ncbi.json --taxid 1386 --complete_genome --chromosome --scaffold --contig
```

This python script provide the information about the number of assemblies for each level. 

