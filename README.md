## Manuscript title 

The workflow is designed to identification of taxonomically restriced genes (TRGs) in Bacillus. There are several python scripts in term of finding TRG, which include downloading all the information (eg:proteome,genome and annotation) from NCBI as well finding the TRGs.

### Requiremnts
Python (https://www.python.org/) >= 3.3

ETE toolkit (http://etetoolkit.org/)

## General workflow

1. The information related to Bacterial genomes available in GenBank and RefSeq is organized into easily parsable JSON files, which includes information such as updated taxonomical classification (at all taxonomical levels), URLs for the sequence and annotation files.s. 

```
python 01-ncbi_json_genome_assembly.py --taxid 1386 
```
This python script takes taxonomic identifier as an arguments and provide the output in JSON format (taxid:2 for whole bacteria)

2. To get the information about assemblies level for all the genomes queried taxonomy ID from the JSON file

```
python 02-ncbi_json_summary.py  --json ncbi.json --taxid 1386 --complete_genome --chromosome --scaffold --contig
```
This python script provides the information about the number of assemblies for each level. 

3. To download all the sequences for all the genomes/proteome for queried taxonomic ID. 

```
python 03-ncbi_json_download.py --json ncbi.json --taxid 1386 --filetype protein.faa.gz --outdir ncbi_protein --complete_genome  --chromosome --scaffold --contig
```
4. To make this pipeline more efficient, protein sequences of all the downloded genomes under the queried genus were clustered into species file based on dentical sequences (generated one protein fasta file for respective species which include the all protein sequences from all strains under the same species)
``` 
python 04-ncbi_cluster_proteins.py --indir ncbi_protein --taxid 1386 --outdir ncbi_protein_cluster
```

This python script provides clustered fasta files for each species for queried genus taxid in ncbi_protein_cluster folder 

5. To efficient the BLAST search for identification, one fasta file split into multi fasta files 
```
python 05-ncbi_split_proteins.py --indir ncbi_protein_cluster --taxid 1386 --outdir ncbi_protein_clustesplit --nseq 1000 
``` 

This python script splits 1000 number of sequences per file. 

6. After splitting the sequences into multiple files, perform BLAST search and parse the blast results 

``` 
python 06_finding_trg.py  --fastadir ncbi_protein_clustesplit --taxid 1386 --blastdir ncbi_blast --evalue 10
```
 



