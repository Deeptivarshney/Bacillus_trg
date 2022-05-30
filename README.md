## Manuscript title 

The workflow is designed to identify taxonomically restriced genes (TRGs) in Bacillus. There are several python scripts in term of finding TRG, which include downloading all the information (eg:proteome,genome and annotation) from NCBI as well finding the TRGs.

### Requiremnts
Python (https://www.python.org/) >= 3.3

ETE toolkit (http://etetoolkit.org/)

## General workflow

1. Before downloding the NCBI data, the information related to Bacterial genomes available in GenBank and RefSeq is organized into easily parsable JSON files, which includes information such as updated taxonomical classification (at all taxonomical levels), URLs for the sequence(s) and annotation file(s).

For example : To fetch the information of all Bacteria available in NCBI, run following command : 

```
python 01-ncbi_json_genome_assembly.py --taxid 2 
```
This python script takes taxonomic identifier as an arguments and provide the output in JSON format (taxid:2 for whole bacteria)

2. To get the information about assemblies level for all the genomes queried taxonomy ID from the JSON file

```
python 02-ncbi_json_summary.py  --json ncbi.json --taxid 1386 --complete_genome --chromosome --scaffold --contig
```
This python script provides the information about the number of available assemblies for each level for Bacillus (taxid : 1386). 


3. To download the sequences for required file type for available genomes by using JSON File. 

```
python 03-ncbi_json_download.py --json ncbi.json --taxid 1386 --filetype protein.faa.gz --outdir ncbi_protein --complete_genome  --chromosome --scaffold --contig
```
This python scripts gives you thechoice of filetype to download.

4. To make this pipeline more efficient, protein sequences of all the downloded genomes under the queried genus were clustered into species file based on dentical sequences (generated one protein fasta file for respective species which include the all protein sequences from all strains under the same species)
``` 
python 04-ncbi_cluster_proteins.py --indir ncbi_protein --taxid 1386 --outdir ncbi_protein_cluster
```
This python script provides one clustered fasta for sequences and log.json file for sequence id and assembly id for the each protein for every species under queried genus taxid in ncbi_protein_cluster (output folder)

5. To efficient the BLAST search for identification, one fasta file split into multi fasta files 
```
python 05-ncbi_split_proteins.py --indir ncbi_protein_cluster --taxid 1386 --outdir ncbi_protein_clustesplit --nseq 1000 
``` 

This python script splits 1000 number of sequences per file. 

6. After splitting the sequences into multiple files, perform BLAST search and parse the blast results 

``` 
python 06_finding_trg.py  --fastadir ncbi_protein_clustesplit --taxid 1386 --blastdir ncbi_blast --evalue 10
```
As mentioned in the command line for taxid (1386 : Bacillus), this python script generates the number of sequences for those have no homology at provided evalue (10) on Genus level. 


7. Reciprocal Blast search was also implemented this pipeline, A High scorer (bit score) hit was considered as a best hit in both search. After performing forword blast search, add the strainID for each hit and then can be parsed by this script

```
python 07_forwrd_reciprocal_hit_parse.py --indir forwrd_blast_result --outdir blastp_fwd_hit_result

```
This python script gives the tab seperated file for each blast result which contains following information : 
``` 
queryid,sub_genus_taxid,sub_species_taxid,sub_genome,subject_besthit_ids
```
8. After both side of BLAST searches, reciprocal hits can be detected by running this script :

```
python 08_reciprocal_hits.py --qgff bacillus_gff --fblast blastp_fwd_hit_result --rblast reverse_blastp_results --taxid 1386 --outdir reciprocal_hits

```







