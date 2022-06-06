## Identification of taxonomically restricted genes (TRGs) in Bacillus 

The workflow is designed to identify taxonomically restricted genes (TRGs) in *Bacillus*. The repository contains several python scripts for downloading and handling sequence data from NCBI, parsing BLAST output, and finding the TRGs at the species and genus levels.

### Requiremnts
1. Python (https://www.python.org/) >= 3.3
2. ETE toolkit (http://etetoolkit.org/)
3. BLAST+ 

## General workflow

1. `01-ncbi_json_genome_assembly.py` provides information related to bacterial genomes available in GenBank and RefSeq and stores it in JSON format. The information includes full taxonomic classification and URLs of the sequence(s) and annotation file(s) at the NCBI FTP server.

For example, to fetch the information related to all Bacteria (`taxid : 2`) available in GenBank and RefSeq, run the following command: 

```
python 01-ncbi_json_genome_assembly.py --taxid 2 
```

Output will be a JSON format which stored all the information for all available genomes (ex:b.cereus.json) 
*Note: Provided example JSON is just to show the format of all information for only 2 genomes for one Bacillus species (taxid:1396)

2. `02-ncbi_json_summary.py` provides information about genome assemblies available for a given bacterial taxonomic unit (e.g., *Bacillus*).

For example, to obtain information on the number of available genomes of *Bacillus* (`taxid: 1386`), run the following command:

```
python --json b.cereus.json --taxid 1386 --complete_genome --chromosome --scaffold --contig
```

This will provide the number of available genomes in all the assembly levels(complete_genome,chromosome,scaffold and contig)


3. `03-ncbi_json_download.py` downloads the sequence data (protein or genome) of a given taxonomic unit of Bacteria.

For example, to download all proteins of *Bacillus*, run the following command: 

```
python 03-ncbi_json_download.py --json b.cereus.json --taxid 1386 --filetype protein.faa.gz --outdir protein_seq --complete_genome  --chromosome --scaffold --contig
```

In output `ncbi_protein` directory will be created with following structure. For example : `protein_seq/1386/1396/GCF_XXXXX/GCF_XXXXX.protein.gz`


4. `04-ncbi_cluster_proteins.py` clusters identical protein sequences within species. To make this pipeline more efficient, protein sequences of all the downloaded genomes under the queried genus are clustered into species files based on identical sequences (generated one protein fasta file for respective species which include the all protein sequences from all strains under the same species).

``` 
python --indir ncbi_protein --taxid 1386 --outdir protein_seq_cluster
```
This python script provides one clustered fasta for sequences and a log.json file for `sequence id` and `assembly id` for each protein for every species under queried genus taxid in `protein_seq_cluster` (output folder)

The output directory will look like this : `protein_seq_cluster/1386/1392/1392.protein.faa` and  `protein_seq_cluster/1386/1392/log.json` 


5. `05-ncbi_split_proteins.py` splits fasta files into multiple smaller fasta files that can be used as queries to BLAST searches.

For example, to split fasta files of Bacillus (`taxid: 1386`) into smaller fasta files (1000 sequences per file), run the following command:

```
python 05-ncbi_split_proteins.py --indir ncbi_protein_cluster --taxid 1386 --outdir ncbi_protein_clustesplit --nseq 1000 
```

The ouput directory will look like this : `protein_seq_split/1386/1392/1392.protein.faa.001`, `protein_seq_split/1386/1396/1396.protein.faa.00*` and so on according to the number of sequences. 

After that, BLAST analysis can be perfomed for all the splitted fasta against the NCBI bacterial proteome by using following command. For example :
 
```blastp -query protein_seq_split/1386/1396/1396.protein.faa.00* -db ncbi_whole_bacteria.protein.faa -outfmt 6 -evalue 10 -num_threads 28 -num_alignments 500 -out protein_seq_split_blast_result/1396.protein.faa.00*```

* Note : `outfmt 6` format stores the output in tabular output in following directory : ```protein_seq_split_blast_result/1396.protein.faa.00*```

6. `python 06_finding_trg.py` parses BLAST output to identify TRGs. 
For example, to find the TRGs in *Bacillus*, run the following command: 

``` 
python 06_finding_trg.py  --fastadir ncbi_protein_clustesplit --taxid 1386 --blastdir ncbi_blast --evalue 10
```
As mentioned in the command line for taxid (1386: Bacillus), this python script generates the number of sequences for those that have no homology at provided value (10) on the Genus level. 


7. Reciprocal Blast search has also implemented this pipeline, A High scorer (bit score) hit was considered as the best hit in both searches. After performing the forward blast search, add the strain ID for each hit, and then can be parsed by this script

```
python 07_forwrd_reciprocal_hit_parse.py --indir forwrd_blast_result --outdir blastp_fwd_hit_result

```
This python script gives the tab a separated file for each blast result which contains the following information : 
``` 
queryid,sub_genus_taxid,sub_species_taxid,sub_genome,subject_besthit_ids
```
8. After both sides of BLAST searches, reciprocal hits can be detected by running this script :

```
python 08_reciprocal_hits.py --qgff bacillus_gff --fblast blastp_fwd_hit_result --rblast reverse_blastp_results --taxid 1386 --outdir reciprocal_hits

```