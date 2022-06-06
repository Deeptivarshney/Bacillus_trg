## Identification of taxonomically restricted genes (TRGs) in Bacillus 

The workflow is designed to identify taxonomically restricted genes (TRGs) in *Bacillus*. The repository contains several python scripts for downloading and handling sequence data from NCBI,parsing BLAST output, and finding the TRGs at the species and genus levels.

### Requiremnts
1. [Python](https://www.python.org/) >= 3.3
2. [ETE toolkit](http://etetoolkit.org/)
3. [BLAST+ 2.8.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
4. [FASTA-36.3.8g](https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta36/)

## General workflow

1. `01-ncbi_json_genome_assembly.py` provides information related to bacterial genomes available in GenBank and RefSeq and stores it in JSON format. The information includes full taxonomic
    classification and URLs of the sequence(s) and annotation file(s) at the NCBI FTP server. <br /> For example, to fetch the information related to all Bacteria (`taxid : 2`) available in GenBank and RefSeq, run the following command:   <br />
    ```
    python 01-ncbi_json_genome_assembly.py --taxid 2  
    ```
    Output will be a JSON format which stored all the information for all available genomes (ex:b.cereus.json)  <br />
    > Note : Provided example JSON is just to show the format of all information for only 2 genomes for one Bacillus species (taxid:1396) [b.cereus.json](./example/b.cereus.json) . <br />
    <br />
2. `02-ncbi_json_summary.py` provides information about genome assemblies available for a given bacterial taxonomic unit (e.g., *Bacillus*). <br /> For example, to obtain information on the 
    number of available genomes of *Bacillus* (`taxid: 1386`), run the following command: <br />
    ```
    python --json b.cereus.json --taxid 1386 --complete_genome --chromosome --scaffold --contig  
    ```
    This will provide the number of available genomes in all the assembly levels (complete_genome, chromosome, scaffold and contig)  <br />
    <br />
3. `03-ncbi_json_download.py` downloads the sequence data (protein or genome) of a given taxonomic unit of Bacteria.  <br /> For example, to download all proteins of *Bacillus*, run the 
    following command:  <br />
    ```
    python 03-ncbi_json_download.py --json b.cereus.json --taxid 1386 --filetype protein.faa.gz --outdir protein_seq --complete_genome  --chromosome --scaffold --contig
    ```
    In output `ncbi_protein` directory will be created with following structure. For example : [example/protein_seq/1386/1396/GCF_002199365.1_ASM219936v1](./example/protein_seq/1386/1396/GCF_002199365.1_ASM219936v1/)  <br />
    <br />
4. After downloading `*.gz` files, unzip all the files by following command: <br />
    `gzip -d ./example/protein_seq/1386/*/*/*.gz`

5. `04-ncbi_cluster_proteins.py` clusters identical protein sequences within species. To make this pipeline more efficient, protein sequences of all the downloaded genomes under the queried genus
    are clustered into species files based on identical sequences (generated one protein fasta file for respective species which include the all protein sequences from all strains under the same species). <br />
    ``` 
    python --indir ncbi_protein --taxid 1386 --outdir protein_seq_cluster
    ```
    This python script provides one clustered fasta for sequences and a log.json file for `sequence id` and `assembly id` for each protein for every species under queried genus taxid in `protein_seq_cluster` (output folder). The output directory will look like this : [1396.protein.faa](./example/protein_seq_cluster/1386/1396/1396.protein.faa) and  [log.json](./example/protein_seq_cluster/1386/1396/log.json)  <br />
    <br />
6. `05-ncbi_split_proteins.py` splits fasta files into multiple smaller fasta files that can be used as queries to BLAST searches. <br />
    For example, to split fasta files of Bacillus (`taxid: 1386`) into smaller fasta files (1000 sequences per file), run the following command:
    ```
    python 05-ncbi_split_proteins.py --indir protein_seq_cluster --taxid 1386 --outdir protein_seq_split --nseq 1000 
    ```
    The output splitted files will look like this :[1396.protein.faa.001](./example/protein_seq_split/1386/1396/1396.protein.faa.001), [1396.protein.faa.002](./example/protein_seq_split/1386/1396/1396.protein.faa.002) etc. and so on according to the number of sequences.  <br />
    <br />
7. After that, BLAST analysis can be perfomed for all the splitted fasta against the NCBI bacterial proteome by using following command. For example :  <br />
    ```
    blastp -query protein_seq_split/1386/1396/1396.protein.faa.00* \
    -db ncbi_whole_bacteria.protein.faa \
    -outfmt 6 -evalue 10 -num_threads 28 -num_alignments 500 \
    -out protein_seq_split_blast_result/1396.protein.faa.00*
    ```
    > Note : `outfmt 6` format stores the output in tabular output in following directory : ```protein_seq_split_blast_result/1396.protein.faa.00*```  <br />
    <br />
8. `python 06_finding_trg.py` parses BLAST output to identify TRGs. <br /> Assuming you have Blast result files for all protein fasta files: To find the TRGs in *Bacillus*, run the 
    following command: 
    ``` 
    python 06_finding_trg.py  --fastadir protein_seq_split --taxid 1386 --blastdir protein_seq_split_blast_result --evalue 10
    ```
    As mentioned in the command line for taxid (1386: Bacillus), this python script generates the number of sequences for those that have no homology at provided value (10) on the Genus level. 
    <br />

9. Reciprocal Blast search has also implemented this pipeline, A High scorer (bit score) hit was considered as the best hit in both searches. 
    ```
    blastp -query protein_seq_split/1386/1396/1396.protein.faa.00* \
    -db ncbi_whole_bacteria.protein.faa \
    -outfmt 6 -evalue 10 -num_threads 28 -num_alignments 20 \
    -out forwrd_blast_result/1396.protein.faa.00*
    ```
    Since, to make the pipeline to more efficient, protein sequences of query and database clustered based on identical on species level. After performing the forward blast search, strain ID can be added for each hit, and then can be parsed by this script <br />
    ```
    python 07_forwrd_reciprocal_hit_parse.py --indir forwrd_blast_result --outdir blastp_fwd_hit_result
    ```
    This python script gives the tab a separated file for each blast result which contains the following information : <br />
    
   
    | queryid |sub_genus_taxid |sub_species_taxid | sub_genome | subject_besthit_ids|
    | --- | --- | --- | --- |--- |
    
    * Same way, High scrorer hit can be query to back the query proteome by using the above mentioned blast command.
    
10. After both sides of BLAST searches, reciprocal hits can be detected by running this script :<br />
    ```
     python 08_reciprocal_hits.py --qgff bacillus_gff --fblast blastp_fwd_hit_result --rblast reverse_blastp_results --taxid 1386 --outdir reciprocal_hits
    ```
11. tFASTy search performed for idenitified putative TRGs against the non coding region of other NCBI Bacterial genomes.
    For example : tFASty search for a fasta sequence can be performed by follwing command (assuming `ncbi_genome_assemblies.faa` as database) :

    ```tfasty36 -m 0 -w 500 -C 50 -T 24 -s BL62 seq.faa ncbi_genome_assemblies.faa >seq.tfasty```


