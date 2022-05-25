import os 
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--indir', '--i', '-i', metavar="DIRECTORY", type=str, default="Blast_result",
                    help='input directory for blast result [DEFAULT: %(default)s]')

parser.add_argument('--outdir', '--o', '-o', metavar="DIRECTORY", type=str, default="tfasty_result",
                    help='output directory to save blast hit [DEFAULT: %(default)s]')

args = parser.parse_args()


out_dir = os.path.join(args.outdir, "1386")
os.makedirs(out_dir, exist_ok=True)

Genus_taxid = "1386"

bacillus_result_path = args.indir
if not os.path.exists(bacillus_result_path):
    exit()

for blast_files in os.listdir(bacillus_result_path):
    blast_file_path = os.path.join(bacillus_result_path, blast_files)
    print(blast_files)

    out_blast_file_path = os.path.join(out_dir, blast_files)

    oh = open(out_blast_file_path, "w")
    oh.write("queryid\tsub_genus_taxid\tsub_species_taxid\tsub_genome\tsubject_besthit_ids\n")

    best_hit = {}
    fh = open(blast_file_path)
    for line in fh:
        sl = line.rstrip().split("\t")
        query_id = str(sl[0])
        sub_id = str(sl[1])
        sub_genome_id = sub_id.split("|")[2]
        score = float(sl[-1])
        
        """ take the best hit from every genome outside the genus"""
        combine_key = query_id, sub_genome_id 
        if combine_key not in best_hit:
            best_hit[combine_key] = set()
        best_hit[combine_key].add((sub_id, score))

    fh.close()
    
    for BLAST_QUERY, BLAST_HITS in best_hit.items():
        SORTED_LIST = sorted(BLAST_HITS,key=lambda x: x[1], reverse=True)
        res = float(max(SORTED_LIST, key = lambda i : i[1])[1])
        hits_list = set([v[0] for i, v in enumerate(BLAST_HITS) if v[1] == res])
        oh.write("{}\t{}\t{}\t{}\t{}\n".format(BLAST_QUERY[0], str(hits_list).split("|")[0].replace("{'", ""), str(hits_list).split("_")[1], BLAST_QUERY[1], ",".join(["|".join(i.split("|")[:2]) for i in hits_list])))
    
    oh.close()  
