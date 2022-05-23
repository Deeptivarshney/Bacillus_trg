import os
import argparse


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--qgff', '--gf', '-gf', metavar="DIRECTORY", type=str, default="bacilus_gff",
                    help='input directory with query genus with genome subdirectories for GFF [DEFAULT: %(default)s]')
parser.add_argument('--fblast', '--fb', '-fb', metavar="DIRECTORY", type=str, default="bacillus_fwd_blastphit",
                    help='input directory forwrd blast homologs for each species [DEFAULT: %(default)s]')
parser.add_argument('--rblast', '--rb', '-rb', metavar="DIRECTORY", type=str, default="bacillus_rev_blastphit",
                    help='input directory reverse blast homologs for each species [DEFAULT: %(default)s]')
parser.add_argument('--outdir', '--o', '-o', metavar="DIRECTORY", type=str, default="reciprocal_blast_result",
                    help='output directory to save results [DEFAULT: %(default)s]')
parser.add_argument('--taxid', '-taxid', dest='taxid', type=int,  default=1396,
                    help='Taxonomy ID   [DEFAULT: %(default)s]')
args = parser.parse_args()

QUERY_GFF_PATH = args.qgff

def check_fwd_blast(trg_genome_uniq_id, trg_genome_original_id, trg_genome, other_genome_protein_id):
    combine_id_for_lookup = trg_genome_uniq_id, trg_genome
    if combine_id_for_lookup in homolog:
        for blast_result in homolog[combine_id_for_lookup]:
            if blast_result[-1] == other_genome_protein_id:
                out_file_path.write("{}\t{}\t{}\t{}\n".format(trg_genome_uniq_id, trg_genome_original_id, trg_genome, "\t".join(blast_result[2:])))
            
    


species_id = args.taxid
original_to_uniq = {}

IN_GENUS_PATH = os.path.join(QUERY_GFF_PATH, str(args.taxid))

in_species = os.path.join(IN_GENUS_PATH, species_id)
    
for gffassembly_dir in os.listdir(in_species):
    gffassembly_path = os.path.join(in_species,gffassembly_dir)
    
    gffassembly_name = '{}_genomic_gff.txt'.format(gffassembly_dir)
    gff_assembly_name_path = os.path.join(gffassembly_path, gffassembly_name)
    
    gffopen = open(gff_assembly_name_path)
    for line in gffopen:
        sl = line.rstrip().split("\t")
        uniq_protein_id = sl[0]
        original_id = sl[5]
        
        if gffassembly_dir not in original_to_uniq:
            original_to_uniq[gffassembly_dir]= {}
        if original_id not in original_to_uniq[gffassembly_dir]:
            original_to_uniq[gffassembly_dir][original_id] = uniq_protein_id
         
          
    gffopen.close()


out_path = args.outdir
os.makedirs(out_path, exist_ok= True)

out_file_path = open(os.path.join(out_path, "{}.rbh.out".format(species_id)), "w")

homolog = {}
forwrd_blast_path = args.fblast

for files in os.listdir(forwrd_blast_path):
    if files.startswith("{}".format(species_id)):

        in_filepath = os.path.join(forwrd_blast_path, files)
        fh = open(in_filepath, "r")
        data = fh.readlines()[1:]
        for line in data:
            spl = line.rstrip().split("\t")
            trg_genome_protein_id = spl[0]
            trg_genome = spl[1]

            homo_combine_key = trg_genome_protein_id, trg_genome

            if homo_combine_key not in homolog:
                homolog[homo_combine_key] = []
            homolog[homo_combine_key].append(spl)
        fh.close()


reversed_best_hit_path = args.rblast
for result_file in os.listdir(reversed_best_hit_path):
    result_file_path = os.path.join(reversed_best_hit_path, result_file)
    
    rev_blast_fh = open(result_file_path, "r")
    ignore_first_line = rev_blast_fh.readlines()[1:]

    for blast_line in ignore_first_line:
        blast_ln_spl = blast_line.rstrip().split("\t")
        
        protein_query_id_from_other_genome = blast_ln_spl[0]
        trg_genome_as_subject = blast_ln_spl[1]
        protein_original_id_as_subject_from_trg_genomes = blast_ln_spl[2].split(",")

        for ech_id in protein_original_id_as_subject_from_trg_genomes:
            if ech_id in original_to_uniq[trg_genome_as_subject]:
                result = check_fwd_blast(original_to_uniq[trg_genome_as_subject][ech_id], ech_id, trg_genome_as_subject, protein_query_id_from_other_genome)
out_file_path.close()   

