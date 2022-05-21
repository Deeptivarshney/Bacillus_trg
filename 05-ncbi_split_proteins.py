import argparse
import json
import os
import shutil
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--indir', '--i', '-i', metavar="DIRECTORY", type=str, default="ncbi_bacteria_cluster",
                    help='input directory with genus subdirectories [DEFAULT: %(default)s]')
parser.add_argument('--taxid', '-taxid', dest='taxid', type=int,  default=1386,
                    help='Taxonomy ID of query group of organisms [DEFAULT: %(default)s]')
parser.add_argument('--outdir', '--o', '-o', metavar="DIRECTORY", type=str, default="ncbi_bacteria_split",
                    help='output directory to save results [DEFAULT: %(default)s]')
parser.add_argument('--nseq', '--n', '-n', type=int, default=1000,
                    help='number of sequences per file [DEFAULT: %(default)s]')
args = parser.parse_args()


IN_GENUS_PATH = os.path.join(args.indir, str(args.taxid))
if not os.path.exists(IN_GENUS_PATH):
    print('Input directory {} does not exist.'.format(IN_GENUS_PATH))
    exit()


OUT_GENUS_PATH = os.path.join(args.outdir, str(args.taxid))
if os.path.exists(OUT_GENUS_PATH):
    # Remove existing splitted sequences
    shutil.rmtree(OUT_GENUS_PATH)
os.makedirs(OUT_GENUS_PATH)


N = 0
for species_taxid in os.listdir(IN_GENUS_PATH):
    in_species_path = os.path.join(IN_GENUS_PATH, species_taxid)
    in_filename = '{}.protein.faa'.format(species_taxid)
    in_filename_path = os.path.join(in_species_path, in_filename)

    fh = open(in_filename_path)
    out_species_path = os.path.join(OUT_GENUS_PATH, species_taxid)
    if not os.path.exists(out_species_path):
        os.makedirs(out_species_path)

    n = 1
    out_filename = '{}.{:03}'.format(in_filename, n)
    oh = open(os.path.join(out_species_path, out_filename), 'w')
    for i, seq_record in enumerate(SeqIO.parse(in_filename_path, 'fasta')):
        if i and i % args.nseq == 0:
            oh.close()
            n += 1
            out_filename = '{}.{:03}'.format(in_filename, n)
            oh = open(os.path.join(out_species_path, out_filename), 'w')
        oh.write(seq_record.format('fasta'))
    N += n
            
    oh.close()

    fh.close()
print('Number of files: {}'.format(N))