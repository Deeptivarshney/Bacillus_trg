import argparse
import json
import os

from Bio import SeqIO


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--indir', '--i', '-i', metavar="DIRECTORY", type=str, default="ncbi_bacteria",
                    help='input directory with genus subdirectories [DEFAULT: %(default)s]')
parser.add_argument('--taxid', '-taxid', dest='taxid', type=int,  default=1386,
                    help='Taxonomy ID of query group of organisms [DEFAULT: %(default)s]')
parser.add_argument('--outdir', '--o', '-o', metavar="DIRECTORY", type=str, default="ncbi_bacteria_cluster",
                    help='output directory to save results [DEFAULT: %(default)s]')
args = parser.parse_args()




IN_GENUS_PATH = os.path.join(args.indir, str(args.taxid))
if not os.path.exists(IN_GENUS_PATH):
    print('Input directory {} does not exist.'.format(IN_GENUS_PATH))
    exit()


OUT_GENUS_PATH = os.path.join(args.outdir, str(args.taxid))
if not os.path.exists(OUT_GENUS_PATH):
    os.makedirs(OUT_GENUS_PATH)


N = 0
K = 0
A = 0
S = 0
for species_taxid in os.listdir(IN_GENUS_PATH):
    in_species_path = os.path.join(IN_GENUS_PATH, species_taxid)
    seq2num = {}
    d = {}
    k = 0
    n = 0
    a = 0
    for assembly_dir in os.listdir(in_species_path):
        assembly_path = os.path.join(in_species_path, assembly_dir)
        filename = '{}_protein.faa'.format(assembly_dir)
        filename_path = os.path.join(assembly_path, filename)
        fh = open(filename_path)
        for seq_record in SeqIO.parse(fh, 'fasta'):
            seq = str(seq_record.seq)
            n += 1
            if seq not in seq2num:
                k +=1
                num = 'taxid_{}_{:06}'.format(species_taxid, k)
                seq2num[seq] = num
                d[num] = []
            else:
                num = seq2num[seq]
            d[num].append((seq_record.id, assembly_dir))
        fh.close()
        a += 1

    N += n
    K += k
    A += a
    out_species_path = os.path.join(OUT_GENUS_PATH, species_taxid)
    if not os.path.exists(out_species_path):
        os.makedirs(out_species_path) 

    out_fasta_path = os.path.join(out_species_path, '{}.protein.faa'.format(species_taxid))
    oh = open(out_fasta_path, 'w')
    for seq, num in seq2num.items():
        oh.write('>{}\n{}\n'.format(num, seq))
    oh.close()

    out_log_path = os.path.join(out_species_path, 'log.json')
    oh = open(out_log_path, 'w')
    json.dump(d, oh, indent=2)
    oh.close()
    S += 1
    print('species_taxid:{}\tassemblies:{}\t{:,}\t{:,}\t{:.1f}'.format(species_taxid, a, n, k, (n-k)/n*100))
print('SPECIES:{}\tASSEMBLIES:{}\t{:,}\t{:,}\t{:.1f}'.format(S, A, N, K, (N-K)/N*100))