import os
import argparse


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--fastadir', '--f', '-f', metavar="DIRECTORY", type=str, default="ncbi_bacteria_split",
                    help='input directory for fasta file [DEFAULT: %(default)s]')
parser.add_argument('--taxid', '-taxid', dest='taxid', type=int,  default=1386,
                    help='Taxonomy ID of genus [DEFAULT: %(default)s]')
parser.add_argument('--blastdir', '--b', '-b', metavar="DIRECTORY", type=str, default="ncbi_bacteria_blast_result",
                    help='input directory with blast result [DEFAULT: %(default)s]')
parser.add_argument('--evalue', '--e', '-e', type=float, default=10.0,
                    help='evalue for parsing the blast result [DEFAULT: %(default)s]')                                      
args = parser.parse_args()


GENUS_TAXID = args.taxid
FASTA_DIR = '{}/{}'.format(args.fastadir,GENUS_TAXID)
BLAST_DIR = args.blastdir
EVALUE = args.evalue


n_trg_genus = 0    # Number of genus-level TRG proteins in Bacillus.
n_all_genus = 0    # Number of all proteins in Bacillus.
oh = open('log.genus.txt', 'w')
oh.write('#species_taxid\ttrg\ttotal\tpercent\n')
for species_taxid in os.listdir(FASTA_DIR):
    species_fasta_path = os.path.join(FASTA_DIR, species_taxid)
    n_trg = 0
    n_all = 0
    for fasta_name in os.listdir(species_fasta_path):
        # Get seqids from a splitted FASTA file
        all_proteins = set()        
        fh = open(os.path.join(species_fasta_path, fasta_name))
        for line in fh:
            if line.startswith('>'):
                seqid = line.split()[0].lstrip('>')
                all_proteins.add(seqid)
        fh.close()
        n_all += len(all_proteins)

        # BLAST results for the splitted FASTA file.
        blast_name = '{}_blastp.txt'.format(fasta_name.replace('.faa', ''))
        fh = open(os.path.join(BLAST_DIR, blast_name))
        homologs = set([])
        for line in fh:
            sl = line.split()
            qid = sl[0]
            sid = sl[1]
            evalue = float(sl[10])
            if evalue < EVALUE:
                sid_genus_taxid = sid.split('|')[0]
                if GENUS_TAXID != sid_genus_taxid:
                    homologs.add(qid)
        fh.close()

        # Get TRGs.
        # As difference between a set of all protein ids and
        # a set of proteins that have homologs.
        trgs = all_proteins.difference(homologs)
        n_trg += len(trgs)
    
    oh.write('{}\t{}\t{}\t{:.3f}\n'.format(species_taxid, n_trg, n_all, n_trg/n_all*100))
    n_trg_genus += n_trg
    n_all_genus += n_all
oh.close()
print(n_trg_genus, n_all_genus, n_trg_genus/n_all_genus*100)