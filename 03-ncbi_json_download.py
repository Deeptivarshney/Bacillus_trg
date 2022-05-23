import argparse
import json
import os
import urllib.request

# More information at:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/565/GCF_000146565.1_ASM14656v1

CHOICES = [
'assembly_report.txt',
'assembly_stats.txt',
'cds_from_genomic.fna.gz',
'feature_count.txt.gz',
'feature_table.txt.gz',
'genomic.fna.gz',
'genomic.gbff.gz',
'genomic.gff.gz',
'protein.faa.gz',
'protein.gpff.gz',
'translated_cds.faa.gz'
]


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--json', '-json', '-j', dest='json', metavar="FILE",
                    help='JSON file from ncbi_json_genome_assembly.py')
parser.add_argument('--taxid', '-taxid', dest='taxid', type=int,  default=2,
                    help='Taxonomy ID of query group of organisms [DEFAULT: %(default)s]')
parser.add_argument("--filetype", '-f', default='protein.faa.gz', choices = CHOICES,
                       help="File to download [DEFAULT: %(default)s]"),
parser.add_argument('--outdir', '--o', '-o', metavar="DIRECTORY", type=str,
                    help='output directory for downloaded files')
parser.add_argument('--complete_genome', '-complete_genome', action='store_true',
                    help='Include complete genomes')
parser.add_argument('--chromosome', '-chromosome', action='store_true')
parser.add_argument('--scaffold', '-scaffold', action='store_true')
parser.add_argument('--contig', '-contig', action='store_true')
args = parser.parse_args()





# Read JSON file as dictionary
fh = open(args.json)
d = json.load(fh)
fh.close()


levels_allowed = []
if args.complete_genome: levels_allowed.append('Complete Genome')
if args.chromosome: levels_allowed.append('Chromosome')
if args.scaffold: levels_allowed.append('Scaffold')
if args.contig: levels_allowed.append('Contig')


for assembly_id, assembly in d.items():
#for assembly_id in ['001375535']:
    #assembly = d[assembly_id]
    if args.taxid in assembly['lineage']:
        if assembly['assembly_level'] in levels_allowed:
            ftp_dir = assembly['ftp_path']
            ftp_filename_core = assembly['ftp_path'].split('/')[-1]
            ftp_filename = '{}_{}'.format(ftp_filename_core, args.filetype)

            out_genus_dir = os.path.join(args.outdir, str(assembly['genus_taxid']))
            out_species_dir = os.path.join(out_genus_dir, str(assembly['species_taxid_updated']))
            out_assembly_dir = os.path.join(out_species_dir, ftp_filename_core)
            file_download_dir = os.path.join(out_assembly_dir, ftp_filename)
            file_download_md5sum_dir = '{}.md5sum'.format(file_download_dir)

            ncbi_md5sum = None
            uh = urllib.request.urlopen(ftp_dir+'/md5checksums.txt')
            for line in uh:
                line = line.decode()
                if args.filetype in line:
                    ncbi_md5sum = line.split()[0]
            uh.close()

            if ncbi_md5sum:

                if not os.path.exists(out_assembly_dir):
                    os.makedirs(out_assembly_dir)
            
                todownload = False
                if not os.path.exists(file_download_md5sum_dir):
                    todownload = True
                else:
                    fh = open(file_download_md5sum_dir)
                    my_md5sum = fh.read().split()[0]
                    fh.close()
                    
                    if my_md5sum != ncbi_md5sum:
                        todownload = True

                if todownload:
                    print('{}: downloading...'.format(assembly['assembly_accession']))
                    os.system('wget {}/{} -P {} --quiet'.format(ftp_dir, ftp_filename, out_assembly_dir))
                    os.system('md5sum {0} > {0}.md5sum'.format(file_download_dir))
                    os.system('gunzip -f {}'.format(file_download_dir))
                else:
                    print('{}: already exists.'.format(assembly['assembly_accession']))
            
            else:
                print('{} does not have {}.'.format(assembly['assembly_accession'], args.filetype))




                
