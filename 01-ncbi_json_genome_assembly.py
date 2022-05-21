import argparse
import json
import urllib.request

from ete3 import NCBITaxa

# More information about genomes:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_assembly_summary.txt

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--taxid', '-taxid', dest='taxid', type=int, default=2,
                    help='Taxonomy ID of query group of organisms [DEFAULT: %(default)s]')
args = parser.parse_args()


URL_ROOT = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_{}.txt'
DATABASES = ['refseq', 'genbank']

ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

assemblies = {}
for DB in DATABASES:
    url = URL_ROOT.format(DB)
    uh = urllib.request.urlopen(url)
    print('Processing: {}'.format(url))
    for line in uh:
        line = line.decode()
        # Skip comment lines.
        if not line.startswith('#'):
            sl = line.strip().split('\t')
            assembly_accession = sl[0]
            assembly_version = int(assembly_accession.split('.')[1])
            assembly_id = assembly_accession[4:assembly_accession.index('.')]
            taxid = int(sl[5].strip())
            species_taxid = int(sl[6].strip())
            organism_name = sl[7].strip()
            infraspecific_name = sl[8].strip()
            version_status = sl[10].strip()
            assembly_level = sl[11].strip()
            asm_name = sl[15].strip()
            seq_rel_date = sl[14]
            ftp_path = sl[19]
            lineage = ncbi.get_lineage(species_taxid)
            # If organism belongs to GROUP_TAXID
            if args.taxid in lineage:
                if assembly_id not in assemblies:
                    names = ncbi.get_taxid_translator(lineage)
                    ranks = ncbi.get_rank(lineage)
                    lineage_names = [names[taxid] for taxid in lineage]
                    lineage_ranks = [ranks[taxid] for taxid in lineage]
                    species_taxid_updated = None
                    genus_taxid = None
                    for tid, name, rank in zip(lineage, lineage_names, lineage_ranks):
                        if rank == 'species':
                            species_taxid_updated = tid
                        elif rank == 'genus':
                            genus_taxid = tid

                    assemblies[assembly_id] = {
                        'assembly_accession': assembly_accession,
                        'assembly_version': assembly_version,
                        'asm_name': asm_name,
                        'assembly_level': assembly_level,
                        'ogranism_name': organism_name,
                        'taxid': taxid,
                        'species_taxid': species_taxid,
                        'species_taxid_updated': species_taxid_updated,
                        'genus_taxid': genus_taxid,
                        'version_status': version_status,
                        'seq_rel_date': seq_rel_date,
                        'ftp_path': ftp_path,
                        'lineage': lineage,
                        'lineage_names': lineage_names,
                        'lineage_ranks': lineage_ranks,
                    }
    uh.close()


# Get name of query taxid
taxid2name = ncbi.get_taxid_translator([args.taxid])
name = taxid2name[args.taxid]
# Save to json
oh = open('{}.{}.json'.format(name, args.taxid), 'w')
json.dump(assemblies, oh, indent=3)
oh.close()
