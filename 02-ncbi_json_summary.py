import argparse
import json

# More information about genomes:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_assembly_summary.txt

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--json', '-json', '-j', dest='json', metavar="FILE",
                    help='JSON file from ncbi_json_genome_assembly.py')
parser.add_argument('--taxid', '-taxid', dest='taxid', type=int,  default=2,
                    help='Taxonomy ID of query group of organisms [DEFAULT: %(default)s]')
parser.add_argument('--complete_genome', '-complete_genome', action='store_true',
                    help='''all chromosomes are gapless and have no runs of 10 or more ambiguous 
                    bases (Ns), there are no unplaced or unlocalized scaffolds, and all the expected 
                    chromosomes are present (i.e. the assembly is not noted as having partial genome 
                    representation). Plasmids and organelles may or may not be included in the assembly 
                    but if present then the sequences are gapless.''')
parser.add_argument('--chromosome', '-chromosome', action='store_true',
                    help='''there is sequence for one or more chromosomes. This could be a completely 
                    sequenced chromosome without gaps or a chromosome containing scaffolds or contigs 
                    with gaps between them. There may also be unplaced or unlocalized scaffolds.'''
                    )
parser.add_argument('--scaffold', '-scaffold', action='store_true',
                    help='''some sequence contigs have been connected across gaps to create scaffolds, 
                    but the scaffolds are all unplaced or unlocalized.''')
parser.add_argument('--contig', '-contig', action='store_true',
                    help='''nothing is assembled beyond the level of sequence contigs''')
args = parser.parse_args()


# Read JSON file as dictionary
fh = open(args.json)
d = json.load(fh)
fh.close()


levels_allowed = {}
if args.complete_genome: levels_allowed['Complete Genome'] = [0, set([])]
if args.chromosome: levels_allowed['Chromosome'] = [0, set([])]
if args.scaffold: levels_allowed['Scaffold'] = [0, set([])]
if args.contig: levels_allowed['Contig'] = [0, set([])]

for assembly_id, assembly in d.items():
    if args.taxid in assembly['lineage']:
        if assembly['assembly_level'] in levels_allowed:
            levels_allowed[assembly['assembly_level']][1].add(assembly['species_taxid_updated'])
            levels_allowed[assembly['assembly_level']][0] += 1

print('{:<20}{:<12}{:<12}'.format('ASSEMBLY_LEVEL', 'N_SPECIES', 'N_GENOMES'))
for level in levels_allowed:
    print('{:<20}{:<12,d}{:<12,d}'.format(level, len(levels_allowed[level][1]), levels_allowed[level][0]))

# Total species and genomes
s = set([])
g = 0
for level in levels_allowed:
    for taxid in levels_allowed[level][1]:
        s.add(taxid)
    g += levels_allowed[level][0]
print('{:<20}{:<12,d}{:<12,d}'.format('IN TOTAL', len(s), g))