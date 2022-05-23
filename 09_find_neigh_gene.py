import os
import argparse
import qgenome
import subprocess



parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('--indir1', '--i1', '-i1', metavar="DIRECTORY", type=str, default="db_bacteria_gff",
                    help='input directory with db genus with genome subdirectories for GFF [DEFAULT: %(default)s]')

parser.add_argument('--taxid', '-taxid', dest='taxid', type=int,  default=1386,
                    help='Taxonomy ID of query genus group of organisms [DEFAULT: %(default)s]')

parser.add_argument('--speciestaxid', '-speciestaxid', dest='speciestaxid', type=int,  default=1455,
                    help='Taxonomy ID of query species group of organisms [DEFAULT: %(default)s]')

parser.add_argument('--strain', '-strain', dest='strain', type=str, default="GCF_001630115.1_ASM163011v1",
                    help='query genome ID  [DEFAULT: %(default)s]')

parser.add_argument('--indir2', '--i2', '-i2', metavar="DIRECTORY", type=str, default="reciprocal_homo",
                    help='input directory with genus subdirectories for homologs for each species [DEFAULT: %(default)s]')
                  
parser.add_argument('--outdir', '--o', '-o', metavar="DIRECTORY", type=str, default="01-Genus_Neighbouring_gene_homolog_TRG_trial_sp",
                    help='output directory to save results [DEFAULT: %(default)s]')

args = parser.parse_args()

def Neigh_homo(left_index, left_gene, DB_STRAIN_TO_ACCESSION, db_genome):
    flag = None
    if left_gene in HOMOLOG[db_genome]:
        homo_protein_id = HOMOLOG[db_genome][left_gene].split("|")[1]
        if db_genome in DB_STRAIN_TO_ACCESSION:
            for accession in DB_STRAIN_TO_ACCESSION[db_genome]:
                if homo_protein_id in DB_STRAIN_TO_ACCESSION[db_genome][accession]:
                    flag = "{}\t{}\t{}".format(db_genome, accession, "\t".join(map(str,DB_STRAIN_TO_ACCESSION[db_genome][accession][homo_protein_id])))
    return(flag)

DIC_Strains_PATH = {}
### GEt the Path for DB Genomes 
NCBI_GFF_PATH = args.indir1 
for ncbi_genus_taxid in os.listdir(NCBI_GFF_PATH):
    if ncbi_genus_taxid == str(args.taxid):
        continue
    NC_GENUS_GFF_PATH = os.path.join(NCBI_GFF_PATH,ncbi_genus_taxid)
    #print(NC_GENUS_GFF_PATH)

    for ncbi_species_taxid in os.listdir(NC_GENUS_GFF_PATH):
        NC_species_path = os.path.join(NC_GENUS_GFF_PATH, ncbi_species_taxid)
        
        for ncbi_strain in os.listdir(NC_species_path):
            NC_strain_path = os.path.join(NC_species_path, ncbi_strain)
            
            db_gffassembly_name = '{}_genomic_gff.txt'.format(ncbi_strain)
            db_gff_assembly_name_path = os.path.join(NC_strain_path, db_gffassembly_name)
            
            if ncbi_strain not in DIC_Strains_PATH:
                DIC_Strains_PATH[ncbi_strain] = db_gff_assembly_name_path


TRG_DATA = {}
genus_trg_data_file = open("genus_trg.out", "r") # text file for trg with their genome strain ID 
data = genus_trg_data_file.readlines()[1:]
for line in data:
    sl= line.rstrip().split("\t")
    trg_id = sl[0]
    trg_original_id = sl[1]
    trg_genome = sl[2]
    trg_species = trg_id.split("_")[1]

    combine_key = trg_species, trg_genome

    if combine_key not in TRG_DATA:
        TRG_DATA[combine_key] = set()
    TRG_DATA[combine_key].add(trg_id)
    
genus_trg_data_file.close()


### out dir ###
OUTPATH = args.outdir
os.makedirs(OUTPATH, exist_ok=True)

query_trg_species, query_trg_genome = str(args.speciestaxid), str(args.strain)


## open the reciprocal homolog file for query species 
reciprocal_homo_path = args.indir2
rbh_file_name_path = os.path.join(reciprocal_homo_path, "{0}.rbh.out".format(str(args.speciestaxid)))

HOMOLOG = {}
LIST_OF_HOMOLOG_GENOME = set()
awk_query_genome_info = "awk '{{if ($3 ==\"{0}\") print $0}}' {1}".format(query_trg_genome, rbh_file_name_path)
awk_output = subprocess.check_output(['bash', '-c', awk_query_genome_info], universal_newlines=True)


grep_lines = awk_output.split("\n")
grep_lines = filter(None,grep_lines)
for line in grep_lines:
    spl = line.rstrip().split("\t")
    queryid = spl[0]
    subject_species = spl[4]
    subject_genome = spl[5]
    subject_id = spl[6]

    LIST_OF_HOMOLOG_GENOME.add(subject_genome)
    

    if subject_genome not in HOMOLOG:
        HOMOLOG[subject_genome] = {}
    if queryid not in HOMOLOG[subject_genome]:
        HOMOLOG[subject_genome][queryid] = subject_id
    
trg_idenitfierlist = TRG_DATA[query_trg_species, query_trg_genome] ### get the list of TRG ids in the query genome 
print(trg_idenitfierlist)
out_species_dir = os.path.join(OUTPATH, query_trg_species)
os.makedirs(out_species_dir, exist_ok=True)

oh = open(os.path.join(out_species_dir, "{}_{}.out".format(query_trg_species, query_trg_genome)), "w")


for homo_genome in LIST_OF_HOMOLOG_GENOME:
    
    sub_strain_path = DIC_Strains_PATH.get(homo_genome, "")
    
    ### open the db gff file  and make dictionary 
    DB_STRAIN_TO_ACCESSION = {}
    if os.path.exists(sub_strain_path):
        db_gff_strain_file = open(sub_strain_path)
        for db_gff_line in db_gff_strain_file:
            sl_db_gff_line = db_gff_line.rstrip().split('\t')
        
            db_accession = str(sl_db_gff_line[1])
            db_protein_id = str(sl_db_gff_line[0])
            db_protein_st, db_protein_end, db_protein_strand = int(sl_db_gff_line[3]), int(sl_db_gff_line[4]), sl_db_gff_line[6].split("|")[6]
            
            if homo_genome not in DB_STRAIN_TO_ACCESSION:
                DB_STRAIN_TO_ACCESSION[homo_genome] = {}
            if db_accession not in DB_STRAIN_TO_ACCESSION[homo_genome]:
                DB_STRAIN_TO_ACCESSION[homo_genome][db_accession] = {}
            if db_protein_id not in DB_STRAIN_TO_ACCESSION[homo_genome][db_accession]:
                DB_STRAIN_TO_ACCESSION[homo_genome][db_accession][db_protein_id] = db_protein_st, db_protein_end, db_protein_strand

        db_gff_strain_file.close()
    
    for trg in trg_idenitfierlist:

        ## Return scaffold id for trg protein 
        query_protein_to_scaffold = qgenome.protein_to_accession(str(args.taxid), query_trg_species, query_trg_genome, trg)
        return_protein_to_scaffold_id = query_protein_to_scaffold[query_trg_genome][trg]
        #################################

        ### return location for trg protein in the query genome 
        protein_location_dic = qgenome.get_assembly_protein_location(str(args.taxid), query_trg_species, query_trg_genome)
        index_for_trg_protein = protein_location_dic[query_trg_genome][return_protein_to_scaffold_id].index(trg)
        #################################################
        ## return cds gene cordinates for trg protein  in the query genome 
        cds_coordinates_dic = qgenome.protein_seq_postion(str(args.taxid), query_trg_species, query_trg_genome)
        trg_protein_postion = cds_coordinates_dic[query_trg_genome][return_protein_to_scaffold_id][trg]

        ###### FOR LEFT HOMO #######

        if int(index_for_trg_protein) >0:
            left_gene_index = int(index_for_trg_protein)-1 
            left_gene = protein_location_dic[query_trg_genome][return_protein_to_scaffold_id][left_gene_index]
            left_gene_cordinates = cds_coordinates_dic[query_trg_genome][return_protein_to_scaffold_id][left_gene]
            
            leftAccnPos = "left", query_trg_genome,return_protein_to_scaffold_id, trg, "\t".join(map(str,trg_protein_postion)), left_gene_index,left_gene,"\t".join(map(str,left_gene_cordinates))
            
            return_hit_for_left = None
            while left_gene_index >=0 and return_hit_for_left == None:
                return_hit_for_left = Neigh_homo(left_gene_index, left_gene, DB_STRAIN_TO_ACCESSION, homo_genome)
                leftAccnPos = "left",query_trg_genome,return_protein_to_scaffold_id,trg, "\t".join(map(str,trg_protein_postion)), left_gene_index,left_gene,"\t".join(map(str,left_gene_cordinates))
                left_gene_index = left_gene_index -1
                left_gene = protein_location_dic[query_trg_genome][return_protein_to_scaffold_id][left_gene_index]
                left_gene_cordinates = cds_coordinates_dic[query_trg_genome][return_protein_to_scaffold_id][left_gene]

            oh.write("{}\t{}\n".format("\t".join(map(str,leftAccnPos)), return_hit_for_left))

        ###### FOR RIGHT HOMO ########
        right_gene = None
        
        if int(index_for_trg_protein) < len(protein_location_dic[query_trg_genome][return_protein_to_scaffold_id]) and int(len(protein_location_dic[query_trg_genome][return_protein_to_scaffold_id]))!=1 :
            right_gene_index = int(index_for_trg_protein)+1
            try:
                right_gene = protein_location_dic[query_trg_genome][return_protein_to_scaffold_id][right_gene_index]
                right_gene_cordinates = cds_coordinates_dic[query_trg_genome][return_protein_to_scaffold_id][right_gene]
                rightAccnPos= "right",query_trg_genome,return_protein_to_scaffold_id, trg, "\t".join(map(str,trg_protein_postion)), right_gene_index, right_gene, "\t".join(map(str,right_gene_cordinates))
                
            except:
                pass
        
            if right_gene!=None:
                return_hit_for_right = None
                while right_gene_index <= len(protein_location_dic[query_trg_genome][return_protein_to_scaffold_id]) and return_hit_for_right ==None:
                    return_hit_for_right = Neigh_homo(right_gene_index, right_gene, DB_STRAIN_TO_ACCESSION, homo_genome)
                    rightAccnPos= "right",query_trg_genome,return_protein_to_scaffold_id, trg, "\t".join(map(str,trg_protein_postion)), right_gene_index, right_gene, "\t".join(map(str,right_gene_cordinates))
                    right_gene_index = right_gene_index + 1
                    
                    if right_gene_index > int(len(protein_location_dic[query_trg_genome][return_protein_to_scaffold_id])):
                        break
                    try:
                        right_gene = protein_location_dic[query_trg_genome][return_protein_to_scaffold_id][right_gene_index]
                        right_gene_cordinates = cds_coordinates_dic[query_trg_genome][return_protein_to_scaffold_id][right_gene]
                    except:
                        break

                oh.write("{}\t{}\n".format("\t".join(map(str,rightAccnPos)), return_hit_for_right))

oh.close()
    
