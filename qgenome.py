import os
import argparse


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--indir', '--i', '-i', metavar="DIRECTORY", type=str, default="bacillus_gff",
                    help='input directory with query genus with genome subdirectories for GFF [DEFAULT: %(default)s]')
args = parser.parse_args()


PATH = args.indir

#### get the location for query protein 
def get_assembly_protein_location(genus, species, strain):
    ##### genus Path #######
    bacteria_gff_path = os.path.join(PATH, genus)
    
    if not os.path.exists(bacteria_gff_path):
        print('Input directory {} does not exist.'.format(bacteria_gff_path))
        exit()
    ### query genome #############
    assembly = strain
    gffassembly_path = os.path.join(bacteria_gff_path, str(species), assembly)
    gffassembly_name = '{}_genomic_gff.txt'.format(assembly)
    gff_assembly_name_path = os.path.join(gffassembly_path, gffassembly_name)

    # ### Define the dictionary 
    protein_location = {}
    fh = open(gff_assembly_name_path, "r")

    if assembly not in protein_location :# and assembly not in protein_position and assembly not in protein_accession:
        protein_location[assembly] = {}
        for lines in fh:
            sl = lines.rstrip().split("\t")
            uniq_id = str(sl[0])
            accession_id = str(sl[1])
            if accession_id not in protein_location[assembly]:
                protein_location[assembly][accession_id] = []
            protein_location[assembly][accession_id].append(uniq_id)

    fh.close()
    
    return(protein_location)

### get the accession id 

def protein_to_accession(genus, species, assembly, qid):
    input_path = os.path.join(PATH, genus)
    gffassembly_path = os.path.join(input_path, str(species), assembly) 
    #print(gffassembly_path)

    gffassembly_name = '{}_genomic_gff.txt'.format(assembly)
    gff_assembly_name_path = os.path.join(gffassembly_path, gffassembly_name)
    #print(gff_assembly_name_path)

    protein_accession = {}

    fh = open(gff_assembly_name_path, "r")
    for lines in fh:
        sl = lines.rstrip().split("\t")
        uniq_id = str(sl[0])
        accession_id = str(sl[1])
        

        if assembly not in protein_accession:
            protein_accession[assembly] = {}
        if uniq_id not in protein_accession[assembly]:
            protein_accession[assembly][uniq_id] = accession_id
    
    fh.close()
    
    return(protein_accession)

### get the protein seq postion for output 
def protein_seq_postion(genus, species, assembly):
    input_path = os.path.join(PATH, genus)
    gffassembly_path = os.path.join(input_path, str(species), assembly) 

    gffassembly_name = '{}_genomic_gff.txt'.format(assembly)
    gff_assembly_name_path = os.path.join(gffassembly_path, gffassembly_name)

    protein_position = {}

    fh = open(gff_assembly_name_path, "r")
    for lines in fh:
        sl = lines.rstrip().split("\t")
        uniq_id = str(sl[0])
        accession_id = str(sl[1])
        st,end, strand = int(sl[3]), int(sl[4]), sl[6].split("|")[6]
        
        if assembly not in protein_position:
            protein_position[assembly] = {}
        if accession_id not in protein_position[assembly]:
            protein_position[assembly][accession_id] = {}
        if uniq_id not in protein_position[assembly][accession_id]:
            protein_position[assembly][accession_id][uniq_id] = st,end, strand
    
    fh.close()
    

    return(protein_position)