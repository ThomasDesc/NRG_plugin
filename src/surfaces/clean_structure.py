# Developed by Natália Teruel
# Najmanovich Research Group
# Cite this work as Surfaces: A software to quantify and visualize interactions within and between proteins and ligands - Teruel, N. F. B., Borges, V. M., & Najmanovich, R. (2023)

#Imports

def get_atoms(line): #take strings after the ' ' and before the ':'
    list_atoms = []
    L = line[5:-2].split(",")
    for item in L:
        if len(item) > 0:
            item = item.strip()
            id = item.index(":")
            list_atoms.append(item[:id])
    return (list_atoms)

def read_atom(line):
    atnum = int(line[6:11])
    attype = line[11:16].strip()
    resnum = line[23:30].strip()
    res = line[17:20].strip()
    chain = line[21:22].strip()
    return (atnum,attype,resnum,res,chain)

def check_line(line,res,atoms): #see if the atom of that line is described in the def file
    atnum,attype,resnum,residue,chain = read_atom(line)
    if residue in res:
        id = res.index(residue)
        if attype in atoms[id]:
            return (True)
        else:
            print ("WARNING: ATOM " + attype + " NOT DEFINED FOR " + residue)
            return (False)
    else:
        print ("WARNING: " + residue + " NOT DEFINED")
        return (False)


def get_args():
    import argparse
    parser = argparse.ArgumentParser(description="the arguments.", add_help=False)
    parser.add_argument("-f", "--pdb_file", action="store")
    parser.add_argument("-def", "--atomtypes_definition", action="store")
    args = parser.parse_args()
    pdb_file_path = args.pdb_file
    atomtypes_definition_path = args.atomtypes_definition
    clean_pdb_file = open("clean_" + pdb_file_path, "w")
    main(pdb_file_path, atomtypes_definition_path, clean_pdb_file)


def main(pdb_file_path, atomtypes_definition_path, clean_pdb_file):

    res = []
    atoms = []
    def_file = open(atomtypes_definition_path, "r")
    Lines = def_file.readlines()
    for line in Lines:
        res.append(line[:3])
        atoms.append(get_atoms(line))
    def_file.close()
    
    #print (res)
    #print (atoms)
    
    pdb_file = open(pdb_file_path, "r")
    Lines = pdb_file.readlines()
    for line in Lines:
        if line[:4] == 'ATOM' or line[:4] == 'HETA':
            if check_line(line,res,atoms):
                clean_pdb_file.write(line)
        if line[:4] == 'CONE':
            clean_pdb_file.write(line)
    clean_pdb_file.write('END')
    pdb_file.close()
    clean_pdb_file.close()
    
    return


if __name__ == '__main__':
    get_args()
