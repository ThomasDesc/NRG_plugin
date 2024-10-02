import re
import os

def mif_to_pml(mif_file: str, output_dir: str = "./", output_tag: str = ""):
    """
    Convert a MIF file to a PML file for visualization in PyMOL.

    Parameters:
    mif_file (str): Path to the MIF file.
    output_dir (str): Directory to save the generated PML file.
    output_tag (str): Optional tag to append to the output file name.

    Returns:
    str: Path to the generated PML file.
    """
    if not mif_file:
        raise ValueError("Error: Missing MIF File")

    # Extract MIF filename and directory
    mif_name = ""
    mif_folder = ""
    if re.search(r"/([a-z0-9_-]+)\.mif$", mif_file, re.IGNORECASE):
        mif_name = re.findall(r"/([a-z0-9_-]+)\.mif$", mif_file, re.IGNORECASE)[0]
        mif_folder = os.path.dirname(mif_file)
    elif re.search(r"^([a-z0-9_-]+)\.mif$", mif_file, re.IGNORECASE):
        mif_name = re.findall(r"^([a-z0-9_-]+)\.mif$", mif_file, re.IGNORECASE)[0]
        mif_folder = "./"
    else:
        raise ValueError("Error: MIF filename must contain a-z, 0-9, _, or - characters.")

    # Initialize data structures
    ca = []
    pseudo = []
    protgrid = []
    probes = []  # Ensure probes is a list of lists to store each probe's information
    probesLab = []
    grid = []
    probesInt = []
    nbpb = 0
    grids = 0
    gride = 3
    pbColors = ["aquamarine", "brightorange", "blue", "red", "limegreen", "lightmagenta"]
    gridLab = ["200", "150", "100", "050"]
    gridColors = ["br2", "palegreen", "hydrogen", "cyan"]
    thinness = 0

    # Read the MIF file and parse its content
    try:
        with open(mif_file, 'r') as file:
            for line in file:
                # Extract different patterns and store them in corresponding data structures
                if re.search(r"^#protein_grid_distance\s+([.0-9-]+)\s+to\s+([.0-9-]+)", line):
                    match = re.search(r"^#protein_grid_distance\s+([.0-9-]+)\s+to\s+([.0-9-]+)", line)
                    thinness = abs(float(match.group(2)) - float(match.group(1)))
                elif re.search(r"^#ATOM", line):
                    match = re.search(r"#ATOM\s+([a-z0-9]+)\s+([0-9]+)\s+([a-z0-9]+)\s+([0-9]+)\s+([a-z0-9]{1})\s+([.0-9-]+)\s+([.0-9-]+)\s+([.0-9-]+)\s+([0-9]{1})\s+([0-9]{1})", line, re.IGNORECASE)
                    if match and match.group(10) == "1" and match.group(3) == "CA":
                        ca.append(f"{match.group(1)} {match.group(2)} {match.group(3)} {match.group(4)} {match.group(5)} {match.group(6)} {match.group(7)} {match.group(8)}")
                elif re.search(r"^#PSEUDO", line):
                    match = re.search(r"#PSEUDO\s+([a-z]+)\s+([.0-9-]+)\s+([.0-9-]+)\s+([.0-9-]+)", line, re.IGNORECASE)
                    if match:
                        pseudo.append(f"{match.group(1)} {match.group(2)} {match.group(3)} {match.group(4)}")
                elif re.search(r"^#probe\[([0-9]+)\]\s+([0-9a-z\.-]+)$", line, re.IGNORECASE):
                    match = re.search(r"^#probe\[([0-9]+)\]\s+([0-9a-z\.-]+)$", line, re.IGNORECASE)
                    probesLab.append(match.group(2))
                    nbpb += 1
                    # Create an empty list for each new probe
                    probes.append([])
                elif re.search(r"^#zip ([0-9]+)$", line, re.IGNORECASE):
                    match = re.search(r"^#zip ([0-9]+)$", line, re.IGNORECASE)
                    grids = gride = int(match.group(1)) if match.group(1) != "-1" else 3
                elif re.search(r"^#ff ([0-9a-z]+)$", line, re.IGNORECASE):
                    match = re.search(r"^#ff ([0-9a-z]+)$", line, re.IGNORECASE)
                    ff = match.group(1)
                elif re.search(r"^#PG\s+([.0-9-]+)\s+([.0-9-]+)\s+([.0-9-]+)$", line):
                    match = re.search(r"^#PG\s+([.0-9-]+)\s+([.0-9-]+)\s+([.0-9-]+)$", line)
                    protgrid.append(f"{match.group(1)} {match.group(2)} {match.group(3)}")
                elif not line.startswith("#") and line.strip():
                    # Store vertex potential interaction and grid presence
                    info = line.split()
                    for i in range(3, len(info) - 4, 3):
                        pbid = (i // 3) - 1
                        g0 = len(info) - 4
                        g1 = len(info) - 3
                        g2 = len(info) - 2
                        g3 = len(info) - 1
                        if int(info[i]) == 1 and len(probes) > pbid:
                            # Append interaction data to the corresponding probe
                            probes[pbid].append([info[0], info[1], info[2], info[g0], info[g1], info[g2], info[g3], info[i + 1], info[i + 2]])

                    for i in range(len(info) - 4, len(info)):
                        gid = i - (len(info) - 4)
                        gp = info[i]
                        bu = info[-1]
                        # Ensure the grid list has enough sublists to accommodate the gid index
                        while len(grid) <= gid:
                            grid.append([])  # Append empty sublists until the grid has enough elements

                        grid[gid].append([info[0], info[1], info[2], bu])

    except FileNotFoundError:
        raise FileNotFoundError(f"Error: Cannot open MIF file {mif_file}")

    # Generate the PML file
    output_filename = os.path.join(output_dir, f"{mif_name}{output_tag}.pml")
    with open(output_filename, 'w') as pml:
        pml.write("feedback disable,all,output\n")
        # Reading PDB and printing to PML
        try:
            with open(os.path.join(mif_folder, f"{mif_name}_cpy.pdb"), 'r') as pdb_file:
                pml.write('cmd.read_pdbstr("""\n')
                for line in pdb_file:
                    pml.write(line.rstrip() + "\\\n")
                pml.write(f'TER \\\n""","{mif_name}")\nset connect_mode,1\n')
        except FileNotFoundError:
            raise FileNotFoundError(f"Error: Cannot open {mif_name}_cpy.pdb in {mif_folder}")

        # Print pseudocenters
        if pseudo:
            pml.write('cmd.read_pdbstr("""\n')
            it = 0
            for p in pseudo:
                s = p.split()
                pml.write("HETATM{:5d}  N   {:3s} A0000    {:8.3f}{:8.3f}{:8.3f}  0.00 10.00           N\\\n".format(it, s[0].upper(), float(s[1]), float(s[2]), float(s[3])))
                it += 1
            pml.write(f'TER \\\n""","pseudocenters")\nhide nonbonded\nset connect_mode,1\n')

        # Print MIF points for each probe and each grid resolution
        for i in range(nbpb):
            if len(probes) > i:
                for g in range(grids, gride + 1):
                    if len(probes[i]) > 0:
                        probesInt.append([i, g])
                        pml.write('cmd.read_pdbstr("""\n')
                        for j in range(0, len(probes[i]), 9):
                            # Validate that probes[i][j+3+g] is an integer
                            if probes[i][j + 3 + g] and isinstance(probes[i][j + 3 + g], (int, float)) and int(probes[i][j + 3 + g]) == 1:
                                pml.write("HETATM{:5d}  N   {:3s} A0000    {:8.3f}{:8.3f}{:8.3f}  0.00 10.00           N\\\n".format(it, probesLab[i], float(probes[i][j][0]), float(probes[i][j][1]), float(probes[i][j][2])))
                                it = min(it + 1, 99999)
                        pml.write(f'TER \\\n""","{probesLab[i]}_{gridLab[g]}")\n')

        # Print grid points
        for i in range(grids, gride + 1):
            if len(grid) > i and len(grid[i]) > 0:
                pml.write('cmd.read_pdbstr("""\n')
                for j in range(0, len(grid[i]), 4):
                    pml.write("HETATM{:5d}  N   {:3s} A0000    {:8.3f}{:8.3f}{:8.3f}  0.00 10.00           N\\\n".format(it, gridLab[i], float(grid[i][j][0]), float(grid[i][j][1]), float(grid[i][j][2])))
                    it = min(it + 1, 99999)
                pml.write(f'TER \\\n""","{gridLab[i]}")\n')

        # Print protein grid points
        if len(protgrid) > 0:
            pml.write('cmd.read_pdbstr("""\n')
            for s in protgrid:
                s = s.split()
                pml.write("HETATM{:5d}  N   PGD A0000    {:8.3f}{:8.3f}{:8.3f}  0.00 10.00           N\\\n".format(0, float(s[0]), float(s[1]), float(s[2])))
            pml.write(f'TER \\\n""","{mif_name}_grid")\n')

        # Final PML configurations and displays
        pml.write("feedback enable,all,output\norient\nshow cartoon, {}\nremove (resn HOH)\n".format(mif_name))
        pml.write("show sticks, HET & {}\ncolor white,{}_grid\n".format(mif_name, mif_name))
        pml.write("show nonbonded,{}_grid\n".format(mif_name))
        for i in range(nbpb):
            if len(probes) > i:
                for g in range(grids, gride + 1):
                    if probesInt[i][1] != probesInt[i][1]:
                        pml.write("show spheres, {}_{}\nset sphere_scale,0.2,{}_{}\nrebuild\n".format(probesLab[i], gridLab[g], probesLab[i], gridLab[g]))
                        pml.write("color {},{}_{}\nhide nonbonded,{}_{}\n\n".format(pbColors[i], probesLab[i], gridLab[g], probesLab[i], gridLab[g]))

        for i in range(3):
            if len(grid) > i and len(grid[i]) > 0:
                pml.write("color {},{}\nshow nonbonded,{}\n".format(gridColors[i], gridLab[i], gridLab[i]))

    return output_filename
