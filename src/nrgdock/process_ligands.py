import os
import numpy as np
import sys
install_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(install_dir)
from src.nrgdock.main_processed_target import get_params_dict
import concurrent.futures
from itertools import repeat
import argparse
import json
from pathlib import Path
import pickle


def load_rad_dict(filepath):
    with open(filepath, 'r') as file:
        loaded_atom_data = json.load(file)
    return loaded_atom_data


def get_radius_number(letter_type, rad_dict):
    letter_type = letter_type.upper().replace(' ', '')
    try:
        atm_info = [rad_dict[letter_type]['type_number'], rad_dict[letter_type]['radius']]
    except KeyError:
        atm_info = [39, 2.00]
    atm_type_num = atm_info[0]
    atm_rad = atm_info[1]
    return atm_type_num, atm_rad


def load_atoms_mol2(filename, rad_dict, save_path, ligand_type='ligand'):
    coord_start = 0
    max_atoms = 0
    n_atoms = 0
    n_molecules = 0
    n_unique_molecules = 0
    same_molec_counter = 1
    molecule_name_list = []
    with open(filename) as f:
        lines = f.readlines()

    for counter, line in enumerate(lines):
        if line.startswith('@<TRIPOS>MOLECULE'):
            n_molecules += 1
            molecule_name = lines[counter+1][0:-1]
            if molecule_name != "\n":
                molec_suffix = "_0"
                if n_molecules > 1 and molecule_name_list[n_molecules-2].split("_")[0] == molecule_name:
                    molec_suffix = f"_{same_molec_counter}"
                    same_molec_counter += 1
                else:
                    same_molec_counter = 1
                molecule_name_list.append(lines[counter+1][0:-1] + molec_suffix)
                if molec_suffix =="_0":
                    n_unique_molecules += 1
            else:
                exit("Error when reading molecule name")
        elif line.startswith("@<TRIPOS>ATOM"):
            coord_start = 1
            if max_atoms < n_atoms:
                max_atoms = n_atoms
            n_atoms = 0
        if coord_start == 1:
            if line.startswith("@<TRIPOS>ATOM") is False:
                if line[0] == "@":
                    coord_start = 0
                elif line.split()[1][0] != "H":
                    n_atoms += 1

    if max_atoms < n_atoms:
        max_atoms = n_atoms

    n_atom_array = np.zeros(n_molecules, dtype=np.int32)
    atoms_xyz = np.full((n_molecules, max_atoms, 3), 9999, dtype=np.float32)
    atoms_type = np.full((n_molecules, max_atoms), -1, dtype=np.int32)
    # atoms_name = np.zeros((n_molecules, max_atoms), dtype=object)

    molecule_counter = -1
    atom_counter = 0
    coord_start = 0
    atom_name_list = []
    temp_atom_name_list = []
    atoms_name_count = {}
    for counter, line in enumerate(lines):
        if line.startswith('@<TRIPOS>MOLECULE'):
            atoms_name_count = {}
            molecule_counter += 1
            if molecule_counter > 0:
                n_atom_array[molecule_counter-1] = atom_counter
                atom_counter = 0
                coord_start = 0
                atom_name_list.append(temp_atom_name_list)
                temp_atom_name_list = []
        if line.startswith('@<TRIPOS>ATOM') and line[0] != "\n":
            coord_start = 1
        if coord_start == 1 and line.startswith('@<TRIPOS>ATOM') is False:
            if line[0] != '@':
                line = line.split()
                if line[5][0] != 'H':
                    atoms_xyz[molecule_counter][atom_counter] = np.array([float(line[2]),
                                                                          float(line[3]),
                                                                          float(line[4])])
                    atom_type = line[5]
                    atoms_type_temp, _ = get_radius_number(atom_type, rad_dict)
                    atoms_type[molecule_counter][atom_counter] = atoms_type_temp
                    atm_name = atom_type.split(".")[0]
                    if atm_name in atoms_name_count:
                        atoms_name_count[atm_name] += 1
                    else:
                        atoms_name_count[atm_name] = 1
                    # atoms_name[molecule_counter][atom_counter] = f"{atm_name}{atoms_name_count[atm_name]}"
                    temp_atom_name_list.append(f"{atm_name}{atoms_name_count[atm_name]}")
                    atom_counter += 1
            else:
                coord_start = 0
    n_atom_array[-1] = atom_counter
    atom_name_list.append(temp_atom_name_list)
    np.save(os.path.join(save_path, f"{ligand_type}_atom_xyz"), atoms_xyz)
    np.save(os.path.join(save_path, f"{ligand_type}_atom_type"), atoms_type)
    # np.save(os.path.join(save_path, f"{ligand_type}_atom_name"), atoms_name)
    # np.save(os.path.join(save_path, f"{ligand_type}_molecule_name"), molecule_names)
    with open(os.path.join(save_path, f"{ligand_type}_molecule_name.pkl"), 'wb') as f:
        pickle.dump(molecule_name_list, f)
    with open(os.path.join(save_path, f"{ligand_type}_atom_name.pkl"), 'wb') as f:
        pickle.dump(atom_name_list, f)
    np.save(os.path.join(save_path, f"{ligand_type}_atoms_num_per_ligand"), n_atom_array)
    return n_unique_molecules


def get_suffix(conf_num):
    suffix = ""
    if conf_num != 0:
        suffix = f"_{conf_num}_conf"
    return suffix


def preprocess_ligands_one_target(rad_dict, conf_num, target_path, command, ligand_path, custom_output_path,
                                  ligand_type='ligand'):
    print(f"Preprocessing ligands for {os.path.basename(target_path)}")

    if command == "single_file":
        ligand_list = [ligand_path]
    elif command == 'folder':
        filenames = next(os.walk(target_path), (None, None, []))[2]
        ligand_list = []
        for filename in filenames:
            if filename.endswith("ligands.mol2") or filename.endswith("final.mol2"):
                ligand_list.append(os.path.join(target_path, filename))
    else:
        exit("Error: user specified something other than single_file or folder")

    for ligand in ligand_list:
        suffix = get_suffix(conf_num)
        np_array_savepath = os.path.join(target_path, f"preprocessed_ligands{suffix}")
        if custom_output_path is not None:
            np_array_savepath = custom_output_path
        print("saving to: ", np_array_savepath)
        if not os.path.exists(np_array_savepath):
            os.makedirs(np_array_savepath)
        n_unique_molecules = load_atoms_mol2(ligand, rad_dict, np_array_savepath, ligand_type)
        if command == 'single_file':
            return n_unique_molecules, np_array_savepath


def get_args():
    target = None
    target_path = None
    ligand_path = None
    custom_output_path = None

    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest='command')

    folder = subparser.add_parser('folder')
    folder.add_argument("-p", '--target_path', required=True, default=None, type=str,
                        help='Path to folder containing target folders')
    folder.add_argument("-t", '--target', default=None, type=str, help='Name of folder of specific target')

    single_file = subparser.add_parser('single_file')
    single_file.add_argument("-l", '--ligand_path', required=True, type=str, help='Path to ligand file')
    single_file.add_argument("-o", '--output_path', type=str, help='Custom output path')

    args = parser.parse_args()
    command = args.command

    if command == "folder":
        target_path = args.target_path
        target = args.target

    if command == "single_file":
        ligand_path = args.ligand_path
        target_path = os.path.dirname(ligand_path)
        custom_output_path = args.output_path

    main(command, target, target_path, ligand_path, custom_output_path)


def main(command, target, target_path, ligand_path, custom_output_path):
    root_software_path = Path(__file__).resolve().parents[1]
    os.chdir(root_software_path)
    config_file = "./deps/config.txt"
    params_dict = get_params_dict(config_file)
    rad_dict = load_rad_dict("./deps/atom_type_radius.json")
    conf_num = params_dict["CONFORMER_NUMBER"]
    target_path = os.path.normpath(target_path)

    if command == "single_file":
        preprocess_ligands_one_target(rad_dict, conf_num, target_path, command, ligand_path, custom_output_path)
    elif command == "folder":
        ligand_path = []
        if target is not None:
            final_target_path = [os.path.join(target_path, target)]
        else:
            final_target_path = sorted(os.listdir(target_path))
            for counter, target in enumerate(final_target_path):
                final_target_path[counter] = os.path.join(target_path, target)
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(preprocess_ligands_one_target, repeat(rad_dict), repeat(conf_num), final_target_path, repeat(command), repeat(ligand_path), repeat(custom_output_path))
        # preprocess_ligands_one_target(rad_dict, conf_num, final_target_path[0], command, ligand_path, custom_output_path, ligand_type='active')


if __name__ == "__main__":
    get_args()
