import os
import shutil
from pymol import cmd
from ligand_atomtypes import add_pdb
from clean_structure import main as clean_structure
from surface_cont_lig import main as surface_cont_lig
from pymol_image_surfaces_lig import generate_session


def process_result_flexaid(flexaid_result_file, output):
    with open(flexaid_result_file, 'r') as t1:
        with open(output, 'w') as t2:
            text = t1.readlines()
            for line in text:
                if 'REMARK' not in line:
                    if 'LIG  9999' in line:
                        a_name = line[12:17].split()[0] + line[9:11] + ' ' * (
                                    5 - len(line[12:17].split()[0] + line[9:11]))
                        new_line = line[:12] + a_name + line[17:21] + 'L' + line[22:]
                        t2.write(new_line)
                    else:
                        t2.write(line)


def load_surfaces_result(form, surfaces_output_path):
    cmd.delete('all')
    result_path = form.surfaces_load_result_text.text()
    result_base_name = os.path.basename(result_path)
    dst = os.path.join(surfaces_output_path, result_base_name)
    print(result_base_name)
    print(dst)
    shutil.copy(result_path, dst)
    form.simulate_folder_path.setText(os.path.dirname(dst))
    cmd.load(dst)
    form.surface_select_result.addItems([os.path.splitext(os.path.basename(dst))[0]])


def create_ligand_file(pdb_file_name, lig_path):
    with open(pdb_file_name, "r") as f:
        lines = f.readlines()
    lig_pdb_file = open(lig_path+".pdb", "w")
    lig_name = os.path.basename(lig_path)
    for line in lines:
        if line[:4] == 'ATOM' or line[:4] == 'HETA':
            res = line[17:20].strip()
            if res == lig_name:
                lig_pdb_file.write(line)
        if line[:4] == 'CONE':
            lig_pdb_file.write(line)
    lig_pdb_file.write('END')
    lig_pdb_file.close()
    return


def retrieve_flexaid_result(flexaid_simulation_folder):
    cmd.delete('all')
    files = os.listdir(flexaid_simulation_folder)
    for file in files:
        if file.endswith('.pdb') and not file.endswith('INI.pdb'):
            file_path = os.path.join(flexaid_simulation_folder, file)
            if os.path.isfile(file_path):
                cmd.load(file_path)


def run_run_surfaces(selected_result, surfaces_output_path, flexaid_simulation_folder, main_folder_path, vcon_path):
    def_file = os.path.join(main_folder_path, "surfaces_defs", 'AMINO_FlexAID.def')
    flexaid_dat_path = os.path.join(main_folder_path, "surfaces_defs", 'FlexAID.dat')
    color_rgb_path = os.path.join(main_folder_path, "surfaces_defs", 'color_rgb.txt')
    open_def_file = open(def_file, "r")
    flexaid_result_file = os.path.join(flexaid_simulation_folder, selected_result + '.pdb')
    processed_result_path = os.path.join(surfaces_output_path, os.path.basename(flexaid_result_file)[:-4] + '_processed.pdb')
    process_result_flexaid(flexaid_result_file, processed_result_path)
    ligand_file_name = os.path.join(os.path.dirname(processed_result_path), 'LIG')
    create_ligand_file(processed_result_path, ligand_file_name)
    custom_def_path = os.path.join(surfaces_output_path, f'custom_{os.path.basename(def_file)}')
    custom_def_file = open(custom_def_path, 'w')
    add_pdb(ligand_file_name + '.pdb', open_def_file, custom_def_file)
    open_def_file.close()
    custom_def_file.close()
    cleaned_file_path = os.path.join(os.path.dirname(processed_result_path), f"cleaned_{os.path.basename(processed_result_path)}")
    clean_pdb_file = open(cleaned_file_path, "w")
    clean_structure(processed_result_path, custom_def_path, clean_pdb_file)
    clean_pdb_file.close()
    vcon_out_file = os.path.join(os.path.dirname(surfaces_output_path), 'vcon_file.txt')
    csv_path = os.path.join(surfaces_output_path, 'csv_output.csv')
    list_file_path, image_file_path = surface_cont_lig(cleaned_file_path, 'ABC', 'LIG', csv_path, custom_def_path, flexaid_dat_path, vcon_path, vcon_out_file)
    generate_session(cleaned_file_path, image_file_path, list_file_path, color_rgb_path)