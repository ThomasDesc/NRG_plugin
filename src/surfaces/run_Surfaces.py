import os
import matplotlib.pyplot as plt
from src.surfaces.ligand_atomtypes import add_pdb
from src.surfaces.clean_structure import main as clean_structure
from src.surfaces.surface_cont_lig import main as surface_cont_lig
from src.surfaces.pymol_image_surfaces_lig import generate_session
from pymol import cmd
from general_functions import get_residue_info
import pandas as pd


def process_result_flexaid(flexaid_result_file, output):
    with open(flexaid_result_file, 'r') as t1:
        with open(output, 'w') as t2:
            text = t1.readlines()
            for line in text:
                if 'REMARK' not in line:
                    if 'LIG  9999' in line:
                        a_name = str(int(line[12:17].split()[0])) + line[9:11] + ' ' * (
                                    5 - len(line[12:17].split()[0] + line[9:11]))
                        new_line = line[:12] + a_name + line[17:21] + 'L' + line[22:]
                        t2.write(new_line)
                    else:
                        t2.write(line)


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
    lig_pdb_file.close()
    return


def flex_res(target_file):
    with open(target_file, "r") as f:
        texto=f.readlines()
        for line in texto:
            if 'LIG  9999' in line:
                return 1
    return 0


def get_chains_from_object(object_name):
    chains = cmd.get_chains(object_name)
    str_chain=''
    for chain in chains:
        str_chain=str_chain+str(chain)
    return str_chain


def load_surfaces(form, temp_path, main_folder_path, binary_folder_path, binary_suffix):
    vcon_binary_path = os.path.join(binary_folder_path, f'vcon{binary_suffix}')
    target = form.surface_select_result.currentText(),
    lig = form.surface_select_lig.currentText(),
    target_2 = form.surface_select_result_2.currentText(),
    lig_2 = form.surface_select_lig_2.currentText(),

    if lig!='None':
        if target_2 == 'None':
            target_file = temp_path + '/{}.pdb'.format(target)
            cmd.save(target_file, target)
            ligands=get_residue_info(lig)
            for ligand in ligands:
                target_chain=get_chains_from_object(target)
                interac_dic=run_surfaces_lig(target_file, target_chain, ligand[0], temp_path, main_folder_path, vcon_binary_path)
            interac_dic = dict(sorted(interac_dic.items(), key=lambda item: abs(item[1]),reverse=True))
            plot_interactive_table(interac_dic.keys(),interac_dic.values())
        else:
            cf_dic={}
            target_file = temp_path + '/{}.pdb'.format(target)
            cmd.save(target_file, target)
            ligands = get_residue_info(lig)
            for ligand in ligands:
                target_chain = get_chains_from_object(target)
                interac_dic = run_surfaces_lig(target_file, target_chain, ligand[0], temp_path, main_folder_path,
                                               vcon_path)
                interac_dic = dict(sorted(interac_dic.items(), key=lambda item: abs(item[1]), reverse=True))
                cf_ref=sum(interac_dic.values())
            #plot_interactive_table(interac_dic.keys(), interac_dic.values())
            for state in range(cmd.count_states(target_2)):
                output_file = temp_path + '/{}_{}.pdb'.format(target_2, state+1)
                cmd.save(output_file, target_2, state=state + 1)
                target_chain = get_chains_from_object(target_2)
                ligands = get_residue_info(lig)
                for ligand in ligands:
                    target_chain = get_chains_from_object(target)
                    interac_dic = run_surfaces_lig(output_file, target_chain, ligand[0], temp_path, main_folder_path,
                                               vcon_path)
                    cf_dic[os.path.basename(output_file) + '_' + ligand[0]] = sum(interac_dic.values())-cf_ref
            plt.bar(cf_dic.keys(), cf_dic.values())
            plt.xticks(rotation=90)
            plt.title('CF difference per state')
            plt.ylabel('deltaCF')
            plt.xlabel('State-Ligand')
            plt.tight_layout()
            plt.show()


def run_surfaces_lig(target_file,target_chain,lig, temp_path, main_folder_path, vcon_path):
    surfaces_output_path=temp_path+'/Surfaces'
    def_file = os.path.join(main_folder_path, "surfaces_defs", 'AMINO_FlexAID.def')
    flexaid_dat_path = os.path.join(main_folder_path, "surfaces_defs", 'FlexAID.dat')
    color_rgb_path = os.path.join(main_folder_path, "surfaces_defs", 'color_rgb.txt')
    open_def_file = open(def_file, "r")
    if flex_res(target_file):
        process_result_flexaid(target_file, target_file)
    target_file=target_file
    ligand_file_name = os.path.join(os.path.dirname(target_file), lig)
    create_ligand_file(target_file, ligand_file_name)
    custom_def_path = os.path.join(surfaces_output_path, f'custom_{os.path.basename(def_file)}')
    custom_def_file = open(custom_def_path, 'w')
    add_pdb(ligand_file_name + '.pdb', open_def_file, custom_def_file)
    open_def_file.close()
    custom_def_file.close()
    cleaned_file_path = os.path.join(os.path.dirname(target_file), f"cleaned_{os.path.basename(target_file)}")
    clean_pdb_file = open(cleaned_file_path, "w")
    clean_structure(target_file, custom_def_path, clean_pdb_file)
    clean_pdb_file.close()
    os.rename(cleaned_file_path,target_file[:-4]+'_.pdb')
    cleaned_file_path=target_file[:-4]+'_.pdb'
    vcon_out_file = os.path.join(os.path.dirname(surfaces_output_path), 'vcon_file.txt')
    csv_path = os.path.join(surfaces_output_path, 'csv_output.csv')
    list_file_path, image_file_path = surface_cont_lig(cleaned_file_path,target_chain, lig, csv_path, custom_def_path, flexaid_dat_path, vcon_path, vcon_out_file)
    generate_session(cleaned_file_path, image_file_path, list_file_path, color_rgb_path)
    interact_dic={}
    print(list_file_path)
    with open(list_file_path,'r') as f:
        texto=f.readlines()
        for line in texto:
            if line.split(',')[0][3:-1] in interact_dic.keys():
                interact_dic[line.split(',')[0][3:-1]]+=float(line.split(',')[-1][:-1])
            else:
                interact_dic[line.split(',')[0]] = float(line.split(',')[-1][:-1])
    return interact_dic


def plot_interactive_table(labels, values):
    # Create a DataFrame
    df = pd.DataFrame({'Res': labels, 'Value': values})

    # Create a figure and axis object
    fig, ax = plt.subplots()  # Adjust figure size based on the number of rows

    # Hide axes
    ax.axis('tight')
    ax.axis('off')

    # Create the table
    table = ax.table(cellText=df.values, colLabels=df.columns, cellLoc='center', loc='center')

    # Customize the table appearance (optional)
    table.auto_set_font_size(False)
    #table.set_fontsize(12)
    table.auto_set_column_width(col=[0, 1])
    plt.show()


