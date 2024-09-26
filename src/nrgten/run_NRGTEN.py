from modeller.automodel import *
from nrgten.encom import ENCoM
from nrgten.atom import Atom
from pymol import cmd
from src.surfaces.ligand_atomtypes import add_pdb
from src.surfaces.run_Surfaces import create_ligand_file, flex_res, process_result_flexaid,compare_residues
from src.nrgten.model_ensemble import model_states
import matplotlib.pyplot as plt
import os
import numpy as np


def remove_selection_and_save(object_name, selection, output_file):
    cmd.create("object_without_selection", f"{object_name} and not ({selection})")
    cmd.save(output_file, "object_without_selection")
    cmd.delete("object_without_selection")


def flex_aid_matrix(main_folder_path):
    matrix = np.array([np.zeros(41)] * 41)
    flexaid_dat_path =  os.path.join(main_folder_path, "deps","surfaces", "FlexAID.dat")
    with open(flexaid_dat_path, "r") as t1:
        texto = t1.readlines()
        for line in texto:
            matrix[int(line.split("=")[0].split('-')[0])][int(line.split("=")[0].split('-')[1])] = -float(
                line.split("=")[1])
    return matrix


def encom_model(target_file, main_folder_path, temp_path):
    if flex_res(target_file):
        process_result_flexaid(target_file, target_file)
    list_het = find_het(target_file, temp_path, main_folder_path)
    atype_lists = [os.path.join(main_folder_path, "deps", "nrgten", 'amino_acids.atomtypes')]
    mass_lists = [os.path.join(main_folder_path, "deps", "nrgten", 'amino_acids.masses')]
    for lig in list_het:
        atype_lists.append(os.path.join(temp_path, 'NRGTEN', f'{lig}.atomtypes'))
        mass_lists.append(os.path.join(temp_path,'NRGTEN', f'{lig}.masses'))

    matrix = flex_aid_matrix(main_folder_path)
    model = ENCoM(target_file, interact_mat=matrix, atypes_list=atype_lists, massdef_list=mass_lists)
    return model


def standardize_to_minus1_plus1(data):
    max_abs_value = max(abs(x) for x in data)
    standardized_data = [x / max_abs_value for x in data]
    return standardized_data


def find_het(target_file, temp_path, main_folder_path):
    het_dic = {}
    with open(target_file, "r") as t1:
        texto = t1.readlines()
        for line in texto:
            if 'HETATM' in line:
                het_dic[line[17:20]] = '1'
    for lig in list(het_dic):
        def_file = os.path.join(main_folder_path, "deps","surfaces", 'AMINO_FlexAID.def')
        flexaid_dat_path = os.path.join(main_folder_path, "deps","surfaces", 'FlexAID.dat')
        open_def_file = open(def_file, "r")
        ligand_file_name = os.path.join(os.path.dirname(target_file), lig)
        create_ligand_file(target_file, ligand_file_name)
        custom_def_path = os.path.join(temp_path,'NRGTEN', f'custom_{os.path.basename(def_file)}')
        custom_def_file = open(custom_def_path, 'w')
        add_pdb(ligand_file_name + '.pdb', open_def_file, custom_def_file)
        open_def_file.close()
        custom_def_file.close()
        generate_massfile(ligand_file_name + '.pdb', os.path.join(temp_path, 'NRGTEN', f'{lig}.masses'))
        with open(custom_def_path, "r") as t2:
            texto = t2.readlines()
            definition = texto[-1][:-2] + '\n'
        with open(os.path.join(temp_path,'NRGTEN', f'{lig}.atomtypes'), "w") as t3:
            t3.write(definition)
    return list(het_dic)


def generate_massfile(pdb_filename, mass_filename):
    atoms = []
    with open(pdb_filename) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('ATOM') or line.startswith("HETATM"):
            atoms.append(Atom(line))
    xyz_data = np.zeros((len(atoms), 3))
    for i, atom in enumerate(atoms):
        xyz_data[i] = atom.xyz
    centroid = np.array([np.mean(xyz_data[:, 0]), np.mean(xyz_data[:, 1]), np.mean(xyz_data[:, 2])])
    medoid = None
    mindist = float('Inf')
    other_atom = None
    other_flag = True
    for atom in atoms:
        dist = np.linalg.norm(atom.xyz - centroid)
        if dist < mindist:
            mindist = dist
            medoid = atom
        elif other_flag:
            other_atom = atom
            other_flag = False
    if other_flag:
        other_atom = atoms[0]
    with open(mass_filename, "w") as f:
        f.write("CONNECT: {} -> {}\n".format(medoid.name, other_atom.name))
        f.write("N_MASSES: 1\n")
        f.write("CENTER: {}\n".format(medoid.name))
        f.write("NAME: {}\n".format(medoid.name))


def write_b_factor(target, dyna_sig, temp_path, labels):
    target_file = os.path.join(temp_path,'NRGTEN', f'{target}.pdb')
    b_factor_dict = {}
    with open(os.path.join(temp_path, 'NRGTEN',target+'_dynasig.txt'),'w') as t1:
        for res in range(len(labels)):
            t1.write('{} {}\n'.format(labels[res],dyna_sig[res]))
    for i in range(len(dyna_sig)):
        key = '{}_{}_{}'.format(labels[i].split('|')[0][:3], labels[i].split('|')[2], labels[i].split('|')[1])
        b_factor_dict[key] = dyna_sig[i]
    with open(target_file, 'r') as t1:
        texto = t1.readlines()
        texto_list = []
        for line in texto:
            if 'HETATM' in line or 'ATOM' in line:
                key = '{}_{}_{}'.format(line[17:20], line[21], int(line[22:26]))
                number = b_factor_dict[key]
                lined = f"{number:.2f}"
                lined = lined.rjust(6)
                texto_list.append(line[:60] + lined + line[66:])
            else:
                texto_list.append(line)
    output_path = os.path.join(temp_path,'NRGTEN', f'{target}_dynasig.pdb')
    with open(output_path, 'w') as t2:
        for line in texto_list:
            t2.write(line)
    return b_factor_dict





def dynamical_signature(target, lig, target_2, beta, main_folder_path, temp_path):
    target_file = os.path.join(temp_path,'NRGTEN', f'{target}.pdb')
    if target_2 == 'None':
        cmd.save(target_file, target)
        b_fact_dict = run_dynamical_signature(target_file, beta, main_folder_path, temp_path)[0]
        cmd.disable(target)
        cmd.load(target_file[:-4] + '_dynasig.pdb')
        cmd.spectrum(selection=os.path.basename(target_file[:-4] + '_dynasig.pdb')[:-4], palette='blue_white_red',
                     expression='b')
        if lig != 'None':
            output_file = os.path.join(temp_path,'NRGTEN', f'no_lig_{target}.pdb')
            remove_selection_and_save(target, lig, output_file)
            dyna_ob = run_dynamical_signature(output_file, beta, main_folder_path, temp_path)
            dyna_sig_no_lig = dyna_ob[1]
            model_no_lig = dyna_ob[2]

            # Use os.path.split() and os.path.splitext()
            _, filename = os.path.split(output_file)
            key_base, _ = os.path.splitext(filename)

            for b_factor in range(len(dyna_sig_no_lig)):
                mass = model_no_lig.get_mass_labels()[b_factor]
                key = '{}_{}_{}'.format(mass.split('|')[0][:3], mass.split('|')[2], mass.split('|')[1])
                dyna_sig_no_lig[b_factor] = (b_fact_dict[key] / dyna_sig_no_lig[b_factor]) - 1

            write_b_factor(key_base, dyna_sig_no_lig, temp_path, model_no_lig.get_mass_labels())
            cmd.load(output_file[:-4] + '_dynasig.pdb')
            cmd.spectrum(selection=key_base + '_dynasig', palette='blue_white_red', expression='b', minimum=-1,
                         maximum=1)
    else:
        cmd.save(target_file, target)
        b_fact_dict = run_dynamical_signature(target_file, beta, main_folder_path, temp_path)[1]
        cmd.disable(target)
        cmd.load(target_file[:-4] + '_dynasig.pdb')
        cmd.spectrum(selection=os.path.basename(target_file[:-4] + '_dynasig.pdb')[:-4], palette='blue_white_red',
                     expression='b')
        object_list = []
        diff_list=[]
        for state in range(cmd.count_states(target_2)):
            output_file = os.path.join(temp_path,'NRGTEN' ,f'{target_2}_{state}.pdb')
            cmd.save(output_file, target_2, state=state + 1)

            diff = compare_residues(target_file, output_file)
            diff_list.append(diff)
            os.rename(output_file, os.path.join(temp_path, 'NRGTEN', f'{target_2}_{diff}.pdb'))
            output_file = os.path.join(temp_path, 'NRGTEN', f'{target_2}_{diff}.pdb')

            dyna_ob = run_dynamical_signature(output_file, beta, main_folder_path, temp_path)
            dyna_sig_no_lig = dyna_ob[1]
            model_no_lig = dyna_ob[2]
            _, filename = os.path.split(output_file)
            key_base, _ = os.path.splitext(filename)

            for b_factor in range(len(dyna_sig_no_lig)):
                dyna_sig_no_lig[b_factor] = (b_fact_dict[b_factor] / dyna_sig_no_lig[b_factor]) - 1
            dyna_sig_no_lig = standardize_to_minus1_plus1(dyna_sig_no_lig)
            plt.plot(dyna_sig_no_lig, label=diff)
            write_b_factor(key_base, dyna_sig_no_lig, temp_path, model_no_lig.get_mass_labels())
            cmd.load(os.path.join(temp_path,'NRGTEN', f'{key_base}_dynasig.pdb'), f'{target_2}_dynasigdif_{diff}')
            object_list.append(f'{target_2}_dynasigdif_{diff}')
        for state in diff_list:
            cmd.spectrum(selection=f'{target_2}_dynasigdif_{state}', palette='blue_white_red', expression='b',
                         minimum=-1, maximum=1)
        create_group(f'{target_2}_dynasigdif', object_list)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.savefig(os.path.join(target_file[:-4]+'.png'), dpi=300, bbox_inches='tight')
        plt.show()


def create_group(group_name, object_list):
    members = ', '.join(object_list)
    cmd.group(group_name, members)
    return 0


def run_dynamical_signature(target_file, beta, main_folder_path, temp_path):
    _, filename = os.path.split(target_file)
    target, _ = os.path.splitext(filename)
    model = encom_model(target_file, main_folder_path, temp_path)
    dyna_sig = model.compute_bfactors_boltzmann(beta=float(beta))
    b_fact_dict = write_b_factor(target, dyna_sig, temp_path, model.get_mass_labels())
    return [b_fact_dict, dyna_sig, model]


def conformational_ensemble(target, modes_list, step, max_conf, max_disp, opt, main_folder_path, temp_path,form):
    modes_list = list(map(int, modes_list.split(',')))
    target_file = os.path.join(temp_path,'NRGTEN', f'{target}.pdb')
    cmd.save(target_file, target)
    model = encom_model(target_file, main_folder_path, temp_path)
    ensemble_path = os.path.join(temp_path,'NRGTEN', f"{target}_conf_ensemble.pdb")
    model.build_conf_ensemble(modes_list, ensemble_path, step=float(step), max_displacement=float(max_disp),
                              max_conformations=int(max_conf))
    if opt:
        model_states(ensemble_path, target, temp_path, main_folder_path,form)

    else:
        cmd.load(ensemble_path)
        cmd.show('cartoon', f'{target}_conf_ensemble')



