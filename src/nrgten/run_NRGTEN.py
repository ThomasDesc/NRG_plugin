import shutil
from modeller import Environ, Model, Alignment
from modeller.automodel import *
from modeller.selection import selection
from nrgten.encom import ENCoM
from nrgten.atom import Atom
from pymol import cmd
from ligand_atomtypes import add_pdb
from clean_structure import main as clean_structure
import os
import numpy as np
from run_Surfaces import create_ligand_file
import matplotlib.pyplot as plt


def model_states(input_file,target,temp_path,main_folder_path):
    model_data = []
    current_model = []
    model_dir=temp_path+'/'+target+'_models'
    os.mkdir(model_dir)
    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('MODEL'):
                if current_model:
                    model_data.append(current_model)
                current_model = [line]
            elif line.startswith('ENDMDL'):
                current_model.append(line)
                model_data.append(current_model)
                current_model = []
            else:
                if current_model:
                    current_model.append(line)
    if current_model:
        model_data.append(current_model)

    for i, model in enumerate(model_data):
        with open(f'{model_dir}/model_{i + 1}.pdb', 'w') as output_file:
            output_file.writelines(model)
        code = f'{model_dir}/model_{i+1}'
        cleaned_file_path=code+'_clean.pdb'
        clean_pdb_file = open(cleaned_file_path, "w")
        clean_structure(code+'.pdb',main_folder_path+'/surfaces_defs/AMINO_FlexAID.def',clean_pdb_file)
        clean_pdb_file.close()
        model_state(cleaned_file_path, model_dir)
        cmd.load(cleaned_file_path,f'{target}_ensemble',state=i+1)
        cmd.show('cartoon',f'{target}_ensemble')
        if i == 10:
            break

def model_state(cleaned_file_path,model_dir):
    pdb_path=cleaned_file_path

    # Step 1: Load the environment
    env = Environ()

    # Step 2: Read in the sequence of the target structure (pdb.pdb) without specifying chain
    mdl = Model(env, file=pdb_path)  # Use all chains
    aln = Alignment(env)
    aln.append_model(mdl, align_codes='pdb_target', atom_files=pdb_path)

    # Step 3: Append the template structure (pdb.pdb) without specifying chain
    mdl = Model(env, file=pdb_path)  # Use all chains
    aln.append_model(mdl, align_codes='pdb_template', atom_files=pdb_path)

    # Step 4: Align the target and template sequences
    aln.align2d()

    # Step 5: Write the alignment file for Modeller
    aln.write(file='alignment.ali', alignment_format='PIR')

    # Step 6: Define a custom automodel class
    class MyModel(automodel):
        def select_atoms(self):
            return selection(self)

    # Specify the directory to save the models
    output_dir = 'output_models'  # Replace with your desired output directory

    # Create the directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Step 7: Initialize the model building
    a = MyModel(env,
                alnfile='alignment.ali',  # Alignment filename
                knowns='pdb_template',  # Code of the template
                sequence='pdb_target')  # Code of the target sequence

    a.starting_model = 1
    a.ending_model = 1
    a.very_fast()
    a.make()
    os.rename('pdb_target.B99990001.pdb',cleaned_file_path)



def remove_selection_and_save(object_name, selection, output_file):
    cmd.create("object_without_selection", f"{object_name} and not ({selection})")
    cmd.save(output_file, "object_without_selection")
    cmd.delete("object_without_selection")


def flex_aid_matrix(main_folder_path):
    matrix=np.array([np.zeros(41)]*41)
    with open(main_folder_path+"/surfaces_defs/FlexAID.dat","r") as t1:
        texto=t1.readlines()
        for line in texto:
            matrix[int(line.split("=")[0].split('-')[0])][int(line.split("=")[0].split('-')[1])]=-float(line.split("=")[1])
    return(matrix)

def encom_model(target_file, main_folder_path,temp_path):
    list_het=find_het(target_file,temp_path, main_folder_path)
    atype_lists=[main_folder_path+'/nrgten_defs/amino_acids.atomtypes']
    mass_lists=[main_folder_path+'/nrgten_defs/amino_acids.masses']
    for lig in list_het:
        atype_lists.append(temp_path+'/{}.atomtypes'.format(lig))
        mass_lists.append(temp_path+'/{}.masses'.format(lig))

    matrix = flex_aid_matrix(main_folder_path)
    model = ENCoM(target_file,interact_mat=matrix,atypes_list=atype_lists,massdef_list=mass_lists)
    return model


def standardize_to_minus1_plus1(data):
    max_abs_value = max(abs(x) for x in data)

    standardized_data = [x / max_abs_value for x in data]
    return standardized_data

def find_het(target_file,temp_path,main_folder_path):
    het_dic={}
    with open(target_file, "r") as t1:
        texto=t1.readlines()
        for line in texto:
            if 'HETATM' in line:
                het_dic[line[17:20]]='1'
    for lig in list(het_dic):
        def_file = os.path.join(main_folder_path, "surfaces_defs", 'AMINO_FlexAID.def')
        flexaid_dat_path = os.path.join(main_folder_path, "surfaces_defs", 'FlexAID.dat')
        open_def_file = open(def_file, "r")
        ligand_file_name = os.path.join(os.path.dirname(target_file), lig)
        create_ligand_file(target_file, ligand_file_name)
        custom_def_path = os.path.join(temp_path, f'custom_{os.path.basename(def_file)}')
        custom_def_file = open(custom_def_path, 'w')
        add_pdb(ligand_file_name + '.pdb', open_def_file, custom_def_file)
        open_def_file.close()
        custom_def_file.close()
        generate_massfile(ligand_file_name+'.pdb', temp_path+'/{}.masses'.format(lig))
        with open(custom_def_path, "r") as t2:
            texto=t2.readlines()
            definition=texto[-1][:-2]+'\n'
        with open(temp_path+"/{}.atomtypes".format(lig),"w") as t3:
            t3.write(definition)
    return(list(het_dic))

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

def write_b_factor(target,dyna_sig,temp_path,labels):
    target_file = temp_path + '/{}.pdb'.format(target)
    b_factor_dict={}
    for i in range(len(dyna_sig)):
        key='{}_{}_{}'.format(labels[i].split('|')[0][:3],labels[i].split('|')[2],labels[i].split('|')[1])
        b_factor_dict[key]=dyna_sig[i]
    with open(target_file,'r') as t1:
        texto=t1.readlines()
        texto_list=[]
        for line in texto:
            if 'HETATM' in line or 'ATOM' in line:
                key='{}_{}_{}'.format(line[17:20], line[21], int(line[23:26]))
                number=b_factor_dict[key]
                lined=f"{number:.2f}"
                lined=lined.rjust(6)
                texto_list.append(line[:60]+lined+line[66:])
            else:
                texto_list.append(line)
    with open(target_file[:-4]+'_dynasig.pdb','w') as t2:
        for line in texto_list:
            t2.write(line)
    return(b_factor_dict)

def refresh_dropdown_NRG(dropdown_to_refresh, output_box, filter_for='', no_warning=False):
    list_pymol_objects = cmd.get_names('all')
    list_pymol_objects.append('None')
    if filter_for != '':
        list_pymol_objects = [x for x in list_pymol_objects if filter_for in x]
    if len(list_pymol_objects) == 0 and no_warning is False:
        output_message(output_box, 'No objects found', 'warning')
    dropdown_to_refresh.clear()
    dropdown_to_refresh.addItems(list_pymol_objects)
    if len(list_pymol_objects) > 0:
        dropdown_to_refresh.setCurrentText(list_pymol_objects[-1])

def dynamical_signature(target,lig,target_2,beta,main_folder_path,temp_path):
    if target_2 == 'None':
        target_file = temp_path + '/{}.pdb'.format(target)
        cmd.save(target_file, target)
        b_fact_dict=run_dynamical_signature(target_file,beta,main_folder_path,temp_path)[0]
        cmd.disable(target)
        cmd.load(target_file[:-4] + '_dynasig.pdb')
        cmd.spectrum(selection=os.path.basename(target_file[:-4] + '_dynasig.pdb')[:-4],
                            palette='blue_white_red', expression='b')
        if lig!='None':
            output_file = temp_path + '/no_lig_' + target + '.pdb'
            remove_selection_and_save(target, lig, output_file)
            dyna_ob=run_dynamical_signature(output_file,beta,main_folder_path,temp_path)
            dyna_sig_no_lig=dyna_ob[1]
            model_no_lig=dyna_ob[2]
            for b_factor in range(len(dyna_sig_no_lig)):
                mass = model_no_lig.get_mass_labels()[b_factor]
                key = '{}_{}_{}'.format(mass.split('|')[0][:3], mass.split('|')[2], mass.split('|')[1])
                dyna_sig_no_lig[b_factor] = (b_fact_dict[key] / dyna_sig_no_lig[b_factor]) - 1
            write_b_factor(output_file[:-4].split('/')[-1], dyna_sig_no_lig, temp_path,
                               model_no_lig.get_mass_labels())
            cmd.load(output_file[:-4] + '_dynasig.pdb')
            cmd.spectrum(selection=output_file[:-4].split('/')[-1] + '_dynasig', palette='blue_white_red',
                             expression='b', minimum=-1, maximum=1)
    else:
        target_file = temp_path + '/{}.pdb'.format(target)
        cmd.save(target_file, target)
        b_fact_dict = run_dynamical_signature(target_file, beta, main_folder_path, temp_path)[1]
        cmd.disable(target)
        cmd.load(target_file[:-4] + '_dynasig.pdb')
        cmd.spectrum(selection=os.path.basename(target_file[:-4] + '_dynasig.pdb')[:-4],
                     palette='blue_white_red', expression='b')
        object_list=[]
        for state in range(cmd.count_states(target_2)):
            output_file=temp_path + '/{}_{}.pdb'.format(target_2,state)
            cmd.save(output_file, target_2, state=state+1)
            dyna_ob = run_dynamical_signature(output_file, beta, main_folder_path, temp_path)
            dyna_sig_no_lig = dyna_ob[1]
            model_no_lig = dyna_ob[2]
            for b_factor in range(len(dyna_sig_no_lig)):
                dyna_sig_no_lig[b_factor] = (b_fact_dict[b_factor] / dyna_sig_no_lig[b_factor]) - 1
            dyna_sig_no_lig=standardize_to_minus1_plus1(dyna_sig_no_lig)
            plt.plot(dyna_sig_no_lig, label=str(state+1))
            write_b_factor(output_file[:-4].split('/')[-1], dyna_sig_no_lig, temp_path,
                           model_no_lig.get_mass_labels())
            cmd.load(output_file[:-4] + '_dynasig.pdb',target_2+'_dynasigdif_'+str(state+1))
            object_list.append(target_2+'_dynasigdif_'+str(state+1))
        for state in range(cmd.count_states(target_2)):
            cmd.spectrum(selection=target_2+'_dynasigdif_'+str(state+1), palette='blue_white_red',
                         expression='b', minimum=-1, maximum=1)
        create_group(target_2+'_dynasigdif', object_list)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.show()

def create_group(group_name, object_list):
    members = ', '.join(object_list)
    cmd.group(group_name, members)
    return 0



def run_dynamical_signature(target_file,beta,main_folder_path,temp_path):
    target=target_file.split('/')[-1][:-4]
    model=encom_model(target_file, main_folder_path,temp_path)
    dyna_sig=model.compute_bfactors_boltzmann(beta=float(beta))
    b_fact_dict=write_b_factor(target, dyna_sig,temp_path, model.get_mass_labels())
    return [b_fact_dict,dyna_sig,model]


def conformational_ensemble(target,modes_list,step,max_conf,max_disp,opt,main_folder_path,temp_path):
    modes_list=list(map(int,modes_list.split(',') ))
    target_file = temp_path + '/{}.pdb'.format(target)
    cmd.save(target_file, target)
    model = encom_model(target_file, main_folder_path, temp_path)
    model.build_conf_ensemble(modes_list, temp_path+"/{}_conf_ensemble.pdb".format(target),step=float(step), max_displacement=float(max_disp),max_conformations=int(max_conf))
    if opt:
        model_states(temp_path+"/{}_conf_ensemble.pdb".format(target),target,temp_path,main_folder_path)
    else:
        cmd.load(temp_path+"/{}_conf_ensemble.pdb".format(target))
        cmd.show('cartoon','{}_conf_ensemble'.format(target))

