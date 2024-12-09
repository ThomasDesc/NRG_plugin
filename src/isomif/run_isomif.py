from pymol import cmd
import os
import sys
from general_functions import output_message
import subprocess
import time

def get_residue_string(selection_name):
    residue_info = {'resn': '', 'resi': '', 'chain': ''}
    cmd.iterate(selection_name, 'residue_info.update({"resn": resn, "resi": resi, "chain": chain})', space=locals())
    residue_string = f"{residue_info['resn']}{residue_info['resi']}{residue_info['chain']}-"

    return residue_string

def run_cleft_lig(target,target_file,lig, get_cleft, ISOMIF_res):
    lig_string=get_residue_string(lig)
    command=f'{get_cleft} -p {target_file} -o {os.path.join(ISOMIF_res,target)} -s -a {lig_string}'
    print(command)
    os.system(command)
    for file in os.listdir(ISOMIF_res):
        if target in  file and lig_string in file and '_sph' in file:
            return os.path.join(ISOMIF_res,file)


def run_mif_pre(isomif_temp_result_path, target, ligand_name, cleft_name, output_box, temp_path, mif_binary_path,
                mifview_binary_path, python_version, isomif_deps_path):
    target_file = os.path.join(isomif_temp_result_path, f'{target}.pdb')
    cmd.save(target_file, target)
    cleft_save_path = os.path.join(isomif_temp_result_path, f'{cleft_name}.pdb')
    cmd.save(cleft_save_path, cleft_name)
    ligand_residue_name = 'None'
    if ligand_name != "None":
        ligand_residue_name = get_residue_string(ligand_name)
    run_mif(target, output_box, temp_path, cleft_save_path, mif_binary_path, mifview_binary_path, isomif_temp_result_path,
            python_version, ligand_residue_name, isomif_deps_path)


def mif_plot(form, binary_folder_path, binary_suffix, install_dir):
    output_box = form.output_box
    python_version = f'{sys.version.split(".")[0]}.{sys.version_info[1]}'

    mif_binary_path = os.path.join(binary_folder_path, f'mif{binary_suffix}')
    isomif_binary_path = os.path.join(binary_folder_path, f'isomif{binary_suffix}')
    mifView_binary_path = os.path.join(install_dir,'src','isomif', f'mifView.py')
    isoMifView_binary_path = os.path.join(install_dir,'src','isomif', f'isoMifView.py')
    isomif_deps_path = os.path.join(install_dir,'deps','isomif')
    temp_path = form.temp_line_edit.text()
    isomif_temp_result_path=os.path.join(temp_path,'ISOMIF')

    target_1 = form.ISOMIF_select_target_object_1.currentText()
    cleft_name_1 = form.ISOMIF_select_cleft_object_1.currentText()
    ligand_name_1 = form.ISOMIF_select_ligand_object_1.currentText()

    target_2 = form.ISOMIF_select_target_object_2.currentText()
    cleft_name_2 = form.ISOMIF_select_cleft_object_2.currentText()
    ligand_name_2 = form.ISOMIF_select_ligand_object_2.currentText()

    if target_1:
        if cleft_name_1:
            run_mif_pre(isomif_temp_result_path, target_1, ligand_name_1, cleft_name_1, output_box, temp_path,
                        mif_binary_path, mifView_binary_path, python_version, isomif_deps_path)
    # if target_2 != "None":
    #     if cleft_name_2 != "None":
    #         run_mif_pre(isomif_temp_result_path, target_2, ligand_name_2, cleft_name_2, output_box, temp_path,
    #                     mif_binary_path, mifView_binary_path, python_version, isomif_deps_path)
    #
    #         run_isomif(target_1, target_2, isomif_binary_path, isoMifView_binary_path, isomif_temp_result_path, python_version)
    #         with open(os.path.join(isomif_temp_result_path,f'iso_{target_1}_h_match_{target_2}_h.isomif'),'r') as t1:
    #             texto=t1.readlines()
    #             for line in texto:
    #                 if 'REMARK CLIQUE' in line:
    #                     output_message(form.output_box,f'ISOMIF results: {line[15:]}','valid')


def run_isomif(target,target_2, isomif_binary_path, isoMifView_binary_path, ISOMIF_res, python_version):

    command_isomif = f'{isomif_binary_path} -p1 {os.path.join(ISOMIF_res,target+"_h.mif")} -p2 {os.path.join(ISOMIF_res,target_2+"_h.mif")} -o {os.path.join(ISOMIF_res,"iso_")} -c 2'
    print(command_isomif)
    os.system(command_isomif)

    isomif_file = os.path.join(ISOMIF_res,f'iso_{target}_h_match_{target_2}_h.isomif')
    command_isomifView = [f'python{python_version}',  isoMifView_binary_path, '-m', isomif_file, '-o', os.path.join(ISOMIF_res,"view_"), '-g', '2']
    print(command_isomifView)
    _process = subprocess.Popen(command_isomifView)
    while _process.poll() is None:
        time.sleep(0.1)


    # os.system(command_isomifView)

    initial_objects = set(cmd.get_object_list())
    cmd.load(os.path.join(ISOMIF_res,f'view_{target}_h_{target_2}_h.pml'))
    new_objects = set(cmd.get_object_list()) - initial_objects
    cmd.group(f'isomif_{target}_{target_2}'," ".join(new_objects))
    cmd.group('isomif_results',f'isomif_{target}_{target_2}')


def run_mif(target, output_box, temp_path, cleft_file_path, mif_binary_path, mifView_binary_path, ISOMIF_res, python_version, lig_str, deps):


    target_file=os.path.join(temp_path,'ISOMIF',f'{target}.pdb')

    output_message(output_box, 'Running ISOMIF...', 'valid')

    cmd.create(f'{target}_h', target)
    cmd.h_add(f'{target}_h')
    cmd.save(target_file[:-4]+'_h.pdb',f'{target}_h')
    cmd.delete(f'{target}_h')

    command_mif = f'{mif_binary_path} -p {target_file[:-4] + "_h.pdb"} -g {cleft_file_path} -o {ISOMIF_res} -s 1 -dp {deps}'
    if lig_str != 'None':
        command_mif += f' -l {lig_str}'
        if lig_str[-2:] == '9-':
             command_mif += '-'
    print(command_mif)

    os.system(command_mif)

    command_view=f'python{python_version} {mifView_binary_path} -m {target_file[:-4]}_h.mif -o {ISOMIF_res}'
    os.system(command_view)
    cmd.load(target_file[:-4]+'_h.pml')
    cmd.delete(target+'_h')

    try:
        cmd.set_name('neg_100', f'{target}_neg_100')
        cmd.group(f'{target}_isomif', f'{target}_neg_100')
    except:
        print('no neg_100')
    try:
        cmd.set_name('don_100', f'{target}_don_100')
        cmd.group(f'{target}_isomif', f'{target}_don_100')
    except:
        print('no don_100')
    try:
        cmd.set_name('acc_100', f'{target}_acc_100')
        cmd.group(f'{target}_isomif', f'{target}_acc_100')
    except:
        print('no acc_100')
    try:
        cmd.set_name('pos_100', f'{target}_pos_100')
        cmd.group(f'{target}_isomif', f'{target}_pos_100')
    except:
        print('no pos_100')
    try:
        cmd.set_name('arm_100', f'{target}_arm_100')
        cmd.group(f'{target}_isomif', f'{target}_arm_100')
    except:
        print('no arm_100')
    try:
        cmd.set_name('hyd_100', f'{target}_hyd_100')
        cmd.group(f'{target}_isomif', f'{target}_hyd_100')
    except:
        print('no hyd_100')
    try:
        cmd.set_name('100', f'{target}_100')
        cmd.group(f'{target}_isomif', f'{target}_100')
        cmd.disable(f'{target}_100')
    except:
        print('no 100')

    cmd.group('isomif_results', f'{target}_isomif')

    #mif_to_pml(target_file[:-4]+'_h.mif',ISOMIF_res,"")

