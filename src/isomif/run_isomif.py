from src.surfaces.run_Surfaces import flex_res, process_result_flexaid
from pymol import cmd
import os
from src.isomif.MifView import mif_to_pml
import shutil


def get_residue_string(selection_name):
    # Dictionary to store residue information
    residue_info = {'resn': '', 'resi': '', 'chain': ''}

    # Use cmd.iterate to get the residue name, number, and chain for the single residue in the selection
    cmd.iterate(selection_name, 'residue_info.update({"resn": resn, "resi": resi, "chain": chain})', space=locals())

    # Format the result as 'resn + resi + chain', e.g., 'ARG123A'
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



def mif_plot(form, outputbox, binary_folder_path, binary_suffix, operating_system):
    processlig_binary_path = mif_binary_path = os.path.join(binary_folder_path,f'Process_Ligand{binary_suffix}')
    mif_binary_path = os.path.join(binary_folder_path, f'mif{binary_suffix}')
    isomif_binary_path = os.path.join(binary_folder_path, f'isomif{binary_suffix}')
    mifView_binary_path = os.path.join(binary_folder_path, f'mifView{binary_suffix}')
    isoMifView_binary_path = os.path.join(binary_folder_path, f'isoMifView{binary_suffix}')
    get_cleft_bineary_path=os.path.join(binary_folder_path, f'GetCleft{binary_suffix}')
    temp_path = form.temp_line_edit.text()
    ISOMIF_res=os.path.join(temp_path,'ISOMIF')
    lig_name = form.ISOMIF_select_cleft.currentText()
    target = form.ISOMIF_select_target.currentText()
    lig_name_2 = form.ISOMIF_select_cleft_1.currentText()
    target_2 = form.ISOMIF_select_target_1.currentText()
    if form.ISOMIF_select_target:
        if form.ISOMIF_select_cleft:
            target_file = os.path.join(temp_path, 'ISOMIF', f'{target}.pdb')
            cmd.save(target_file, target)

            cleft_name=os.path.join(ISOMIF_res,lig_name+'.pdb')
            cmd.save(os.path.join(ISOMIF_res,lig_name+'.pdb'),lig_name)
            #shutil.copyfile(os.path.join(temp_path,'GetCleft','Clefts',lig_name+'.pdb'),cleft_name)

            #cleft_name=run_cleft_lig(target,target_file,lig_name,get_cleft_bineary_path,ISOMIF_res)
            run_mif(target, form, temp_path, cleft_name, mif_binary_path, mifView_binary_path, ISOMIF_res)
    if target_2!="None":
        if lig_name_2!="None":
            target_file_2 = os.path.join(temp_path, 'ISOMIF', f'{target}.pdb')
            cmd.save(target_file_2, target_2)

            cleft_name_2 = os.path.join(ISOMIF_res, lig_name_2+'.pdb')
            cmd.save(os.path.join(ISOMIF_res,lig_name_2+'.pdb'),lig_name_2)

            #shutil.copyfile(os.path.join(temp_path, 'GetCleft', 'Clefts', lig_name_2+'.pdb'), cleft_name_2)

            #cleft_name_2 = run_cleft_lig(target_2, target_file_2, lig_name_2, get_cleft_bineary_path, ISOMIF_res)
            run_mif(target_2, form, temp_path, cleft_name_2, mif_binary_path, mifView_binary_path, ISOMIF_res)
            run_isomif(target,target_2,cleft_name, cleft_name_2,form,temp_path, isomif_binary_path, isoMifView_binary_path, ISOMIF_res)

def run_isomif(target,target_2,cleft_name, cleft_name_2,form,temp_path, isomif_binary_path, isoMifView_binary_path, ISOMIF_res):

    command_isomif = f'{isomif_binary_path} -p1 {os.path.join(ISOMIF_res,target+"_h.mif")} -p2 {os.path.join(ISOMIF_res,target_2+"_h.mif")} -o {os.path.join(ISOMIF_res,"iso_")} -c 2'
    os.system(command_isomif)

    print(command_isomif)

    isomif_file=os.path.join(ISOMIF_res,f'iso_{target}_h_match_{target_2}_h.isomif')
    command_isomifView=f'{isoMifView_binary_path} -m {isomif_file} -o {os.path.join(ISOMIF_res,'view_')} -g 2'

    print(command_isomifView)
    os.system(command_isomifView)


def run_mif(target,form,temp_path,cleft_file,mif_binary_path,mifView_binary_path, ISOMIF_res):
            target_file=os.path.join(temp_path,'ISOMIF',f'{target}.pdb')

            form.output_box.append(f'Running ISOMIF...')

            cmd.create(f'{target}_h', target)
            cmd.h_add(f'{target}_h')
            cmd.save(target_file[:-4]+'_h.pdb',f'{target}_h')
            cmd.delete(f'{target}_h')

            command_mif=f'{mif_binary_path} -p {target_file[:-4]+'_h.pdb'} -g {cleft_file} -o {ISOMIF_res} -s 1'

            print(command_mif)

            os.system(command_mif)

            command_view=f'{mifView_binary_path} -m {target_file[:-4]}_h.mif -o {ISOMIF_res}'
            os.system(command_view)
            cmd.load(target_file[:-4]+'_h.pml')
            cmd.delete(target+'_h')

            cmd.set_name('neg_100',f'{target}_neg_100')
            cmd.group(f'{target}_isomif', f'{target}_neg_100')

            cmd.set_name('don_100',f'{target}_don_100')
            cmd.group(f'{target}_isomif', f'{target}_don_100')

            cmd.set_name('acc_100',f'{target}_acc_100')
            cmd.group(f'{target}_isomif', f'{target}_acc_100')

            cmd.set_name('pos_100',f'{target}_pos_100')
            cmd.group(f'{target}_isomif', f'{target}_pos_100')

            cmd.set_name('arm_100',f'{target}_arm_100')
            cmd.group(f'{target}_isomif', f'{target}_arm_100')

            cmd.set_name('hyd_100',f'{target}_hyd_100')
            cmd.group(f'{target}_isomif', f'{target}_hyd_100')

            cmd.set_name('100', f'{target}_100')
            cmd.group(f'{target}_isomif', f'{target}_100')

            cmd.group('isomif_results', f'{target}_isomif')

            #mif_to_pml(target_file[:-4]+'_h.mif',ISOMIF_res,"")

