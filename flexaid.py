import os
from pymol import cmd
import shutil
import subprocess


def count_flex(ligand_inp_file_path):
    with open('ligands/{}.inp'.format(ligand_inp_file_path), 'r') as t:
        texto = t.readlines()
        count = 0
        for line in texto:
            if line[0] == 'F':
                count += 1
    return count


def write_config(target,cleft,mol_id, save_res, Path_to_workdir):
    with open("template_CONFIG.inp", "r") as t1:
        texto = t1.readlines()
        count = 0
    with open("config/{}_CONFIG.inp".format(mol_id), "w") as t2:
        for line in texto:
            if count == 1:
                t2.write('INPLIG {}/ligands/{}.inp\n'.format(Path_to_workdir, mol_id))
            elif count== 0:
                t2.write('PDBNAM {}/target/{}.inp.pdb\n'.format(Path_to_workdir,target.split('/')[-1][:-4]))
            elif count==3:
                t2.write("RNGOPT LOCCLF {}\n".format(cleft))
            elif count==6:
                for flex in range(count_flex(mol_id)):
                    t2.write('OPTIMZ 9999 - {}\n'.format(str(flex + 1)))
                t2.write(line)
            elif count== 7:
                t2.write('STATEP "{}"\n'.format(Path_to_workdir))
            elif count == 18:
                t2.write('MAXRES {}\n'.format(save_res))
            else:
                t2.write(line)
            count += 1


def process_ligand(process_ligand_path, input_path, process_ligand_output, istarget=False):
    if istarget:
        process_ligand_command = f'"{process_ligand_path}" -f "{input_path}" -target'
    else:
        process_ligand_command = f'"{process_ligand_path}" -f "{input_path}" -o "{os.path.join(process_ligand_output, "ligand")}"'
    print(process_ligand_command)
    subprocess.run(process_ligand_command, shell=True)
    if istarget:
        print('removing: ', os.path.splitext(input_path)[0] + '.mol2.tmp')
        os.remove(os.path.splitext(input_path)[0] + '.mol2.tmp')
        print('\nto: ', os.path.join(process_ligand_output, 'flexaid_target.inp.pdb'))
        shutil.move(os.path.splitext(input_path)[0] + '.inp.pdb', os.path.join(process_ligand_output, 'flexaid_target.inp.pdb'))


def run_flexaid(flexaid_output_path, form, cleft_save_path, process_ligand_path, flexaid_path):
    process_ligand_output = os.path.join(flexaid_output_path, 'process_ligand_output')
    if os.path.isdir(process_ligand_output):
        shutil.rmtree(process_ligand_output)
    os.mkdir(process_ligand_output)
    target_name = form.flexaid_select_target.currentText()
    target_save_path = os.path.join(flexaid_output_path, 'flexaid_target.pdb')
    cmd.save(target_save_path, target_name)
    ligand_name = form.flexaid_select_ligand.currentText()
    ligand_save_path = os.path.join(flexaid_output_path, 'flexaid_ligand.pdb')
    cmd.save(ligand_save_path, ligand_name)
    binding_site_name = form.flexaid_select_binding_site.currentText()
    shutil.copy(os.path.join(cleft_save_path, binding_site_name + '.pdb'), os.path.join(flexaid_output_path, 'binding_site_sph_.pdb'))
    process_ligand(process_ligand_path, target_save_path, process_ligand_output, istarget=True)
    process_ligand(process_ligand_path, ligand_save_path, process_ligand_output)
    print('Done')
