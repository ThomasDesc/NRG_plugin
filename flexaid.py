import os
from pymol import cmd
import shutil
import subprocess


def count_flex(ligand_inp_file_path):
    with open(ligand_inp_file_path, 'r') as t:
        texto = t.readlines()
        count = 0
        for line in texto:
            if line[0] == 'F':
                count += 1
    return count


def write_config(target_inp_path, cleft, ligand_inp_path, max_results, flexaid_output_path):
    with open("template_CONFIG.inp", "r") as t1:
        lines = t1.readlines()
    config_file_output_path = os.path.join(flexaid_output_path, 'config.inp')
    with open(config_file_output_path, "w") as t2:
        for line_counter, line in enumerate(lines):
            if line_counter == 0:
                t2.write(f'PDBNAM {target_inp_path}\n')
            if line_counter == 1:
                t2.write(f'INPLIG {ligand_inp_path}\n')
            elif line_counter == 3:
                t2.write(f"RNGOPT LOCCLF {cleft}\n")
            elif line_counter == 6:
                for flex in range(count_flex(ligand_inp_path)):
                    t2.write(f'OPTIMZ 9999 - {str(flex + 1)}\n')
                t2.write(line)
            elif line_counter == 7:
                t2.write(f'STATEP "{flexaid_output_path}"\n')
            elif line_counter == 18:
                t2.write(f'MAXRES {max_results}\n')
            else:
                t2.write(line)


def process_ligand(process_ligand_path, input_path, process_ligand_output, istarget=False):
    # file_name = os.path.splitext(os.path.basename(input_path))[0]
    if istarget:
        process_ligand_command = f'"{process_ligand_path}" -f "{input_path}" -target'
    else:
        process_ligand_command = f'"{process_ligand_path}" -f "{input_path}" --atom_index 90000 -ref'
    print(process_ligand_command)
    subprocess.run(process_ligand_command, shell=True)
    # if istarget:
    #     print('removing: ', os.path.splitext(input_path)[0] + '.mol2.tmp')
    #     os.remove(os.path.splitext(input_path)[0] + '.mol2.tmp')
    #     print('\nto: ', os.path.join(process_ligand_output, 'flexaid_target.inp.pdb'))
    #     shutil.move(os.path.splitext(input_path)[0] + '.inp.pdb', os.path.join(process_ligand_output, f"{file_name}.inp.pdb"))


def run_flexaid(flexaid_output_path, form, cleft_save_path, process_ligand_path, flexaid_path):
    process_ligand_output = os.path.join(flexaid_output_path, 'process_ligand_output')
    max_results = 10
    if os.path.isdir(process_ligand_output):
        shutil.rmtree(process_ligand_output)
    binding_site_path = os.path.join(flexaid_output_path, 'binding_site_sph_.pdb')
    os.mkdir(process_ligand_output)
    target_name = form.flexaid_select_target.currentText()
    target_save_path = os.path.join(flexaid_output_path, 'flexaid_target.pdb')
    cmd.save(target_save_path, target_name)
    ligand_name = form.flexaid_select_ligand.currentText()
    ligand_save_path = os.path.join(flexaid_output_path, 'flexaid_ligand.pdb')
    cmd.save(ligand_save_path, ligand_name)
    binding_site_name = form.flexaid_select_binding_site.currentText()
    shutil.copy(os.path.join(cleft_save_path, binding_site_name + '.pdb'), binding_site_path)
    process_ligand(process_ligand_path, target_save_path, process_ligand_output, istarget=True)
    process_ligand(process_ligand_path, ligand_save_path, process_ligand_output)
    target_inp_path = os.path.splitext(target_save_path)[0] + '.inp.pdb'
    ligand_inp_path = os.path.splitext(target_save_path)[0] + '.inp'
    write_config(target_inp_path, binding_site_path, ligand_inp_path, max_results, flexaid_output_path)
    print('Done')
