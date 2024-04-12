import os
from pymol import cmd
import shutil
import subprocess
import datetime
import thread_test


def submit_command(getcleft_command):
    print('submitting command')
    subprocess.run(getcleft_command, shell=True)


def count_flex(ligand_inp_file_path):
    with open(ligand_inp_file_path, 'r') as t:
        texto = t.readlines()
        count = 0
        for line in texto:
            if line[0] == 'F':
                count += 1
    return count


def write_config(target_inp_path, cleft, ligand_inp_path, max_results, flexaid_output_path, flexaid_result_path):
    with open(os.path.join(os.path.dirname(__file__), 'template_config.inp'), "r") as t1:
        lines = t1.readlines()
    config_file_output_path = os.path.join(flexaid_output_path, 'config.inp')
    flexaid_deps_path = os.path.join(os.path.dirname(__file__), 'flexaid_deps')
    matrix_path = os.path.join(flexaid_deps_path, 'MC_5p_norm_P10_M2_2.dat')
    with open(config_file_output_path, "w") as t2:
        for line in lines:
            if line.startswith('PDBNAM'):
                t2.write(f'PDBNAM {target_inp_path}\n')
            elif line.startswith('INPLIG'):
                t2.write(f'INPLIG {ligand_inp_path}\n')
            elif line.startswith('RNGOPT LOCCLF'):
                t2.write(f"RNGOPT LOCCLF {cleft}\n")
            elif line.startswith('OPTIMZ 9999 - 0'):
                t2.write(line)
                for flex in range(count_flex(ligand_inp_path)):
                    t2.write(f'OPTIMZ 9999 - {str(flex + 1)}\n')
            elif line.startswith('STATEP'):
                t2.write(f'STATEP {flexaid_result_path}\n')
            elif line.startswith('IMATRX'):
                t2.write(f'IMATRX {matrix_path}\n')
            elif line.startswith('DEPSPA'):
                t2.write(f'DEPSPA {flexaid_deps_path}\n')
            elif line.startswith('MAXRES'):
                t2.write(f'MAXRES {max_results}\n')
            else:
                t2.write(line)
    return config_file_output_path


def process_ligand(process_ligand_path, input_path, istarget=False):
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


def get_simulation_settings(form):
    setting_dictionary = {'number_chromosomes': form.input_num_chromosomes.text(),
                          'number_generations': form.input_num_generations.text()}
    return setting_dictionary


def edit_ga(ga_template_path, ga_write_path, setting_dictionary):
    with open(ga_template_path, "r") as ga_template:
        lines = ga_template.readlines()
    with open(ga_write_path, "w") as ga_write:
        for line in lines:
            if line.startswith('NUMCHROM'):
                ga_write.write(f'NUMCHROM {setting_dictionary["number_chromosomes"]}\n')
            elif line.startswith('NUMGENER'):
                ga_write.write(f'NUMGENER {setting_dictionary["number_generations"]}\n')
            else:
                ga_write.write(line)


def run_flexaid(flexaid_output_path, form, cleft_save_path, process_ligand_path, flexaid_path, simulation_folder_path):
    max_results = 10
    setting_dictionary = get_simulation_settings(form)
    date = datetime.datetime.now()
    date_time_str = date.strftime("%d-%m-%y-%I-%M-%S")
    flexaid_result_path = os.path.join(simulation_folder_path, date_time_str)
    os.mkdir(flexaid_result_path)
    flexaid_result_name_path = os.path.join(flexaid_result_path, "RESULT")
    binding_site_path = os.path.join(flexaid_output_path, 'binding_site_sph_.pdb')
    target_name = form.flexaid_select_target.currentText()
    target_save_path = os.path.join(flexaid_output_path, 'flexaid_target.pdb')
    cmd.save(target_save_path, target_name)
    ligand_name = form.flexaid_select_ligand.currentText()
    ligand_save_path = os.path.join(flexaid_output_path, 'flexaid_ligand.pdb')
    cmd.save(ligand_save_path, ligand_name)
    binding_site_name = form.flexaid_select_binding_site.currentText()
    shutil.copy(os.path.join(cleft_save_path, binding_site_name + '.pdb'), binding_site_path)
    process_ligand(process_ligand_path, target_save_path, istarget=True)
    process_ligand(process_ligand_path, ligand_save_path)
    target_inp_path = os.path.splitext(target_save_path)[0] + '.inp.pdb'
    ligand_inp_path = os.path.splitext(ligand_save_path)[0] + '.inp'
    config_file_path = write_config(target_inp_path, binding_site_path, ligand_inp_path, max_results, flexaid_output_path, flexaid_result_path)
    ga_path = os.path.join(flexaid_output_path, 'ga_inp.dat')
    edit_ga(os.path.join(os.path.dirname(__file__), 'ga_inp.dat'), ga_path, setting_dictionary)
    flexaid_command = f'"{flexaid_path}" "{config_file_path}" "{ga_path}" "{flexaid_result_name_path}"'
    with open('/Users/thomasdescoteaux/Documents/NRGSuite_Qt/' + 'flex_cmd.txt', 'w') as f:
        f.write(flexaid_command)
    print(flexaid_command)
    form.output_box.append(f'Please wait...Running Flexaid with command: \n{flexaid_command}')
    thread_test.yo(flexaid_command)
    # worker.finished.connect(lambda: load_show_cleft(cleft_save_path, color_list, form.output_box, pymol_object))

