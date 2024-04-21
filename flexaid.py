import os
from pymol import cmd
from pymol.Qt import QtGui
from pymol.Qt import QtWidgets
import shutil
import subprocess
import datetime
import thread_test
import time
import general_functions


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
    matrix_path = os.path.join(flexaid_deps_path, 'MC_st0r5.2_6.dat')
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
            elif line.startswith('TEMPOP'):
                t2.write(f'TEMPOP {os.path.join(flexaid_result_path, "temp")}\n')
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


def toggle_buttons(form, true_false):
    form.flexaid_button_pause.setEnabled(true_false)
    form.flexaid_button_abort.setEnabled(true_false)
    form.flexaid_button_stop.setEnabled(true_false)


def colour_specific_cell(table_widget, data):
    num_column = table_widget.columnCount()
    row = int(data[1]) - 1
    for column_counter, column in enumerate(range(num_column)):
        if column_counter == 0:
            item = QtWidgets.QTableWidgetItem()
            item.setBackground(QtGui.QColor(data[column_counter]))
        else:
            item = QtWidgets.QTableWidgetItem()
            item = QtWidgets.QTableWidgetItem(str(data[column_counter]))
        table_widget.setItem(row, column, item)


def receive_list(table_list):
    colour_specific_cell(*table_list)


def run_flexaid_worker(command, form, simulation_folder, hex_colour_list, max_generations):
    worker = thread_test.WorkerThread(command, simulation_folder, form.flexaid_result_table, hex_colour_list, max_generations)
    time.sleep(1)
    worker.start()
    worker.table_signal_received.connect(receive_list)
    worker.current_generation_signal_received.connect(form.flexaid_progress.setValue)
    worker.generation_str_signal_received.connect(form.generation_label.setText)
    worker.finished.connect(worker.quit)
    worker.finished.connect(lambda: toggle_buttons(form, False))
    worker.finished.connect(lambda: load_show_flexaid_result(simulation_folder))
    worker.finished.connect(lambda: form.flexaid_progress.setValue(max_generations))
    worker.finished.connect(lambda: form.generation_label.setText(f'Generation: {max_generations}/{max_generations}'))
    # worker.finished.connect(lambda: load_show_cleft(cleft_save_path, color_list, form.output_box, pymol_object))


def pause_simulation(form):
    simulation_path = form.simulate_folder_path.text()
    with open(os.path.join(simulation_path, '.pause'), 'a'):
        pass
    form.flexaid_button_pause.setText('Resume')


def abort_simulation(form):
    simulation_path = form.simulate_folder_path.text()
    with open(os.path.join(simulation_path, '.abort'), 'a'):
        pass
    form.flexaid_button_start.setEnabled(True)
    form.flexaid_tab.setCurrentIndex(1)


def stop_simulation(form):
    simulation_path = form.simulate_folder_path.text()
    with open(os.path.join(simulation_path, '.stop'), 'a'):
        pass
    form.flexaid_button_start.setEnabled(True)


def resume_simulation(form):
    simulation_path = form.simulate_folder_path.text()
    os.remove(os.path.join(simulation_path, '.pause'))
    form.flexaid_button_pause.setText('Pause')


def load_show_flexaid_result(result_path):
    result_files = next(os.walk(result_path), (None, None, []))[2]
    result_files = sorted(result_files)
    cmd.disable('everything')
    for file in result_files:
        if file.startswith('RESULT') and not file.endswith('INI.pdb') and file.endswith('.pdb'):
            file_path = os.path.join(result_path, file)
            cmd.load(file_path)


def pause_resume_simulation(form):
    if form.flexaid_button_pause.text() == 'Pause':
        pause_simulation(form)
    elif form.flexaid_button_pause.text() == 'Resume':
        resume_simulation(form)


def run_flexaid(flexaid_output_path, form, cleft_save_path, process_ligand_path, flexaid_path, simulation_folder_path, hex_colour_list, tmp_path):
    if form.flexaid_button_start.text() == 'Start':
        max_results = 10
        setting_dictionary = get_simulation_settings(form)
        max_generations = int(setting_dictionary['number_generations'])
        form.flexaid_progress.setMaximum(max_results)
        date = datetime.datetime.now()
        date_time_str = date.strftime("%d-%m-%y-%I-%M-%S")
        flexaid_result_path = os.path.join(simulation_folder_path, date_time_str)
        form.simulate_folder_path.setText(flexaid_result_path)
        os.mkdir(flexaid_result_path)
        os.mkdir(os.path.join(flexaid_result_path, 'temp'))
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
        toggle_buttons(form, True)
        flexaid_command = f'"{flexaid_path}" "{config_file_path}" "{ga_path}" "{flexaid_result_name_path}"'
        # with open(os.path.join(tmp_path, 'flex_cmd.txt'), 'w') as f:
        #     f.write(flexaid_command)
        print(flexaid_command)
        form.output_box.append(f'Please wait...Running Flexaid with command: \n{flexaid_command}')
        form.flexaid_tab.setCurrentIndex(2)
        run_flexaid_worker(flexaid_command, form, flexaid_result_path, hex_colour_list, max_generations)

