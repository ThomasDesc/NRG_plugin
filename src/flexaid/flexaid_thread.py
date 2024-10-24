from PyQt5.QtCore import pyqtSignal, QThread
import general_functions
import subprocess
import os
from PyQt5.QtWidgets import QTableWidgetItem
from PyQt5.QtGui import QColor
from pymol import cmd
import datetime
import time
#TODO: abort and stop from flexaid.py

def load_color_list(color_list_path):
    with open(color_list_path, 'r') as file:
        color_list = [line.strip() for line in file]
    return color_list


class FlexAIDManager:
    def __init__(self, form, binary_folder_path, binary_suffix, install_dir):
        self.form = form
        self.flexaid_temp_path = os.path.join(self.form.temp_line_edit.text(), 'FlexAID')
        self.flexaid_simulation_folder = os.path.join(self.flexaid_temp_path, 'Simulation')
        os.makedirs(self.flexaid_temp_path, exist_ok=True)
        os.makedirs(self.flexaid_simulation_folder, exist_ok=True)
        self.binary_folder_path = binary_folder_path
        self.binary_suffix = binary_suffix
        self.install_dir = install_dir
        self.rmsd = form.flexaid_ligref_checkBox.isChecked()
        self.run_specific_simulate_folder_path = os.path.join(self.flexaid_simulation_folder, datetime.datetime.now().strftime("%d-%m-%y-%I-%M-%S"))
        self.update_file_path = os.path.join(self.run_specific_simulate_folder_path, ".update")
        self.simulation_settings = self.get_simulation_settings()
        self.max_generations = int(self.simulation_settings['number_generations'])
        self.form.flexaid_progress.setValue(0)
        self.form.flexaid_progress.setMaximum(self.max_generations)
        self.form.generation_label.setText(f'Generation: 0/{self.max_generations}')

        self.flexaid_deps_path = os.path.join(install_dir, 'deps', 'flexaid')
        self.hex_color_list = load_color_list(os.path.join(self.flexaid_deps_path, 'hex_colors.txt'))

    def start_run(self):
        self.start_file_updater()
        self.start_flexaid_thread()

    def start_file_updater(self):
        self.file_updater_thread = FileUpdaterThread(self.update_file_path, self.hex_color_list)
        self.file_updater_thread.update_signal.connect(self.perform_updates)
        self.file_updater_thread.start()

    def perform_updates(self, current_generation, data):
        self.form.flexaid_progress.setValue(current_generation)
        self.form.generation_label.setText(f'Generation: {current_generation}/{self.max_generations}')
        self.colour_specific_cell(data)

    def start_flexaid_thread(self):
        target_name = self.form.flexaid_select_target.currentText()
        target_save_path = os.path.join(self.flexaid_temp_path, 'flexaid_target.pdb')
        cmd.save(target_save_path, target_name)
        ligand_name = self.form.flexaid_select_ligand.currentText()
        ligand_save_path = os.path.join(self.flexaid_temp_path, 'flexaid_ligand.pdb')
        cmd.save(ligand_save_path, ligand_name)
        binding_site_name = self.form.flexaid_select_binding_site.currentText()
        binding_site_path = os.path.join(self.flexaid_temp_path, 'binding_site_sph_.pdb')
        cmd.save(binding_site_path, binding_site_name)
        self.toggle_buttons(True)
        general_functions.output_message(self.form.output_box, f"=========== FlexAID ===========", 'valid')
        general_functions.output_message(self.form.output_box, 'Running FlexAID...', 'valid')
        self.form.flexaid_tab.setCurrentIndex(2)
        self.form.flexaid_tab.setTabEnabled(2, True)
        self.worker = FlexAIDWorkerThread(self.binary_folder_path, self.binary_suffix, self.install_dir,
                                          self.run_specific_simulate_folder_path, self.flexaid_temp_path,
                                          target_save_path, ligand_save_path, binding_site_path, self.simulation_settings)
        self.worker.finished.connect(self.flexaid_end_processes)
        self.worker.start()

    def get_simulation_settings(self):
        setting_dictionary = {'number_chromosomes': self.form.input_num_chromosomes.text(),
                              'number_generations': self.form.input_num_generations.text()}
        return setting_dictionary

    def colour_specific_cell(self, data):
        table_widget = self.form.flexaid_result_table
        num_column = table_widget.columnCount()
        row = int(data[1]) - 1
        for column_counter, column in enumerate(range(num_column)):
            if column_counter == 0:
                item = QTableWidgetItem()
                item.setBackground(QColor(data[column_counter]))
            else:
                item = QTableWidgetItem(str(data[column_counter]))
            table_widget.setItem(row, column, item)

    def flexaid_end_processes(self):
        general_functions.output_message(self.form.output_box, f"=========== FlexAID ===========", 'valid')
        self.worker.quit()
        self.file_updater_thread.quit()
        self.toggle_buttons(False)
        self.load_show_flexaid_result()
        self.form.flexaid_progress.setValue(self.max_generations)
        self.form.generation_label.setText(f'Generation: {self.max_generations}/{self.max_generations}')
        if self.rmsd:
            self.show_rmsd()

    def toggle_buttons(self, true_false):
        self.form.flexaid_button_pause.setEnabled(true_false)
        self.form.flexaid_button_abort.setEnabled(true_false)
        self.form.flexaid_button_stop.setEnabled(true_false)

    def load_show_flexaid_result(self):
        result_files = next(os.walk(self.run_specific_simulate_folder_path), (None, None, []))[2]
        result_files = sorted(result_files)
        cmd.disable('everything')
        for file in result_files:
            if file.startswith('RESULT') and not file.endswith('INI.pdb') and file.endswith('.pdb'):
                file_path = os.path.join(self.run_specific_simulate_folder_path, file)
                cmd.load(file_path)
                cmd.group('flexaid_results', os.path.basename(file_path[:-4]))

    def show_rmsd(self):
        table = self.form.flexaid_result_table
        last_column = table.columnCount() - 1
        rmsd_list = []
        for i in range(5):
            with open(os.path.join(self.run_specific_simulate_folder_path, f'RESULT_{i}.pdb'), 'r') as t1:
                lines = t1.readlines()
                for line in lines:
                    if 'RMSD to' in line:
                        rmsd_list.append(line.split()[1])
                        break
        for row in range(5):
            table.setItem(row, last_column, QTableWidgetItem(rmsd_list[row]))


class FileUpdaterThread(QThread):
    update_signal = pyqtSignal(int, tuple)

    def __init__(self, update_file_path, hex_colour_list):
        super(FileUpdaterThread, self).__init__()
        self.update_file_path = update_file_path
        self.hex_colour_list = hex_colour_list
        self.current_generation = 0
        self.running = True

    def run(self):
        while self.running:
            if os.path.exists(self.update_file_path):
                self.read_update(self.hex_colour_list)
                os.remove(self.update_file_path)
                self.current_generation += 1
            time.sleep(0.1)

    def read_update(self, hex_colour_list, num_results=5):
        number_color_list = general_functions.create_number_list(num_results, len(hex_colour_list))
        with open(self.update_file_path, "r") as f:
            for line_counter, line in enumerate(f):
                if line_counter > 1 and self.current_generation % 10 == 0:
                    line = line.split()
                    top_number = int(line[0]) + 1
                    cf = line[-5]
                    fitness = line[-1]
                    rmsd = 'N/A'
                    data = (hex_colour_list[number_color_list[top_number - 1]], top_number, cf, fitness, rmsd)
                    self.update_signal.emit(self.current_generation, data)

    def stop(self):
        self.running = False


class FlexAIDWorkerThread(QThread):
    finished = pyqtSignal()

    def __init__(self, binary_folder_path, binary_suffix, install_dir, flexaid_result_path, flexaid_temp_path,
                 target_save_path, ligand_save_path, binding_site_path, setting_dictionary):
        super(FlexAIDWorkerThread, self).__init__()
        self.binary_folder_path = binary_folder_path
        self.binary_suffix = binary_suffix
        self.install_dir = install_dir
        self.max_results = 10
        self.flexaid_result_path = flexaid_result_path
        self.flexaid_temp_path = flexaid_temp_path
        self.target_save_path = target_save_path
        self.ligand_save_path = ligand_save_path
        self.binding_site_path = binding_site_path
        self.setting_dictionary = setting_dictionary

    @staticmethod
    def process_ligand(process_ligand_path, input_path, istarget=False):
        if istarget:
            process_ligand_command = f'"{process_ligand_path}" -f "{input_path}" -target'
        else:
            process_ligand_command = f'"{process_ligand_path}" -f "{input_path}" --atom_index 90000 -ref'
        print(process_ligand_command)
        subprocess.run(process_ligand_command, shell=True)

    @staticmethod
    def count_flex(ligand_inp_file_path):
        with open(ligand_inp_file_path, 'r') as t:
            texto = t.readlines()
            count = 0
            for line in texto:
                if line[0] == 'F':
                    count += 1
        return count

    def write_config(self, target_inp_path, cleft, ligand_inp_path, max_results, flexaid_output_path,
                     flexaid_result_path, flexaid_deps_path):
        with open(os.path.join(flexaid_deps_path, 'template_config.inp'), "r") as t1:
            lines = t1.readlines()
        config_file_output_path = os.path.join(flexaid_output_path, 'config.inp')
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
                    for flex in range(self.count_flex(ligand_inp_path)):
                        t2.write(f'OPTIMZ 9999 - {str(flex + 1)}\n')
                elif line.startswith('STATEP'):
                    t2.write(f'STATEP {flexaid_result_path}\n')
                elif line.startswith('RMSDST'):
                    t2.write(f'RMSDST {ligand_inp_path[:-4]}_ref.pdb\n')
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

    @staticmethod
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

    def run(self):
        flexaid_binary_path = os.path.join(self.binary_folder_path, f'FlexAID{self.binary_suffix}')
        process_ligand_path = os.path.join(self.binary_folder_path, f'Process_ligand{self.binary_suffix}')
        flexaid_deps_path = os.path.join(self.install_dir, 'deps', 'flexaid')
        os.mkdir(self.flexaid_result_path)
        os.mkdir(os.path.join(self.flexaid_result_path, 'temp'))
        flexaid_result_name_path = os.path.join(self.flexaid_result_path, "RESULT").replace('\\', '/')

        self.process_ligand(process_ligand_path, self.target_save_path, istarget=True)
        self.process_ligand(process_ligand_path, self.ligand_save_path)
        target_inp_path = os.path.splitext(self.target_save_path)[0] + '.inp.pdb'
        ligand_inp_path = os.path.splitext(self.ligand_save_path)[0] + '.inp'
        config_file_path = self.write_config(target_inp_path, self.binding_site_path, ligand_inp_path, self.max_results,
                                             self.flexaid_temp_path, self.flexaid_result_path,
                                             flexaid_deps_path).replace('\\', '/')

        ga_path = os.path.join(self.flexaid_temp_path, 'ga_inp.dat').replace('\\', '/')
        self.edit_ga(os.path.join(flexaid_deps_path, 'ga_inp.dat'), ga_path, self.setting_dictionary)
        flexaid_command = f'{flexaid_binary_path} "{config_file_path}" "{ga_path}" "{flexaid_result_name_path}"'
        print(flexaid_command)
        os.system(flexaid_command)
        self.finished.emit()
