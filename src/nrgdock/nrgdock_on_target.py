from tokenize import group

from PyQt5.QtCore import pyqtSignal, QThread, QObject
import os
import subprocess
import sys
import pandas as pd
import glob
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from PyQt5.QtGui import QStandardItemModel, QStandardItem, QMovie
from pymol import cmd
import multiprocessing
import general_functions
import re
# TODO: run on more than 1 bd site
# TODO: better folder management (store each run in a folder named by target, binding site, ligand set)
# TODO: do not preprocess if new run uses same target and binding site


# print objects in a group: cmd.get_object_list('(NRGDock_3cqw_sph_2)')
def get_group_of_object(object_name):
    all_objects = cmd.get_names("all")  # Get all objects and groups in PyMOL
    for pymol_object in all_objects:
        if pymol_object != 'NRGDock' and cmd.get_type(pymol_object) == "object:group":  # Check if it is a group
            group_members = cmd.get_object_list(f'({pymol_object})')
            if object_name in group_members:
                return pymol_object
    return None

def show_ligand_from_table(table_object, binding_site, ligand_set):
    selected_indexes = table_object.selectionModel().selectedIndexes()
    ligands_to_show = []
    for index in selected_indexes:
        cell_text = index.data()
        column = index.column()
        if column == 0:
            ligands_to_show.append(cell_text)
    if ligands_to_show:
        ligands_to_show = [f"{item}_{binding_site}_{ligand_set.replace(' ', '_')}" for item in ligands_to_show]
        group_name = get_group_of_object(ligands_to_show[0])
        if not group_name:
            print('no group found')
            return
        else:
            group_objects = cmd.get_object_list(f'({group_name})')
            for group_object in group_objects:
                if group_object not in ligands_to_show:
                    cmd.disable(group_object)
                else:
                    cmd.enable(group_object)


class NRGDockManager:
    def __init__(self, form, install_dir, ligand_set_folder_path, model):
        super().__init__()
        self.form = form
        self.install_dir = install_dir
        self.ligand_set_folder_path = ligand_set_folder_path
        self.model = model
        self.deps_path = os.path.join(self.install_dir, 'deps', 'nrgdock')
        self.config_path = os.path.join(self.deps_path, 'config.txt')
        self.update_config()

    def initialise_progress_bar(self):
        self.form.nrgdock_progress.show()
        self.form.nrgdock_progress.setEnabled(True)
        self.form.nrgdock_progress_label.setText('Screening progress: 0%')
        self.form.nrgdock_progress_bar.setValue(0)
        self.form.nrgdock_progress_bar.setMaximum(100)
        self.form.nrgdock_button_cancel.setEnabled(True)

    def start_loading_gif(self):
        self.label_size = self.form.nrgdock_loading_gif.size()
        self.movie = QMovie(os.path.join(self.install_dir, 'assets', 'loading.gif'))
        self.form.nrgdock_loading_gif.setMovie(self.movie)
        self.movie.setScaledSize(self.label_size)
        self.movie.start()
        self.form.nrgdock_loading_gif.show()

    def update_config(self):
        with open(self.config_path, 'r') as f:
            lines = f.readlines()
        for line_counter, line in enumerate(lines):
            if line.startswith('ROTATIONS_PER_AXIS'):
                number_of_rotations = self.form.nrgdock_ligand_rotations.text().strip()
                if not number_of_rotations.isdigit():
                    general_functions.output_message(self.form.output_box, 'Number of rotations must be an intiger. Defaulting to 9', 'warning')
                    number_of_rotations = 9
                elif not 1 <= int(number_of_rotations) <= 32:
                    general_functions.output_message(self.form.output_box, 'Number of rotations must be between 1 and 32 inclusively. Defaulting to 9', 'warning')
                    number_of_rotations = 9
                lines[line_counter] = f"{line.split(' ')[0]} {number_of_rotations}\n"
        with open(self.config_path, 'w') as out_file:
            out_file.writelines(lines)

    def run_nrgdock(self):
        nrgdock_output_path = os.path.join(self.form.temp_line_edit.text(), 'NRGDock')
        n_poses_to_save = self.form.nrgdock_top_n_poses.text()
        nrgdock_start_ligand = int(self.form.nrgdock_start_ligand.text())
        ligand_set_name = self.form.nrgdock_select_ligand.currentText().replace(' ', '_')
        target_name = self.form.nrgdock_select_target.currentText()
        cpu_usage_target = int(self.form.nrgdock_cpu_usage_target.currentText()[:-1])
        if cpu_usage_target == 100:
            number_of_cores = multiprocessing.cpu_count() + 4
        else:
            number_of_cores = round(multiprocessing.cpu_count() * (cpu_usage_target / 100))
        binding_site_name = self.form.nrgdock_select_binding_site.currentText()
        if target_name == '':
            general_functions.output_message(self.form.output_box, 'No target object selected', 'warning')
            return
        if binding_site_name == '':
            general_functions.output_message(self.form.output_box, 'No binding site object selected', 'warning')
            return

        self.initialise_progress_bar()
        self.start_loading_gif()
        general_functions.disable_run_mutate_buttons(self.form, disable=True)
        self.nrgdock_thread = NRGDockThread(self.ligand_set_folder_path, self.install_dir, nrgdock_output_path,
                                          n_poses_to_save, nrgdock_start_ligand, ligand_set_name,
                                          target_name, binding_site_name, number_of_cores)
        self.nrgdock_thread.message_signal.connect(self.handle_message_signal)
        self.nrgdock_thread.screen_progress_signal.connect(self.handle_screen_progress_signal)
        self.nrgdock_thread.update_table_signal.connect(self.update_nrgdock_result_table)
        self.nrgdock_thread.finished_signal.connect(self.handle_thread_finished)
        self.nrgdock_thread.start()

    def handle_message_signal(self, message):
        general_functions.output_message(self.form.output_box, message, 'valid')

    def handle_screen_progress_signal(self, value):
        current_value = self.form.nrgdock_progress_bar.value()
        if value > current_value:
            self.form.nrgdock_progress_bar.setValue(value)
            self.form.nrgdock_progress_label.setText(f'Screening progress: {value}%')

    def update_nrgdock_result_table(self, csv_result):
        self.model.clear()
        df = pd.read_csv(csv_result)
        df = df[['Name', 'CF']]
        self.model.setHorizontalHeaderLabels(df.columns.tolist())
        for index, row in df.iterrows():
            item_list = []
            for data in row:
                item = QStandardItem(str(data))
                item_list.append(item)
            self.model.appendRow(item_list)
        self.form.nrgdock_result_table.setModel(self.model)
        self.form.nrgdock_result_table.resizeColumnsToContents()
        self.form.NRGDock_tabs.setTabEnabled(2, True)
        self.form.NRGDock_tabs.setCurrentIndex(2)

    def handle_thread_finished(self):
        self.nrgdock_thread.stop()
        self.nrgdock_thread.quit()
        self.nrgdock_thread.wait()
        self.nrgdock_thread = None
        self.form.nrgdock_progress.hide()
        general_functions.disable_run_mutate_buttons(self.form, enable=True)
        self.movie.stop()
        self.form.nrgdock_loading_gif.hide()
        self.form.nrgdock_button_cancel.setDisabled(True)


class NRGDockThread(QThread):
    message_signal = pyqtSignal(str)
    screen_progress_signal = pyqtSignal(int)
    update_table_signal = pyqtSignal(str)
    finished_signal = pyqtSignal()

    def __init__(self, ligand_set_folder_path, install_dir, nrgdock_output_path, nrgdock_top_n_poses,
                 nrgdock_start_ligand, ligand_set_name, target_name, binding_site_name, number_of_cores):
        super().__init__()
        self.ligand_set_folder_path = ligand_set_folder_path
        self.install_dir = install_dir
        self.nrgdock_output_path = nrgdock_output_path
        self.step = 50
        self.nrgdock_top_n_poses = nrgdock_top_n_poses
        self.starting_ligand = nrgdock_start_ligand
        self.ligand_set_name = ligand_set_name
        self.target_name = target_name
        self.binding_site_name = binding_site_name
        self.number_of_cores = number_of_cores
        self.is_running = True
        self.folder_prep()

    def folder_prep(self):
        self.deps_path = os.path.join(self.install_dir, 'deps', 'nrgdock')
        self.config_path = os.path.join(self.deps_path, 'config.txt')
        self.nrgdock_target_folder = os.path.join(self.nrgdock_output_path, self.binding_site_name)
        self.ligand_poses_folder = os.path.join(self.nrgdock_output_path, 'ligand_poses', self.binding_site_name)
        self.docking_result_folder = os.path.join(self.nrgdock_output_path, 'results', self.binding_site_name)
        self.binding_site_folder_path = os.path.join(self.nrgdock_target_folder, 'get_cleft')
        os.makedirs(self.nrgdock_target_folder, exist_ok=True)
        os.makedirs(self.ligand_poses_folder, exist_ok=True)
        os.makedirs(self.docking_result_folder, exist_ok=True)
        os.makedirs(self.binding_site_folder_path, exist_ok=True)
        self.ligand_path = os.path.join(self.ligand_set_folder_path, self.ligand_set_name, 'preprocessed_ligands_1_conf')
        self.total_number_ligands = len(np.load(os.path.join(self.ligand_path, 'ligand_atom_type.npy')))
        self.target_file_path = os.path.join(self.nrgdock_target_folder, 'receptor.mol2')
        self.binding_site_file_path = os.path.join(self.binding_site_folder_path, 'receptor_sph_1.pdb')
        cmd.save(self.target_file_path, self.target_name)
        cmd.save(self.binding_site_file_path, self.binding_site_name)
        self.processes = []
        self.executor = None
        cmd.hide('everything', self.binding_site_name)
        cmd.show('mesh', self.binding_site_name)
        cmd.zoom(self.binding_site_name, buffer=3, complete=1)

    def stop(self):
        self.is_running = False

        if self.executor:
            self.executor.shutdown(wait=False)

        if self.process_target_process and self.process_target_process.poll() is None:
            self.process_target_process.terminate()

        for process in self.processes:
            if process.poll() is None:
                process.terminate()
                process.wait()


    @staticmethod
    def merge_csv(folder):
        csv_output_path = os.path.join(folder, "nrgdock_result.csv")
        csv_files = glob.glob(os.path.join(folder, "*.csv"))
        merged_df = pd.concat((pd.read_csv(file) for file in csv_files), ignore_index=True)
        filtered_df = merged_df[merged_df['CF'] != 100000000]
        sorted_df = filtered_df.sort_values(by='CF')
        sorted_df.to_csv(csv_output_path, index=False)
        for file in csv_files:
            os.remove(file)
        top_10_names = sorted_df['Name'].head(20).tolist()
        return top_10_names, csv_output_path

    @staticmethod
    def manage_poses(top_n_name_list, ligand_poses_folder, binding_site_name, target_name, ligand_set_name):
        ligand_files = glob.glob(os.path.join(ligand_poses_folder, "*.pdb"))
        top_n_name_set = set(top_n_name_list)
        for file in ligand_files:
            file_name = os.path.splitext(os.path.basename(file))[0]
            if file_name not in top_n_name_set:
                os.remove(file)
        nrg_result_group_name = f"nrg_{binding_site_name}_{ligand_set_name}"
        for name in top_n_name_list:
            pymol_object_name = f"{name}_{binding_site_name}_{ligand_set_name}"
            file_path = os.path.join(ligand_poses_folder, f"{name}.pdb")
            if os.path.exists(file_path):
                cmd.load(file_path, pymol_object_name)
                cmd.group(nrg_result_group_name, pymol_object_name)
        cmd.group('NRGDock', nrg_result_group_name)
        cmd.zoom(binding_site_name, buffer=4, complete=1)

    def run(self):
        self.message_signal.emit("=========== NRGDock ===========")
        self.message_signal.emit(f"Processing Target")
        self.process_target_process = subprocess.Popen([sys.executable,
                                                        os.path.join(self.install_dir, 'src', 'nrgdock', 'process_target.py'),
                                                        '-p', self.nrgdock_output_path,
                                                        '-t', self.binding_site_name,
                                                        '-o',
                                                        '-d', os.path.join(self.install_dir, 'deps', 'nrgdock')])
        while self.is_running and self.process_target_process.poll() is None:
            self.msleep(100)

        if self.is_running:
            self.message_signal.emit(f"Screening has started")
            completed_tasks = 0
            total_tasks = int(((self.total_number_ligands - self.starting_ligand) / self.step)) + 1
            self.commands = self.make_commands(self.nrgdock_target_folder, self.ligand_path, self.config_path, self.total_number_ligands)
            self.executor = ThreadPoolExecutor(max_workers=self.number_of_cores)
            futures = {self.executor.submit(self.run_command, command): command for command in self.commands}
            for future in as_completed(futures):
                if not self.is_running :
                    break
                try:
                    future.result()
                    completed_tasks += 1
                    progress_percentage = int((completed_tasks / total_tasks) * 100)
                    self.screen_progress_signal.emit(progress_percentage)
                except Exception as e:
                    print(f"Error occurred: {e}")

        if self.is_running:
            self.message_signal.emit('Screening has finished')
            top_n_name_list, csv_output_path = self.merge_csv(self.docking_result_folder)
            ligand_pose_path = os.path.join(self.nrgdock_output_path, 'ligand_poses',  self.binding_site_name)
            self.manage_poses(top_n_name_list, ligand_pose_path, self.binding_site_name, self.target_name, self.ligand_set_name)
            self.finished_signal.emit()
            self.message_signal.emit("=========== END NRGDock ===========")
            self.update_table_signal.emit(csv_output_path)

    def run_command(self, command):
        process = subprocess.Popen(command)
        self.processes.append(process)
        process.wait()


    def make_commands(self, nrgdock_target_folder, ligand_path, config_path, total_number_ligands):
        command_list = []
        for current_ligand_number in range(self.starting_ligand, total_number_ligands, self.step):
            last_ligand = min(current_ligand_number + self.step, total_number_ligands)
            temp_command = [sys.executable,
             os.path.join(self.install_dir, 'src', 'nrgdock', 'nrgdock.py'),
             '-p', nrgdock_target_folder,
             '-t', 'ligand',
             '-s', str(current_ligand_number),
             '-e', str(last_ligand),
             '-l', ligand_path,
             '-si', '2',
             '-c', config_path,
             '-te', self.nrgdock_output_path,
             '-mp']
            command_list.append(temp_command)
        return command_list


