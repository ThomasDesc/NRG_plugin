from PyQt5.QtCore import pyqtSignal, QThread, QObject
import os
from pymol import cmd
import subprocess
import sys
import pandas as pd
import glob
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from PyQt5.QtGui import QStandardItemModel, QStandardItem, QMovie
import multiprocessing
import general_functions
# TODO: run on more than 1 bd site
# TODO: load own ligands (generate library from smiles)

class NRGDockRunner:
    def __init__(self, form, install_dir, ligand_set_folder_path):
        super().__init__()
        self.form = form
        self.install_dir = install_dir
        self.ligand_set_folder_path = ligand_set_folder_path

    def initialise_progress_bar(self):
        self.form.nrgdock_progress.setEnabled(True)
        self.form.nrgdock_progress_label.setText('Screening progress: 0%')
        self.form.nrgdock_progress_bar.setValue(0)
        self.form.nrgdock_progress_bar.setMaximum(100)

    def start_loading_gif(self):
        self.label_size = self.form.nrgdock_loading_gif.size()
        self.movie = QMovie(os.path.join(self.install_dir, 'assets', 'loading.gif'))
        self.form.nrgdock_loading_gif.setMovie(self.movie)
        self.movie.setScaledSize(self.label_size)
        self.movie.start()

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
        self.nrgdock_thread = NRGDockManager(self.ligand_set_folder_path, self.install_dir, nrgdock_output_path,
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
        df = pd.read_csv(csv_result)
        df = df[['Name', 'CF']]
        model = QStandardItemModel()
        model.setHorizontalHeaderLabels(df.columns.tolist())
        for index, row in df.iterrows():
            item_list = []
            for data in row:
                item = QStandardItem(str(data))
                item_list.append(item)
            model.appendRow(item_list)
        self.form.nrgdock_result_table.setModel(model)
        self.form.nrgdock_result_table.resizeColumnsToContents()
        self.form.NRGDock_tabs.setTabEnabled(2, True)
        self.form.NRGDock_tabs.setCurrentIndex(2)

    def handle_thread_finished(self):
        general_functions.disable_run_mutate_buttons(self.form, enable=True)
        self.movie.stop()
        self.form.nrgdock_loading_gif.hide()


class NRGDockManager(QThread):
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
        self.step = 100
        self.nrgdock_top_n_poses = nrgdock_top_n_poses
        self.starting_ligand = nrgdock_start_ligand
        self.ligand_set_name = ligand_set_name
        self.target_name = target_name
        self.binding_site_name = binding_site_name
        self.number_of_cores = number_of_cores

    @staticmethod
    def run_subprocess(current_ligand_number, last_ligand, install_dir, nrgdock_target_folder, ligand_path, config_path,
                       nrgdock_output_path, ligand_number):
        if last_ligand > ligand_number:
            last_ligand = ligand_number
        # Run the subprocess
        subprocess.run([sys.executable,
                        os.path.join(install_dir, 'src', 'nrgdock', 'main_processed_target.py'),
                        '-p', nrgdock_target_folder,
                        '-t', 'ligand',
                        '-s', str(current_ligand_number),
                        '-e', str(last_ligand),
                        '-l', ligand_path,
                        '-si', '2',
                        '-c', config_path,
                        '-te', nrgdock_output_path,
                        '-mp'],
                       check=True)

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
    def manage_poses(top_n_name_list, ligand_poses_folder, binding_site_name):
        ligand_files = glob.glob(os.path.join(ligand_poses_folder, "*.pdb"))
        top_n_name_set = set(top_n_name_list)
        for file in ligand_files:
            file_name = os.path.splitext(os.path.basename(file))[0]
            if file_name not in top_n_name_set:
                os.remove(file)
        for name in top_n_name_list:
            file_path = os.path.join(ligand_poses_folder, f"{name}.pdb")
            if os.path.exists(file_path):
                cmd.load(file_path)
                cmd.group(f'NRGDock_{binding_site_name}', name)

    def run(self):
        deps_path = os.path.join(self.install_dir, 'deps', 'nrgdock')
        config_path = os.path.join(deps_path, 'config.txt')
        nrgdock_target_folder = os.path.join(self.nrgdock_output_path, self.binding_site_name)
        ligand_poses_folder = os.path.join(self.nrgdock_output_path, 'ligand_poses', self.binding_site_name)
        docking_result_folder = os.path.join(self.nrgdock_output_path, 'results', self.binding_site_name)
        os.makedirs(nrgdock_target_folder, exist_ok=True)
        os.makedirs(ligand_poses_folder, exist_ok=True)
        os.makedirs(docking_result_folder, exist_ok=True)

        ligand_path = os.path.join(self.ligand_set_folder_path, self.ligand_set_name, 'preprocessed_ligands_1_conf')
        total_number_ligands = len(np.load(os.path.join(ligand_path, 'ligand_atom_type.npy')))
        target_file_path = os.path.join(nrgdock_target_folder, 'receptor.mol2')
        cmd.save(target_file_path, self.target_name)
        binding_site_folder_path = os.path.join(nrgdock_target_folder, 'get_cleft')
        os.makedirs(binding_site_folder_path, exist_ok=True)
        binding_site_file_path = os.path.join(binding_site_folder_path, 'receptor_sph_1.pdb')
        cmd.save(binding_site_file_path, self.binding_site_name)
        cmd.hide('everything', self.binding_site_name)
        cmd.show('mesh', self.binding_site_name)

        self.message_signal.emit(f"=========== NRGDock ===========")
        self.message_signal.emit(f"NRGDock: Processing Target")
        subprocess.run([sys.executable, os.path.join(self.install_dir, 'src', 'nrgdock', 'process_target.py'),
                        '-p', self.nrgdock_output_path, '-t', self.binding_site_name, '-o', '-d',
                        os.path.join(self.install_dir, 'deps', 'nrgdock')], check=True)
        self.message_signal.emit(f"NRGDock: Screening has started")

        completed_tasks = 0
        total_tasks = int(((total_number_ligands - self.starting_ligand) / self.step)) + 1
        with ThreadPoolExecutor(self.number_of_cores) as executor:
            futures = []
            for current_ligand_number in range(self.starting_ligand, total_number_ligands, self.step):
                last_ligand = min(current_ligand_number + self.step, total_number_ligands)
                futures.append(executor.submit(
                    self.run_subprocess,
                    current_ligand_number,
                    last_ligand,
                    self.install_dir,
                    nrgdock_target_folder,
                    ligand_path,
                    config_path,
                    self.nrgdock_output_path,
                    total_number_ligands
                ))
            for future in as_completed(futures):
                try:
                    future.result()
                    completed_tasks += 1
                    progress_percentage = int((completed_tasks / total_tasks) * 100)
                    self.screen_progress_signal.emit(progress_percentage)
                except Exception as e:
                    print(f"Error occurred: {e}")
        self.message_signal.emit('NRGDock: Screening has finished')
        top_n_name_list, csv_output_path = self.merge_csv(docking_result_folder)
        self.manage_poses(top_n_name_list, os.path.join(self.nrgdock_output_path, 'ligand_poses',  self.binding_site_name), self.binding_site_name)
        self.finished_signal.emit()
        self.message_signal.emit(f"=========== END NRGDock ===========")
        self.update_table_signal.emit(csv_output_path)
