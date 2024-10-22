from PyQt5.QtCore import pyqtSignal, QThread, QObject
import os
from pymol import cmd
from PyQt5.QtGui import QStandardItemModel, QStandardItem
import subprocess
import sys
import pandas as pd
import glob
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
# TODO: run on more than 1 bd site
# TODO: load own ligands (generate library from smiles)

def run_subprocess(current_ligand_number, last_ligand, install_dir, nrgdock_target_folder, ligand_path, config_path, nrgdock_output_path, ligand_number):
    """Function to run subprocess command."""
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
                    '-te', nrgdock_output_path],
                   check=True)

class WorkerThread(QThread):
    progress_signal = pyqtSignal(str)

    def __init__(self, ligand_set_folder_path, install_dir, nrgdock_output_path, nrgdock_top_n_poses,
                 nrgdock_start_ligand, ligand_set_name, target_name, binding_site_name):
        super().__init__()
        self.ligand_set_folder_path = ligand_set_folder_path
        self.install_dir = install_dir
        self.nrgdock_output_path = nrgdock_output_path
        self.nrgdock_top_n_poses = nrgdock_top_n_poses
        self.nrgdock_start_ligand = nrgdock_start_ligand
        self.ligand_set_name = ligand_set_name
        self.target_name = target_name
        self.binding_site_name = binding_site_name
        self.number_of_cores = multiprocessing.cpu_count()-2

    def run(self):
        nrgdock_instance = NRGDockManager(self.ligand_set_folder_path, self.install_dir, self.nrgdock_output_path,
                                          self.nrgdock_top_n_poses, self.nrgdock_start_ligand, self.ligand_set_name,
                                          self.target_name, self.binding_site_name, self.number_of_cores)
        nrgdock_instance.progress_signal.connect(self.forward_progress_signal)
        nrgdock_instance.run_nrgdock()

    def forward_progress_signal(self, message):
        self.progress_signal.emit(message)


class NRGDockManager(QObject):
    progress_signal = pyqtSignal(str)

    def __init__(self, ligand_set_folder_path, install_dir, nrgdock_output_path, nrgdock_top_n_poses,
                 nrgdock_start_ligand, ligand_set_name, target_name, binding_site_name, number_of_cores):
        super().__init__()
        self.ligand_set_folder_path = ligand_set_folder_path
        self.install_dir = install_dir
        self.nrgdock_output_path = nrgdock_output_path
        self.step = 100
        # to add to main function:
        self.nrgdock_top_n_poses = nrgdock_top_n_poses
        self.starting_ligand = nrgdock_start_ligand
        self.ligand_set_name = ligand_set_name
        self.target_name = target_name
        self.binding_site_name = binding_site_name
        self.number_of_cores = number_of_cores

    def run_nrgdock(self):
        deps_path = os.path.join(self.install_dir, 'deps', 'nrgdock')
        config_path = os.path.join(deps_path, 'config.txt')
        nrgdock_target_folder = os.path.join(self.nrgdock_output_path, 'target')
        if not os.path.exists(nrgdock_target_folder):
            os.mkdir(nrgdock_target_folder)
        ligand_path = os.path.join(self.ligand_set_folder_path, self.ligand_set_name, 'preprocessed_ligands_1_conf')
        total_number_ligands = len(np.load(os.path.join(ligand_path, 'ligand_atom_type.npy')))
        target_file_path = os.path.join(nrgdock_target_folder, 'receptor.mol2')
        cmd.save(target_file_path, self.target_name)
        binding_site_folder_path = os.path.join(nrgdock_target_folder, 'get_cleft')
        if not os.path.isdir(binding_site_folder_path):
            os.mkdir(binding_site_folder_path)
        binding_site_file_path = os.path.join(binding_site_folder_path, 'receptor_sph_1.pdb')
        cmd.save(binding_site_file_path, self.binding_site_name)
        nrgdock_result_folder = os.path.join(self.nrgdock_output_path, 'results')
        if not os.path.exists(nrgdock_result_folder):
            os.mkdir(nrgdock_result_folder)
        self.progress_signal.emit(f"NRGDock: Processing Target")
        subprocess.run([sys.executable, os.path.join(self.install_dir, 'src', 'nrgdock', 'process_target.py'),
                        '-p', self.nrgdock_output_path, '-t', 'target', '-o', '-d',
                        os.path.join(self.install_dir, 'deps', 'nrgdock')], check=True)
        self.progress_signal.emit(f"NRGDock: Screening has started")
        with ThreadPoolExecutor(self.number_of_cores) as executor:
            futures = []
            for current_ligand_number in range(self.starting_ligand, total_number_ligands, self.step):
                last_ligand = current_ligand_number + self.step
                # Submit the subprocess task to the executor for parallel execution
                futures.append(executor.submit(
                    run_subprocess,
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
                except Exception as e:
                    print(f"Error occurred: {e}")

        top_n_name_list, csv_output_path = merge_csv(os.path.join(nrgdock_result_folder, 'target'))
        manage_poses(top_n_name_list, os.path.join(self.nrgdock_output_path, 'ligand_poses', self.target_name))
        # update_nrgdock_result_table(csv_output_path, self.form)
        # self.form.NRGDock_tabs.setTabEnabled(2, True)


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


def manage_poses(top_n_name_list, ligand_poses_folder):
    ligand_files = glob.glob(os.path.join(ligand_poses_folder, "*.pdb"))
    for file in ligand_files:
        file_name = os.path.splitext(os.path.basename(file))[0]
        if file_name not in top_n_name_list:
            os.remove(file)
        else:
            cmd.load(file)


def update_nrgdock_result_table(csv_result, form):
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
    form.nrgdock_result_table.setModel(model)
    form.nrgdock_result_table.resizeColumnsToContents()
    form.NRGDock_tabs.setCurrentIndex(2)


