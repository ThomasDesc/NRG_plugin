import os
from pymol import cmd
import sys
import subprocess
from PyQt5.QtWidgets import QApplication
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nrgdock', 'src'))
try:
    import scipy
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'scipy'])
try:
    import numba
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'numba'])
try:
    import pandas
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'pandas'])
from process_target import main as process_target
from main_processed_target import main as nrgdock_main
import pandas as pd
import glob
import numpy as np
# TODO: run on more than 1 bd site
# TODO: load own ligands (generate library from smiles)


def process_ligands():
    print()


def merge_csv(folder):
    csv_files = glob.glob(os.path.join(folder, "*.csv"))
    merged_df = pd.concat((pd.read_csv(file) for file in csv_files), ignore_index=True)
    filtered_df = merged_df[merged_df['CF'] != 100000000]
    sorted_df = filtered_df.sort_values(by='CF')
    sorted_df.to_csv(os.path.join(folder, "nrgdock_result.csv"), index=False)
    for file in csv_files:
        os.remove(file)
    top_10_names = sorted_df['Name'].head(10).tolist()
    return top_10_names


def manage_poses(top_n_name_list, ligand_poses_folder):
    ligand_files = glob.glob(os.path.join(ligand_poses_folder, "*.pdb"))
    for file in ligand_files:
        file_name = os.path.splitext(os.path.basename(file))[0]
        if file_name not in top_n_name_list:
            os.remove(file)
            print(f"Deleted file: {file}")
        else:
            cmd.load(file)


def run_nrgdock(form, nrgdock_output_path, ligand_set_folder_path, main_folder_path):

    config_path = os.path.join(main_folder_path, 'nrgdock', 'deps', 'config.txt')
    nrgdock_target_folder = os.path.join(nrgdock_output_path, 'target')
    if not os.path.exists(nrgdock_target_folder):
        os.mkdir(nrgdock_target_folder)
    ligand_path = os.path.join(ligand_set_folder_path, form.nrgdock_select_ligand.currentText().replace(' ', '_'), 'preprocessed_ligands_1_conf')
    ligand_number = len(np.load(os.path.join(ligand_path, 'ligand_atom_type.npy')))
    target_name = form.nrgdock_select_target.currentText()
    if target_name == '':
        print('No target object selected')
        return
    else:
        target_file_path = os.path.join(nrgdock_target_folder, 'receptor.mol2')
        cmd.save(target_file_path, target_name)

    binding_site_name = form.nrgdock_select_binding_site.currentText()
    if binding_site_name == '':
        print('No binding site object selected')
        return
    else:
        binding_site_folder_path = os.path.join(nrgdock_target_folder, 'get_cleft')
        if not os.path.isdir(binding_site_folder_path):
            os.mkdir(binding_site_folder_path)
        binding_site_file_path = os.path.join(binding_site_folder_path, binding_site_name + '.pdb')
        cmd.save(binding_site_file_path, binding_site_name)
    nrgdock_result_folder = os.path.join(nrgdock_output_path, 'results')
    if not os.path.exists(nrgdock_result_folder):
        os.mkdir(nrgdock_result_folder)
    form.output_box.append("Processing target...")
    form.output_box.repaint()
    QApplication.processEvents()
    process_target(nrgdock_output_path, ['target'], overwrite=True, run_getcleft=False)
    form.nrgdock_progress.setEnabled(True)
    form.nrgdock_progress_label.setText(f'Generation: 0/{ligand_number}')
    form.nrgdock_progress_bar.setMaximum(ligand_number)
    form.nrgdock_progress_label.repaint()
    form.nrgdock_progress.repaint()
    form.nrgdock_progress_bar.repaint()
    QApplication.processEvents()
    step = 100
    for current_ligand_number in range(0, ligand_number, step):
        last_ligand = current_ligand_number+step
        print(last_ligand)
        if last_ligand > ligand_number:
            last_ligand = ligand_number
        nrgdock_main(config_path, nrgdock_target_folder, 'ligand', current_ligand_number, last_ligand, target_name, None, None, ligand_path, 2, temp_path=nrgdock_output_path)
        form.nrgdock_progress_label.setText(f'Ligand: {last_ligand}/{ligand_number}')
        form.nrgdock_progress_bar.setValue(last_ligand)
        form.nrgdock_progress_label.repaint()
        form.nrgdock_progress_bar.repaint()
        QApplication.processEvents()
    form.output_box.append("Done NRGDock")
    form.output_box.repaint()
    QApplication.processEvents()
    top_n_name_list = merge_csv(os.path.join(nrgdock_result_folder, target_name))
    manage_poses(top_n_name_list, os.path.join(nrgdock_output_path, 'ligand_poses', target_name))