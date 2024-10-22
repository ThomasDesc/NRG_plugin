import os
from pymol import cmd
from PyQt5.QtGui import QStandardItemModel, QStandardItem
import subprocess
import sys
# from src.nrgdock.process_target import main as process_target
from src.nrgdock.main_processed_target import main as nrgdock_main
import pandas as pd
import glob
import numpy as np
# TODO: run on more than 1 bd site
# TODO: load own ligands (generate library from smiles)
# TODO: group results
# TODO load results in order

def process_ligands():
    print()


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


def get_nrgdock_result_model(csv_result, form):
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



step = 100

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


def run_nrgdock(form, nrgdock_output_path, ligand_set_folder_path, install_dir):
    deps_path = os.path.join(install_dir, 'deps', 'nrgdock')
    config_path = os.path.join(deps_path, 'config.txt')
    n_poses_to_save = form.nrgdock_top_n_poses.text()
    starting_ligand = int(form.nrgdock_start_ligand.text())
    print(n_poses_to_save)
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
        binding_site_file_path = os.path.join(binding_site_folder_path, 'receptor_sph_1.pdb')
        cmd.save(binding_site_file_path, binding_site_name)
    nrgdock_result_folder = os.path.join(nrgdock_output_path, 'results')
    if not os.path.exists(nrgdock_result_folder):
        os.mkdir(nrgdock_result_folder)
    # form.output_box.append("Processing target...")
    # form.output_box.repaint()
    # QApplication.processEvents()
    subprocess.run([sys.executable, os.path.join(install_dir, 'src', 'nrgdock', 'process_target.py'),
                      '-p', nrgdock_output_path, '-t', 'target', '-o', '-d', os.path.join(install_dir, 'deps', 'nrgdock')], check=True)
    # process_target(nrgdock_output_path, ['target'], overwrite=True, run_getcleft=False, deps_path=deps_path)
    # form.output_box.append("Running NRGDock...")
    # form.output_box.repaint()

    # To keep:
    # form.nrgdock_progress.setEnabled(True)
    # form.nrgdock_progress_label.setText(f'Molecules docked: 0/{ligand_number}')
    # form.nrgdock_progress_bar.setMaximum(ligand_number)

    # form.nrgdock_progress_label.repaint()
    # form.nrgdock_progress.repaint()
    # form.nrgdock_progress_bar.repaint()
    # QApplication.processEvents()
    step = 100
    import time
    from concurrent.futures import ThreadPoolExecutor, as_completed
    start_time = time.time()
    with ThreadPoolExecutor() as executor:
        futures = []
        for current_ligand_number in range(starting_ligand, ligand_number, step):
            last_ligand = current_ligand_number + step
            # Submit the subprocess task to the executor for parallel execution
            futures.append(executor.submit(
                run_subprocess,
                current_ligand_number,
                last_ligand,
                install_dir,
                nrgdock_target_folder,
                ligand_path,
                config_path,
                nrgdock_output_path,
                ligand_number
            ))

        # Optionally: Wait for all tasks to complete and handle any results or exceptions
        for future in as_completed(futures):
            try:
                future.result()  # This will raise any exception from the subprocess
            except Exception as e:
                print(f"Error occurred: {e}")

        # form.nrgdock_progress_label.setText(f'Ligand: {last_ligand}/{ligand_number}')
        # form.nrgdock_progress_bar.setValue(last_ligand)

        # form.nrgdock_progress_label.repaint()
        # form.nrgdock_progress_bar.repaint()
        # QApplication.processEvents()

    # form.output_box.append("Done NRGDock")

    # form.output_box.repaint()
    # QApplication.processEvents()

    end_time = time.time()
    execution_time = end_time - start_time  # Calculate time difference
    print(f"Execution time: {execution_time} seconds")
    top_n_name_list, csv_output_path = merge_csv(os.path.join(nrgdock_result_folder, 'target'))
    manage_poses(top_n_name_list, os.path.join(nrgdock_output_path, 'ligand_poses', target_name))

    # get_nrgdock_result_model(csv_output_path, form)

    form.NRGDock_tabs.setTabEnabled(2, True)

