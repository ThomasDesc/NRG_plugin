import os
from pymol import cmd
import sys
import subprocess
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nrgdock', 'src'))
try:
    import scipy
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'scipy'])
try:
    import numba
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'numba'])
from process_target import main as process_target
from main_processed_target import main as nrgdock_main

# TODO: run on more than 1 bd site
# TODO: load own ligands (generate library from smiles)


def process_ligands():


def run_nrgdock(form, nrgdock_output_path, ligand_set_folder_path, main_folder_path):
    config_path = os.path.join(main_folder_path, 'nrgdock', 'deps', 'config.txt')
    nrgdock_target_folder = os.path.join(nrgdock_output_path, 'target')
    if not os.path.exists(nrgdock_target_folder):
        os.mkdir(nrgdock_target_folder)
    drugbank_molecule_path = os.path.join(ligand_set_folder_path, 'drugbank_subset_processed')
    target_name = form.nrgdock_select_target.currentText()
    if target_name == '':
        print('No target object selected')
        return
    else:
        target_file_path = os.path.join(nrgdock_target_folder, target_name + '.mol2')
        cmd.save(target_file_path, target_name)

    binding_site_name = form.nrgdock_select_binding_site.currentText()
    if target_name == '':
        print('No binding site object selected')
        return
    else:
        binding_site_folder_path = os.path.join(nrgdock_target_folder, 'get_cleft')
        os.mkdir(binding_site_folder_path)
        binding_site_file_path = os.path.join(binding_site_folder_path, 'binding_site_sph_1.pdb')
        cmd.save(binding_site_file_path, binding_site_name)

    print('Nrgdock output path: ', nrgdock_target_folder)
    print('Nrgdock target path: ', target_name)
    process_target(nrgdock_output_path, ['target'], overwrite=True, run_getcleft=False)
    nrgdock_main(config_path, nrgdock_target_folder, 'ligand', 0, 10, target_name, None, None, drugbank_molecule_path, 2)
