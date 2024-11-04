import os
import shutil
from general_functions import output_message


def delete_ligand_set(ligand_set_name, ligand_set_folder, output_box):
    full_ligand_set_name = ligand_set_name.replace(' ', '_')
    ligand_set_path = os.path.join(ligand_set_folder, full_ligand_set_name)
    print(ligand_set_path)
    if os.path.exists(ligand_set_path):
        shutil.rmtree(ligand_set_path)
        output_message(output_box, f'Deleted ligand set {ligand_set_name}', 'valid')

