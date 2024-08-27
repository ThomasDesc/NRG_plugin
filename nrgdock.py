import os
from pymol import cmd
# TODO: run on more than 1 bd site
# TODO: load own ligands (generate library from smiles)


def run_nrgdock(form, nrgdock_output_path):
    target_name = form.nrgdock_select_target.currentText()
    if target_name == '':
        print('No object selected')
        return
    else:
        target_file_path = os.path.join(nrgdock_output_path, target_name + '.mol2')
        cmd.save(target_file_path, target_name)

    binding_site_name = form.nrgdock_select_binding_site.currentText()
    if target_name == '':
        print('No object selected')
        return
    else:
        binding_site_folder_path = os.path.join(nrgdock_output_path, 'get_cleft')
        os.mkdir(binding_site_folder_path)
        binding_site_file_path = os.path.join(binding_site_folder_path, 'binding_site_sph_1.pdb')
        cmd.save(binding_site_file_path, binding_site_name)

