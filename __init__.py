from __future__ import absolute_import
from __future__ import print_function


def __init_plugin__(app):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('NRGSuite_Qt', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()
    dialog.raise_()


def make_dialog():
    from pymol.Qt import QtGui, QtWidgets
    from pymol.Qt.utils import loadUi

    import os
    import sys
    import shutil
    import subprocess

    install_dir = os.path.dirname(__file__)
    sys.path.append(install_dir)
    import getcleft
    import general_functions
    import flexaid
    import platform
    import spheres
    import run_Surfaces
    import thread_test
    import nrgdock
    dialog = QtWidgets.QDialog()

    OS = platform.system().upper()
    suffix = ''
    # SET the configuration location
    if OS == 'LINUX' or OS == 'BSD':
        folder = 'linux'
    elif OS == 'DARWIN':
        folder = 'mac'
    elif OS == 'WINDOWS' or OS == 'MICROSOFT' or OS == 'WIN32':
        folder = 'win'
        suffix = '.exe'
    else:
        exit('Unknown OS')

    uifile = os.path.join(install_dir, 'nrgdock_widget.ui')
    getcleft_path = os.path.join(install_dir, 'bin', folder, f'GetCleft{suffix}')
    flexaid_path = os.path.join(install_dir, 'bin', folder, f'FlexAID{suffix}')
    vcon_path = os.path.join(install_dir, 'bin', folder, f'vcon')
    process_ligand_path = os.path.join(install_dir, 'bin', folder, f'Process_ligand{suffix}')
    bin_list = [getcleft_path, flexaid_path, vcon_path, process_ligand_path]
    print("Binaries are stored at: ", os.path.join(install_dir, 'bin', folder))
    if folder == 'mac':
        for file in bin_list:
            subprocess.run(["chmod", "755", file])
            result = subprocess.run([file], capture_output=True, text=True)
            if result.returncode != 0 and result.returncode != -11:
                print('Could not run: ', file)
    ligand_set_folder_path = os.path.join(install_dir, 'nrgdock', 'ligand_sets')
    plugin_tmp_output_path = os.path.join(os.path.expanduser('~'), 'Documents', 'NRGSuite_Qt')
    temp_path = os.path.join(plugin_tmp_output_path, 'temp')
    getcleft_output_path = os.path.join(temp_path, 'GetCleft')
    flexaid_output_path = os.path.join(temp_path, 'FlexAID')
    nrgdock_output_path = os.path.join(temp_path, 'NRGDock')
    surfaces_output_path = os.path.join(temp_path, 'Surfaces')
    cleft_save_path = os.path.join(getcleft_output_path, 'Clefts')
    simulation_folder_path = os.path.join(flexaid_output_path, 'Simulation')
    color_list = ['red', 'br8', 'tv_red', 'oxygen', 'iron', 'tv_orange', 'sulfur', 'gold', 'yelloworange', 'neodymium',
                  'limon', 'chartreuse', 'tv_green', 'limegreen', 'teal', 'rhodium', 'slate', 'tv_blue', 'blue',
                  'density']
    hex_colour_list = ['#FF0000', '#E61A33', '#FF1515', '#FF4D4D', '#E06633', '#FF8C26', '#E6C740', '#FFD124',
                       '#FFDE5E', '#C7FFC7', '#BFFF40', '#80FF00', '#33FF33', '#00FF80', '#00BFBF', '#0A7DAB',
                       '#8080FF', '#4D4DFF', '#0000FF', '#1A1A99']

    if os.path.isdir(plugin_tmp_output_path):
        shutil.rmtree(plugin_tmp_output_path)
    os.mkdir(plugin_tmp_output_path)
    os.mkdir(temp_path)
    os.mkdir(getcleft_output_path)
    os.mkdir(cleft_save_path)
    os.mkdir(surfaces_output_path)
    os.mkdir(flexaid_output_path)
    os.mkdir(simulation_folder_path)
    os.mkdir(nrgdock_output_path)
    form = loadUi(uifile, dialog)
    form.stackedWidget.setCurrentIndex(0)
    form.flexaid_tab.setTabEnabled(2, False)

    general_functions.refresh_dropdown(form.cleft_select_object, form.output_box, no_warning=True)
    general_functions.refresh_folder(ligand_set_folder_path, form.nrgdock_select_ligand)
    form.button_getcleft.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(0))
    form.button_partition_cleft.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(1))
    form.button_flexaid.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(2))
    form.button_nrgdock.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(3))
    form.button_Surfaces.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(4))

    form.cleft_partition_button_add.clicked.connect(lambda: spheres.display_sphere(form.cleft_partition_select_object.currentText(), form.cleft_partition_radius_slider, form.partition_sphere_select))

    form.button_hide.clicked.connect(lambda: general_functions.pymol_hide_structures(form))
    form.cleft_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_select_object, form.output_box))
    form.cleft_partition_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_partition_select_object, form.output_box, filter_for='_sph'))
    form.cleft_partition_button_move.clicked.connect(lambda: spheres.move_sphere(form.partition_sphere_select.currentText()))

    form.button_start.clicked.connect(lambda: getcleft.run_getcleft(form, getcleft_path, getcleft_output_path,
                                                                    cleft_save_path, color_list))

    # FlexAID:
    form.flexaid_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_target, form.output_box, exclude='_sph'))
    form.flexaid_ligand_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_ligand, form.output_box, exclude='_sph'))
    form.flexaid_binding_site_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_binding_site, form.output_box, filter_for='_sph'))
    form.flexaid_button_start.clicked.connect(lambda: form.flexaid_tab.setTabEnabled(2, True))
    form.flexaid_retrieve_nrgdock_ligands.clicked.connect(lambda: flexaid.retrieve_nrgdock_ligands(nrgdock_output_path))
    form.flexaid_button_start.clicked.connect(lambda: flexaid.run_flexaid(flexaid_output_path, form, process_ligand_path, flexaid_path, simulation_folder_path, hex_colour_list))

    form.flexaid_button_pause.clicked.connect(lambda: flexaid.pause_resume_simulation(form))
    form.flexaid_button_stop.clicked.connect(lambda: flexaid.stop_simulation(form))
    form.flexaid_button_abort.clicked.connect(lambda: flexaid.abort_simulation(form))

    form.cleft_partition_radius_slider.valueChanged.connect(lambda: spheres.resize_sphere(form.partition_sphere_select.currentText(), form.cleft_partition_radius_slider.value()))
    form.cleft_partition_crop_button.clicked.connect(lambda: spheres.crop_cleft(form.partition_sphere_select.currentText(), form.cleft_partition_radius_slider.value()/100, cleft_save_path, form.cleft_partition_select_object.currentText()))

    # NRGDock:
    form.nrgdock_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.nrgdock_select_target, form.output_box, exclude='_sph'))
    form.nrgdock_binding_site_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.nrgdock_select_binding_site, form.output_box, filter_for='_sph'))
    form.nrgdock_delete_ligand_set_refresh.clicked.connect(lambda: general_functions.refresh_folder(ligand_set_folder_path, form.nrgdock_delete_ligand_set_dropdown))
    form.nrgdock_add_ligandset_button.clicked.connect(lambda: general_functions.folder_browser(form.nrgdock_add_ligand_file_path, ligand_set_folder_path, "Smiles Files (*.smi)"))
    form.nrgdock_button_ligandset_add.clicked.connect(lambda: nrgdock.process_ligands())
    form.nrgdock_ligand_set_refresh.clicked.connect(lambda: general_functions.refresh_folder(ligand_set_folder_path, form.nrgdock_select_ligand))
    form.nrgdock_button_start.clicked.connect(
        lambda: nrgdock.run_nrgdock(form, nrgdock_output_path, ligand_set_folder_path, install_dir))
    form.nrgdock_result_browse_button.clicked.connect(lambda: general_functions.folder_browser(form.nrgdock_result_path, nrgdock_output_path, "CSV file (*.csv)"))

    # Surfaces:
    form.surfaces_refresh_button.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_result, form.output_box, filter_for='RESULT'))
    form.surfaces_run_button.clicked.connect(lambda: run_Surfaces.run_run_surfaces(form.surface_select_result.currentText(), surfaces_output_path, form.simulate_folder_path.text(), install_dir, vcon_path))
    form.surfaces_run_button.clicked.connect(lambda: general_functions.surfaces_enable_buttons(form))
    form.surfaces_retreive_flexaid_result.clicked.connect(lambda: run_Surfaces.retrieve_flexaid_result(form.simulate_folder_path.text()))
    form.surfaces_retreive_flexaid_result.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_result, form.output_box, filter_for='RESULT'))
    # form.class_test.clicked.connect(lambda: getcleft.test_submit_command())
    return dialog
