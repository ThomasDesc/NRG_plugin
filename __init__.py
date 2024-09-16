# from __future__ import absolute_import
# from __future__ import print_function
import os
import sys
import shutil
import subprocess
import importlib.metadata


def __init_plugin__(app):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('NRGSuite_Qt', run_plugin_gui)


dialog = None


def run_plugin_gui():
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()
    dialog.raise_()


def test_binary(binary_folder_path, operating_system):
    all_files = os.listdir(binary_folder_path)
    binary_files = []
    for f in all_files:
        full_path = os.path.join(binary_folder_path, f)
        if os.path.isfile(full_path) and not f.startswith('.'):
            binary_files.append(full_path)
    if operating_system == 'mac':
        for file in binary_files:
            subprocess.run(["chmod", "755", file])
            result = subprocess.run([file], capture_output=True, text=True)
            if result.returncode != 0 and result.returncode != -11:
                print('Could not run: ', file)


def install_package(package, main_folder_path):
    try:
        __import__(package)
        print(f"{package} is already installed.")
    except ImportError:
        if package == 'modeller':
            print('Modeller install not detected. Please install via conda. For now the modeller tab will be unavailable')
        else:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
            if package == 'nrgten':
                distribution = importlib.metadata.distribution(package)
                shutil.copy(os.path.join(main_folder_path, 'deps', 'nrgten', 'amino_acids.atomtypes'),
                            os.path.join(str(distribution.locate_file('')), 'nrgten' 'config' 'amino_acids.atomtypes'))
                shutil.copy(os.path.join(main_folder_path, 'deps', 'nrgten', 'amino_acids.masses'),
                            os.path.join(str(distribution.locate_file('')), 'nrgten' 'config' 'amino_acids.masses'))


def make_dialog():
    from pymol.Qt import QtGui, QtWidgets
    from pymol.Qt.utils import loadUi

    install_dir = os.path.dirname(__file__)
    sys.path.append(install_dir)
    packages = ['nrgten', 'biopython', 'pandas', 'matplotlib', 'colour', 'scipy', 'numpy=2.0', 'numba', 'pandas']
    install_package('nrgten', install_dir)
    install_package('biopython', install_dir)
    install_package('pandas', install_dir)
    install_package('matplotlib', install_dir)
    install_package('colour', install_dir)
    from src.flexaid import flexaid
    from src.getcleft import getcleft
    from src.surfaces import surfaces
    from src.nrgdock import nrgdock
    from src.getcleft import spheres
    import general_functions
    from src.surfaces import run_Surfaces

    import platform
    # import spheres
    # import flexaid_thread
    dialog = QtWidgets.QDialog()

    OS = platform.system().upper()
    binary_suffix = ''
    if OS == 'LINUX' or OS == 'BSD':
        operating_system = 'linux'
    elif OS == 'DARWIN':
        operating_system = 'mac'
    elif OS == 'WINDOWS' or OS == 'MICROSOFT' or OS == 'WIN32':
        operating_system = 'win'
        binary_suffix = '.exe'
    else:
        exit('Unknown operating system')

    uifile = os.path.join(install_dir, 'nrgdock_widget.ui')
    binary_folder_path = os.path.join(install_dir, 'bin', operating_system)
    print('binary path: ', binary_folder_path)
    test_binary(binary_folder_path, operating_system)

    ligand_set_folder_path = os.path.join(install_dir, 'nrgdock_ligand_sets')
    plugin_tmp_output_path = os.path.join(os.path.expanduser('~'), 'Documents', 'NRGSuite_Qt')
    temp_path = os.path.join(plugin_tmp_output_path, 'temp')
    nrgdock_output_path = os.path.join(temp_path, 'NRGDock')
    surfaces_output_path = os.path.join(temp_path, 'Surfaces')
    modeller_save_path = os.path.join(temp_path, 'modeller')

    if os.path.isdir(plugin_tmp_output_path):
        shutil.rmtree(plugin_tmp_output_path)
    os.mkdir(plugin_tmp_output_path)
    os.mkdir(temp_path)
    os.mkdir(surfaces_output_path)
    os.mkdir(nrgdock_output_path)
    os.mkdir(modeller_save_path)
    form = loadUi(uifile, dialog)
    try:
        import modeller
    except ModuleNotFoundError:
        general_functions.output_message(form.output_box, 'Modeller install not detected. '
                                                          'Please install via conda. For now the modeller tab will be '
                                                          'unavailable', 'warning')
        form.button_nrgten.setEnabled(False)
        form.button_modeller.setEnabled(False)
        form.button_nrgten.setStyleSheet("background-color: black; color: white;")
        form.button_modeller.setStyleSheet("background-color: black; color: white;")
        # form.button_nrgten.hide()
        # form.button_modeller.hide()
    else:
        from src.nrgten import run_NRGTEN
        from src.modeller import run_modeller
    form.stackedWidget.setCurrentIndex(0)
    form.flexaid_tab.setTabEnabled(2, False)
    form.NRGDock_settings.setTabEnabled(2, False)
    if operating_system == 'mac':
        form.flexaid_multithread_button.setChecked(True)
    print(form.surface_select_result.currentText())

    general_functions.refresh_dropdown(form.cleft_select_object, form.output_box, no_warning=True)
    general_functions.refresh_folder(ligand_set_folder_path, form.nrgdock_select_ligand)
    form.button_getcleft.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(0))
    form.button_flexaid.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(1))
    form.button_nrgdock.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(2))
    form.button_nrgten.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(3))
    form.button_surfaces.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(4))
    form.button_modeller.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(5))

    # GetCleft
    form.button_hide.clicked.connect(lambda: general_functions.pymol_hide_structures(form))
    form.cleft_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_select_object, form.output_box))
    form.button_start.clicked.connect(lambda: getcleft.run_getcleft(form, binary_folder_path, binary_suffix, temp_path,
                                                                    install_dir))
    # Partition Cleft
    form.cleft_partition_button_add.clicked.connect(lambda: spheres.display_sphere(form.cleft_partition_select_object.currentText(), form.cleft_partition_radius_slider, form.partition_sphere_select, temp_path))
    form.cleft_partition_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_partition_select_object, form.output_box, filter_for='_sph'))
    form.cleft_partition_button_move.clicked.connect(lambda: spheres.move_sphere(form.partition_sphere_select.currentText()))
    form.cleft_partition_radius_slider.valueChanged.connect(lambda: spheres.resize_sphere(form.partition_sphere_select.currentText(), form.cleft_partition_radius_slider.value()))
    form.cleft_partition_crop_button.clicked.connect(lambda: spheres.crop_cleft(form.partition_sphere_select.currentText(), form.cleft_partition_radius_slider.value()/100, temp_path, form.cleft_partition_select_object.currentText()))

    # FlexAID:
    form.flexaid_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_target, form.output_box, exclude='_sph'))
    form.flexaid_ligand_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_ligand, form.output_box, exclude='_sph'))
    form.flexaid_binding_site_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_binding_site, form.output_box, filter_for='_sph'))
    form.flexaid_button_start.clicked.connect(lambda: form.flexaid_tab.setTabEnabled(2, True))
    form.flexaid_retrieve_nrgdock_ligands.clicked.connect(lambda: flexaid.retrieve_nrgdock_ligands(nrgdock_output_path))
    form.flexaid_button_start.clicked.connect(lambda: flexaid.run_flexaid(form, temp_path, binary_folder_path, operating_system, binary_suffix, install_dir))
    form.flexaid_button_pause.clicked.connect(lambda: flexaid.pause_resume_simulation(form))
    form.flexaid_button_stop.clicked.connect(lambda: flexaid.stop_simulation(form))
    form.flexaid_button_abort.clicked.connect(lambda: flexaid.abort_simulation(form))

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
    # form.surfaces_refresh_button.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_result, form.output_box, filter_for='RESULT'))
    # form.surfaces_run_button.clicked.connect(lambda: surfaces.run_run_surfaces(form.surface_select_result.currentText(), surfaces_output_path, form.simulate_folder_path.text(), binary_folder_path, binary_suffix, install_dir))
    # form.surfaces_run_button.clicked.connect(lambda: general_functions.surfaces_enable_buttons(form))
    # form.surfaces_retreive_flexaid_result.clicked.connect(lambda: surfaces.retrieve_flexaid_result(form.simulate_folder_path.text()))
    # form.surfaces_retreive_flexaid_result.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_result, form.output_box, filter_for='RESULT'))
    # form.surfaces_result_browse_button.clicked.connect(
    #     lambda: general_functions.folder_browser(form.surfaces_load_result_text, os.path.join(install_dir, 'result_demo'), "PDB file (*.pdb)"))
    # form.surfaces_load_result_button.clicked.connect(
    #     lambda: surfaces.load_surfaces_result(form, surfaces_output_path))

    # form.class_test.clicked.connect(lambda: getcleft.test_submit_command())

    # surfaces functions

    form.surfaces_refresh_button.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_result, form.output_box))
    form.surfaces_refresh_button.clicked.connect(lambda: run_NRGTEN.refresh_dropdown_NRG(form.surface_select_lig, form.output_box))
    form.surfaces_refresh_button_2.clicked.connect(lambda: run_NRGTEN.refresh_dropdown_NRG(form.surface_select_result_2, form.output_box))
    form.surfaces_refresh_button_2.clicked.connect(lambda: run_NRGTEN.refresh_dropdown_NRG(form.surface_select_lig_2, form.output_box))

    form.surfaces_run_button.clicked.connect(lambda: run_Surfaces.load_surfaces(form, temp_path, install_dir, binary_folder_path, binary_suffix))

    # nrgten functions
    form.NRGten_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.NRGten_select_target, form.output_box))
    form.NRGten_target_refresh.clicked.connect(lambda: run_NRGTEN.refresh_dropdown_NRG(form.NRGten_select_ligand, form.output_box))
    form.NRGten_target_refresh_2.clicked.connect(lambda: run_NRGTEN.refresh_dropdown_NRG(form.NRGten_select_target_2, form.output_box))

    form.NRGten_dynasig_pushButton.clicked.connect(lambda: run_NRGTEN.dynamical_signature(form.NRGten_select_target.currentText(),
                                                                                          form.NRGten_select_ligand.currentText(),
                                                                                          form.NRGten_select_target_2.currentText(),
                                                                                          form.NRGten_dynasig_lineEdit.text(), install_dir, temp_path))
    form.NRGten_conf_ensem_pushButton.clicked.connect(lambda: run_NRGTEN.conformational_ensemble(form.NRGten_select_target.currentText(),
                                                                                                 form.NRGten_modes_lineEdit.text(),
                                                                                                 form.NRGten_step_lineEdit.text(),
                                                                                                 form.NRGten_max_conf_lineEdit.text(),
                                                                                                 form.NRGten_max_dis_lineEdit.text(),
                                                                                                 form.NRGten_optmizestates.isChecked(), install_dir, temp_path))

    # modeller functions
    form.Modeller_target_refresh_1.clicked.connect(lambda: general_functions.refresh_dropdown(form.Modeller_select_target_1, form.output_box))
    form.Modeller_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.Modeller_select_target, form.output_box))
    form.Modeller_pushButton.clicked.connect(lambda: run_modeller.model_mutations(form, temp_path))
    form.Modeller_checkBox_all.clicked.connect(lambda: run_modeller.check_all(form))
    return dialog
