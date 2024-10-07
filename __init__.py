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
    except ImportError as e:
        if package == 'modeller':
            print('Modeller install not detected. Please install via conda. The modeller tab will be unavailable')
        else:
            if package == 'Bio':
                package = 'biopython'
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
    packages = ['nrgten', 'Bio', 'pandas', 'matplotlib', 'colour', 'scipy', 'numpy', 'numba','plotly']
    for package in packages:
        install_package(package, install_dir)
    from src.flexaid import flexaid
    from src.getcleft import getcleft
    from src.nrgdock import nrgdock
    from src.getcleft import spheres
    import general_functions
    from src.surfaces import run_Surfaces
    from src.isomif import run_isomif
    import platform
    dialog = QtWidgets.QDialog()

    OS = platform.system().upper()
    binary_suffix = ''
    if OS == 'LINUX' or OS == 'BSD':
        operating_system = 'linux'
    elif OS == 'DARWIN':
        operating_system = 'mac'
    elif OS == 'WINDOWS' or OS == 'MICROSOFT' or OS == 'WIN32':
        operating_system = 'win'
    else:
        exit('Unknown operating system')

    uifile = os.path.join(install_dir, 'nrgdock_widget.ui')
    form = loadUi(uifile, dialog)
    binary_folder_path = os.path.join(install_dir, 'bin', operating_system)
    print('binary path: ', binary_folder_path)
    test_binary(binary_folder_path, operating_system)

    ligand_set_folder_path = os.path.join(install_dir, 'nrgdock_ligand_sets')
    plugin_tmp_output_path = os.path.join(os.path.expanduser('~'), 'Documents', 'NRGSuite_Qt')
    temp_path = os.path.join(plugin_tmp_output_path, 'temp')
    form.temp_line_edit.setText(temp_path)
    nrgdock_output_path = os.path.join(form.temp_line_edit.text(), 'NRGDock')
    surfaces_output_path = os.path.join(form.temp_line_edit.text(), 'Surfaces')
    modeller_save_path = os.path.join(form.temp_line_edit.text(), 'modeller')
    nrgten_save_path =os.path.join(form.temp_line_edit.text(), 'NRGTEN')
    isomif_save_path = os.path.join(form.temp_line_edit.text(), 'ISOMIF')

    if os.path.isdir(plugin_tmp_output_path):
        shutil.rmtree(plugin_tmp_output_path)
    os.mkdir(plugin_tmp_output_path)
    os.mkdir(form.temp_line_edit.text())
    os.mkdir(surfaces_output_path)
    os.mkdir(nrgdock_output_path)
    os.mkdir(modeller_save_path)
    os.mkdir(nrgten_save_path)
    os.mkdir(isomif_save_path)
    try:
        import modeller
    except ModuleNotFoundError:
        general_functions.output_message(form.output_box, 'Modeller install not detected. '
                                                          'The modeller tab will be unavailable. '
                                                          'Please install via conda.', 'warning')
        form.button_nrgten.setEnabled(False)
        form.button_modeller.setEnabled(False)
        form.button_nrgten.setStyleSheet("background-color: black; color: white;")
        form.button_modeller.setStyleSheet("background-color: black; color: white;")
    else:
        from src.nrgten import run_NRGTEN
        from src.modeller import run_modeller

    # Disable isomif
    #form.button_ISOMIF.setEnabled(False)
    #form.button_ISOMIF.setStyleSheet("background-color: black; color: white;")

    form.stackedWidget.setCurrentIndex(0)
    form.flexaid_tab.setTabEnabled(2, False)
    if operating_system == 'mac':
        form.flexaid_multithread_button.setChecked(True)
    print(form.surface_select_result.currentText())
    form.getcleft_tab_widget.setTabEnabled(2, False)

    general_functions.refresh_dropdown(form.cleft_select_object, form.output_box, no_warning=True)
    general_functions.refresh_folder(ligand_set_folder_path, form.nrgdock_select_ligand)
    form.button_getcleft.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(0))
    form.button_flexaid.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(1))
    form.button_nrgdock.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(2))
    form.button_nrgten.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(3))
    form.button_surfaces.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(4))
    form.button_modeller.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(5))
    form.button_ISOMIF.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(6))

    # save/load
    form.button_save.clicked.connect(lambda: general_functions.show_save_dialog(form,form.temp_line_edit.text()))
    form.button_load.clicked.connect(lambda: general_functions.show_save_dialog(form,form.temp_line_edit.text(),save=0))

    # GetCleft
    form.button_hide.clicked.connect(lambda: general_functions.pymol_hide_structures(form))
    form.cleft_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_select_object, form.output_box))
    form.button_start.clicked.connect(lambda: getcleft.run_getcleft(form, binary_folder_path, binary_suffix, form.temp_line_edit.text(),
                                                                    install_dir))
    # Partition Cleft
    form.cleft_partition_button_add.clicked.connect(lambda: spheres.display_sphere(form.cleft_partition_select_object.currentText(), form.cleft_partition_radius_slider, form.partition_sphere_select, form.temp_line_edit.text()))
    form.cleft_partition_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_partition_select_object, form.output_box, filter_for='_sph'))
    form.cleft_partition_button_move.clicked.connect(lambda: spheres.move_sphere(form.partition_sphere_select.currentText()))
    form.cleft_partition_radius_slider.valueChanged.connect(lambda: spheres.resize_sphere(form.partition_sphere_select.currentText(), form.cleft_partition_radius_slider.value()))
    form.cleft_partition_crop_button.clicked.connect(lambda: spheres.crop_cleft(form.partition_sphere_select.currentText(), form.cleft_partition_radius_slider.value()/100, form.temp_line_edit.text(), form.cleft_partition_select_object.currentText()))

    # FlexAID:
    form.flexaid_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_target, form.output_box, exclude='_sph'))
    form.flexaid_ligand_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_ligand, form.output_box, exclude='_sph'))
    form.flexaid_binding_site_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_binding_site, form.output_box, filter_for='_sph'))
    form.flexaid_button_start.clicked.connect(lambda: form.flexaid_tab.setTabEnabled(2, True))
    form.flexaid_retrieve_nrgdock_ligands.clicked.connect(lambda: flexaid.retrieve_nrgdock_ligands(os.path.join(form.temp_line_edit.text(), 'NRGDock')))
    form.flexaid_button_start.clicked.connect(lambda: flexaid.run_flexaid(form, form.temp_line_edit.text(), binary_folder_path, operating_system, binary_suffix, install_dir))
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
        lambda: nrgdock.run_nrgdock(form, os.path.join(form.temp_line_edit.text(), 'NRGDock'), ligand_set_folder_path, install_dir))
    form.nrgdock_result_browse_button.clicked.connect(lambda: general_functions.folder_browser(form.nrgdock_result_path, os.path.join(form.temp_line_edit.text(), 'NRGDock'), "CSV file (*.csv)"))
    form.nrgdock_load_csv_button.clicked.connect(lambda: nrgdock.get_nrgdock_result_model(form.nrgdock_result_path.text(), form))

    # Surfaces
    form.surfaces_refresh_button.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_result, form.output_box))
    form.surfaces_refresh_button.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_lig, form.output_box, lig=1 ,    add_none=1  ))
    form.surfaces_refresh_button_2.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_result_2, form.output_box,     add_none=1  ))
    form.surfaces_refresh_button_2.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_lig_2, form.output_box,lig=1,     add_none=1  ))
    form.surfaces_run_button.clicked.connect(lambda: run_Surfaces.load_surfaces(form, form.temp_line_edit.text(), install_dir, binary_folder_path, binary_suffix))
    form.surface_select_result_3.currentIndexChanged.connect(lambda: run_Surfaces.load_csv_data(form,os.path.join(os.path.join(form.temp_line_edit.text(),'Surfaces'),form.surface_select_result_3.currentText()+'.txt')))
    form.surface_select_result_4.currentIndexChanged.connect(lambda: run_Surfaces.load_csv_data(form, os.path.join(
        os.path.join(form.temp_line_edit.text(), 'Surfaces'), form.surface_select_result_4.currentText() + '.csv')))
    form.surfaces_refresh_button_3.clicked.connect(lambda:run_Surfaces.refresh_res(form,os.path.join(form.temp_line_edit.text(),'Surfaces')))
    form.surfaces_refresh_button_3.clicked.connect(lambda: run_Surfaces.load_csv_data(form, os.path.join(form.temp_line_edit.text(), 'Surfaces', form.surface_select_result_4.currentText() + '.csv')))
    form.Surfaces_pushButton_2.clicked.connect(lambda: run_Surfaces.read_and_select_residues(os.path.join(form.temp_line_edit.text(),'Surfaces',form.surface_select_result_3.currentText()+'.txt'),form.surface_select_result_3.currentText()[5:-11],num_rows=form.TOPN_lineEdit_2.text()))

    # NRGTEN
    form.NRGten_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.NRGten_select_target, form.output_box))
    form.NRGten_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.NRGten_select_ligand, form.output_box,lig=1,     add_none=1  ))
    form.NRGten_target_refresh_2.clicked.connect(lambda: general_functions.refresh_dropdown(form.NRGten_select_target_2, form.output_box,     add_none=1  ))

    form.NRGten_dynasig_pushButton.clicked.connect(lambda: run_NRGTEN.dynamical_signature(form.NRGten_select_target.currentText(),
                                                                                          form.NRGten_select_ligand.currentText(),
                                                                                          form.NRGten_select_target_2.currentText(),
                                                                                          form.NRGten_dynasig_lineEdit.text(), install_dir, form.temp_line_edit.text()))
    form.NRGten_conf_ensem_pushButton.clicked.connect(lambda: run_NRGTEN.conformational_ensemble(form.NRGten_select_target.currentText(),
                                                                                                 form.NRGten_modes_lineEdit.text(),
                                                                                                 form.NRGten_step_lineEdit.text(),
                                                                                                 form.NRGten_max_conf_lineEdit.text(),
                                                                                                 form.NRGten_max_dis_lineEdit.text(),
                                                                                                 form.NRGten_optmizestates.isChecked(), install_dir, form.temp_line_edit.text(),form))

    # Modeller
    form.Modeller_target_refresh_1.clicked.connect(lambda: general_functions.refresh_dropdown(form.Modeller_select_target_1, form.output_box))
    form.Modeller_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.Modeller_select_target, form.output_box, lig=1 ))
    form.Modeller_pushButton.clicked.connect(lambda: run_modeller.model_mutations(form, form.temp_line_edit.text()))
    form.Modeller_checkBox_all.clicked.connect(lambda: run_modeller.check_all(form))

    #isomif functions
    form.ISOMIF_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.ISOMIF_select_target, form.output_box))
    form.ISOMIF_target_refresh_1.clicked.connect(lambda:general_functions.refresh_dropdown(form.ISOMIF_select_target_1, form.output_box,add_none=1))
    form.ISOMIF_cleft_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.ISOMIF_select_cleft, form.output_box,filter_for='sph'))
    form.ISOMIF_cleft_refresh_1.clicked.connect(lambda: general_functions.refresh_dropdown(form.ISOMIF_select_cleft_1,
                                                                                  form.output_box, filter_for='sph' ,add_none=1))
    form.ISOMIF_pushButton.clicked.connect(lambda: run_isomif.mif_plot(form, form.output_box,binary_folder_path, binary_suffix,operating_system,install_dir))



    return dialog
