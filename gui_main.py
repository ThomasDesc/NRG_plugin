import os
import sys
install_dir = os.path.dirname(__file__)
sys.path.append(install_dir)
import shutil
import subprocess
from src.flexaid import flexaid
from src.getcleft import getcleft
from src.nrgdock import nrgdock
from src.getcleft import spheres
import general_functions
from src.surfaces import run_Surfaces
from src.isomif import run_isomif
from thread_test import main as thread_test
import platform
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi
from functools import partial

try:
    import modeller
except ImportError:
    print('Modeller not installed.')
else:
    from src.nrgten import run_NRGTEN
    from src.modeller import run_modeller


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


class Controller:
    def __init__(self, form, binary_folder_path, binary_suffix, operating_system, ligand_set_folder_path):
        self.form = form
        self.binary_folder_path = binary_folder_path
        self.binary_suffix = binary_suffix
        self.operating_system = operating_system
        self.ligand_set_folder_path = ligand_set_folder_path
        self.getcleftrunner = getcleft.GetCleftRunner()
        self.nrgdockrunner = nrgdock.NRGDockRunner()
        self.setupConnections()

    def setupConnections(self):
        self.form.button_getcleft.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(0))
        self.form.button_flexaid.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(1))
        self.form.button_nrgdock.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(2))
        self.form.button_nrgten.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(3))
        self.form.button_surfaces.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(4))
        self.form.button_modeller.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(5))
        self.form.button_ISOMIF.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(6))

        # save/load
        self.form.button_save.clicked.connect(lambda: general_functions.show_save_dialog(self.form, self.form.temp_line_edit.text()))
        self.form.button_load.clicked.connect(
            lambda: general_functions.show_save_dialog(self.form, self.form.temp_line_edit.text(), save=0))

        # GetCleft
        self.form.button_hide.clicked.connect(lambda: general_functions.pymol_hide_structures(self.form))
        self.form.cleft_button_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.cleft_select_object, self.form.output_box))
        self.form.button_start.clicked.connect(lambda: self.getcleftrunner.run_task(self.form, self.binary_folder_path, self.binary_suffix, install_dir))

        # Partition Cleft
        self.form.cleft_partition_button_add.clicked.connect(
            lambda: spheres.display_sphere(self.form.cleft_partition_select_object.currentText(),
                                           self.form.cleft_partition_radius_slider, self.form.partition_sphere_select,
                                           self.form.temp_line_edit.text()))
        self.form.cleft_partition_button_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.cleft_partition_select_object, self.form.output_box,
                                                       filter_for='_sph'))
        self.form.cleft_partition_button_move.clicked.connect(
            lambda: spheres.move_sphere(self.form.partition_sphere_select.currentText()))
        self.form.cleft_partition_radius_slider.valueChanged.connect(
            lambda: spheres.resize_sphere(self.form.partition_sphere_select.currentText(),
                                          self.form.cleft_partition_radius_slider.value()))
        self.form.cleft_partition_crop_button.clicked.connect(
            lambda: spheres.crop_cleft(self.form.partition_sphere_select.currentText(),
                                       self.form.cleft_partition_radius_slider.value() / 100, self.form.temp_line_edit.text(),
                                       self.form.cleft_partition_select_object.currentText()))

        # FlexAID:
        self.form.flexaid_target_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.flexaid_select_target, self.form.output_box, exclude='_sph'))
        self.form.flexaid_ligand_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.flexaid_select_ligand, self.form.output_box, exclude='_sph'))
        self.form.flexaid_binding_site_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.flexaid_select_binding_site, self.form.output_box,
                                                       filter_for='_sph'))
        self.form.flexaid_button_start.clicked.connect(lambda: self.form.flexaid_tab.setTabEnabled(2, True))
        self.form.flexaid_retrieve_nrgdock_ligands.clicked.connect(
            lambda: flexaid.retrieve_nrgdock_ligands(os.path.join(self.form.temp_line_edit.text(), 'NRGDock')))
        self.form.flexaid_button_start.clicked.connect(
            lambda: flexaid.run_flexaid(self.form, self.form.temp_line_edit.text(), self.binary_folder_path, self.operating_system,
                                        self.binary_suffix, install_dir))
        self.form.flexaid_button_pause.clicked.connect(lambda: flexaid.pause_resume_simulation(self.form))
        self.form.flexaid_button_stop.clicked.connect(lambda: flexaid.stop_simulation(self.form))
        self.form.flexaid_button_abort.clicked.connect(lambda: flexaid.abort_simulation(self.form))

        # NRGDock:
        self.form.nrgdock_target_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.nrgdock_select_target, self.form.output_box, exclude='_sph'))
        self.form.nrgdock_binding_site_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.nrgdock_select_binding_site, self.form.output_box,
                                                       filter_for='_sph'))
        self.form.nrgdock_delete_ligand_set_refresh.clicked.connect(
            lambda: general_functions.refresh_folder(self.ligand_set_folder_path, self.form.nrgdock_delete_ligand_set_dropdown))
        self.form.nrgdock_add_ligandset_button.clicked.connect(
            lambda: general_functions.folder_browser(self.form.nrgdock_add_ligand_file_path, self.ligand_set_folder_path,
                                                     "Smiles Files (*.smi)"))
        self.form.nrgdock_button_ligandset_add.clicked.connect(lambda: nrgdock.process_ligands())
        self.form.nrgdock_ligand_set_refresh.clicked.connect(
            lambda: general_functions.refresh_folder(self.ligand_set_folder_path, self.form.nrgdock_select_ligand))

        # self.form.nrgdock_button_start.clicked.connect(
        #     lambda: nrgdock.run_nrgdock(self.form, os.path.join(self.form.temp_line_edit.text(), 'NRGDock'),
        #                                 self.ligand_set_folder_path, install_dir))
        self.form.nrgdock_button_start.clicked.connect(partial(self.nrgdockrunner.run_task, self.form, self.ligand_set_folder_path, install_dir))
        self.form.nrgdock_result_browse_button.clicked.connect(
            lambda: general_functions.folder_browser(self.form.nrgdock_result_path,
                                                     os.path.join(self.form.temp_line_edit.text(), 'NRGDock'),
                                                     "CSV file (*.csv)"))
        self.form.nrgdock_load_csv_button.clicked.connect(
            lambda: nrgdock.get_nrgdock_result_model(self.form.nrgdock_result_path.text(), self.form))

        # Surfaces
        self.form.surfaces_refresh_button.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.surface_select_result, self.form.output_box))
        self.form.surfaces_refresh_button.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.surface_select_lig, self.form.output_box, lig=1, add_none=1))
        self.form.surfaces_refresh_button_2.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.surface_select_result_2, self.form.output_box, add_none=1))
        self.form.surfaces_refresh_button_2.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.surface_select_lig_2, self.form.output_box, lig=1, add_none=1))
        self.form.surfaces_run_button.clicked.connect(
            lambda: run_Surfaces.load_surfaces(self.form, self.form.temp_line_edit.text(), install_dir, self.binary_folder_path,
                                               self.binary_suffix))
        self.form.surface_select_result_3.currentIndexChanged.connect(lambda: run_Surfaces.load_csv_data(self.form, os.path.join(
            os.path.join(self.form.temp_line_edit.text(), 'Surfaces'), self.form.surface_select_result_3.currentText() + '.txt')))
        self.form.surface_select_result_4.currentIndexChanged.connect(lambda: run_Surfaces.load_csv_data(self.form, os.path.join(
            os.path.join(self.form.temp_line_edit.text(), 'Surfaces'), self.form.surface_select_result_4.currentText() + '.csv')))
        self.form.surfaces_refresh_button_3.clicked.connect(
            lambda: run_Surfaces.refresh_res(self.form, os.path.join(self.form.temp_line_edit.text(), 'Surfaces')))
        self.form.surfaces_refresh_button_3.clicked.connect(lambda: run_Surfaces.load_csv_data(self.form, os.path.join(
            self.form.temp_line_edit.text(), 'Surfaces', self.form.surface_select_result_4.currentText() + '.csv')))
        self.form.Surfaces_pushButton_2.clicked.connect(lambda: run_Surfaces.read_and_select_residues(
            os.path.join(self.form.temp_line_edit.text(), 'Surfaces', self.form.surface_select_result_3.currentText() + '.txt'),
            self.form.surface_select_result_3.currentText()[5:-11], num_rows=self.form.TOPN_lineEdit_2.text()))

        # NRGTEN
        self.form.NRGten_target_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.NRGten_select_target, self.form.output_box))
        self.form.NRGten_target_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.NRGten_select_ligand, self.form.output_box, lig=1, add_none=1))
        self.form.NRGten_target_refresh_2.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.NRGten_select_target_2, self.form.output_box, add_none=1))
        self.form.NRGten_dynasig_pushButton.clicked.connect(
            lambda: run_NRGTEN.dynamical_signature(self.form.NRGten_select_target.currentText(),
                                                   self.form.NRGten_select_ligand.currentText(),
                                                   self.form.NRGten_select_target_2.currentText(),
                                                   self.form.NRGten_dynasig_lineEdit.text(), install_dir,
                                                   self.form.temp_line_edit.text()))
        self.form.NRGten_conf_ensem_pushButton.clicked.connect(
            lambda: run_NRGTEN.conself.formational_ensemble(self.form.NRGten_select_target.currentText(),
                                                       self.form.NRGten_modes_lineEdit.text(),
                                                       self.form.NRGten_step_lineEdit.text(),
                                                       self.form.NRGten_max_conf_lineEdit.text(),
                                                       self.form.NRGten_max_dis_lineEdit.text(),
                                                       self.form.NRGten_optmizestates.isChecked(), install_dir,
                                                       self.form.temp_line_edit.text(), self.form))

        # Modeller
        self.form.Modeller_target_refresh_1.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.Modeller_select_target_1, self.form.output_box))
        self.form.Modeller_target_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.Modeller_select_target, self.form.output_box, lig=1))
        self.form.Modeller_pushButton.clicked.connect(lambda: run_modeller.model_mutations(self.form, self.form.temp_line_edit.text()))
        self.form.Modeller_checkBox_all.clicked.connect(lambda: run_modeller.check_all(self.form))

        # isomif functions
        self.form.ISOMIF_target_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_target, self.form.output_box))
        self.form.ISOMIF_target_refresh_1.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_target_1, self.form.output_box, add_none=1))
        self.form.ISOMIF_cleft_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_cleft, self.form.output_box, filter_for='sph'))
        self.form.ISOMIF_cleft_refresh.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_lig, self.form.output_box, lig=1, add_none=1))
        self.form.ISOMIF_cleft_refresh_1.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_cleft_1,
                                                       self.form.output_box,
                                                       filter_for='sph',
                                                       add_none=1))
        self.form.ISOMIF_cleft_refresh_1.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_lig_1, self.form.output_box, lig=1, add_none=1))

        self.form.ISOMIF_pushButton.clicked.connect(
            lambda: run_isomif.mif_plot(self.form, self.form.output_box, self.binary_folder_path, self.binary_suffix, self.operating_system,
                                        install_dir))


class NRGSuitePlugin(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.form = None
        self.binary_suffix = None
        self.operating_system = None
        # QtWidgets.QApplication.setStyle("Fusion")
        self.get_os()
        self.binary_folder_path = os.path.join(install_dir, 'bin', self.operating_system)
        print('binary path: ', self.binary_folder_path)
        test_binary(self.binary_folder_path, self.operating_system)
        self.load_ui()
        self.get_folders()
        self.manage_dirs()
        self.check_modeller()

        self.form.stackedWidget.setCurrentIndex(0)
        self.form.flexaid_tab.setTabEnabled(2, False)
        if self.operating_system == 'mac':
            self.form.flexaid_multithread_button.setChecked(True)
        self.form.getcleft_tab_widget.setTabEnabled(2, False)

        general_functions.refresh_dropdown(self.form.cleft_select_object, self.form.output_box, no_warning=True)
        general_functions.refresh_folder(self.ligand_set_folder_path, self.form.nrgdock_select_ligand)
        self.controller = Controller(self.form, self.binary_folder_path, self.binary_suffix, self.operating_system, self.ligand_set_folder_path)

    def load_ui(self):
        self.form = loadUi(os.path.join(install_dir, 'nrgdock_widget.ui'), self)

    def get_os(self):
        operating_system = platform.system().upper()
        self.binary_suffix = ''
        if operating_system == 'LINUX' or operating_system == 'BSD':
            self.operating_system = 'linux'
        elif operating_system == 'DARWIN':
            self.operating_system = 'mac'
        elif operating_system == 'WINDOWS' or operating_system == 'MICROSOFT' or operating_system == 'WIN32':
            self.operating_system = 'win'
        else:
            exit('Unknown operating system')

    def get_folders(self):
        self.ligand_set_folder_path = os.path.join(install_dir, 'nrgdock_ligand_sets')
        self.plugin_tmp_output_path = os.path.join(os.path.expanduser('~'), 'Documents', 'NRGSuite_Qt')
        self.temp_path = os.path.join(self.plugin_tmp_output_path, 'temp')
        self.form.temp_line_edit.setText(self.temp_path)
        self.nrgdock_output_path = os.path.join(self.form.temp_line_edit.text(), 'NRGDock')
        self.surfaces_output_path = os.path.join(self.form.temp_line_edit.text(), 'Surfaces')
        self.modeller_save_path = os.path.join(self.form.temp_line_edit.text(), 'modeller')
        self.nrgten_save_path = os.path.join(self.form.temp_line_edit.text(), 'NRGTEN')
        self.isomif_save_path = os.path.join(self.form.temp_line_edit.text(), 'ISOMIF')

    def manage_dirs(self):
        if os.path.isdir(self.plugin_tmp_output_path):
            shutil.rmtree(self.plugin_tmp_output_path)
        os.mkdir(self.plugin_tmp_output_path)
        os.mkdir(self.form.temp_line_edit.text())
        os.mkdir(self.surfaces_output_path)
        os.mkdir(self.nrgdock_output_path)
        os.mkdir(self.modeller_save_path)
        os.mkdir(self.nrgten_save_path)
        os.mkdir(self.isomif_save_path)

    def check_modeller(self):
        if 'modeller' not in sys.modules:
            general_functions.output_message(self.form.output_box, 'Modeller install not detected. '
                                                              'The modeller tab will be unavailable. '
                                                              'Please install via conda.', 'warning')
            self.form.button_nrgten.setEnabled(False)
            self.form.button_modeller.setEnabled(False)
            self.form.button_nrgten.setStyleSheet("background-color: black; color: white;")
            self.form.button_modeller.setStyleSheet("background-color: black; color: white;")