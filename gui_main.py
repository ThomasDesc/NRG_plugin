import os
import sys

install_dir = os.path.dirname(__file__)
sys.path.append(install_dir)
import shutil
import subprocess
from src.flexaid.flexaid import FlexAIDManager, stop_simulation, abort_simulation, pause_resume_simulation, flexaid_show_ligand_from_table
from src.getcleft import getcleft
from src.nrgrank import nrgrank_on_target
from src.getcleft import spheres
import general_functions
from src.surfaces import run_Surfaces
from src.isomif import run_isomif
from src.nrgrank import nrgrank_smiles_management
import platform
from PyQt5.QtWidgets import QWidget, QTableWidget
from PyQt5.QtGui import QStandardItemModel
from PyQt5.uic import loadUi
from src.nrgten import run_NRGTEN
try:
    import modeller
except ImportError:
    print('Modeller not installed.')
else:
    from src.modeller import run_modeller
# TODO: when showing surfaces result hide everything else
# TODO: clickable results in nrgrank table

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
    def __init__(self, form, binary_folder_path, binary_suffix, operating_system, ligand_set_folder_path, color_list):
        self.form = form
        self.binary_folder_path = binary_folder_path
        self.binary_suffix = binary_suffix
        self.operating_system = operating_system
        self.ligand_set_folder_path = ligand_set_folder_path
        self.color_list = color_list
        # self.form.nrgrank_result_table.setSelectionMode(QTableWidget.MultiSelection)
        self.model = QStandardItemModel()
        self.form.nrgrank_result_table.setModel(self.model)
        self.form.flexaid_result_table.setModel(self.model)
        self.setupConnections()

    def setupConnections(self):
        self.form.button_getcleft.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(0))
        self.form.button_nrgrank.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(1))
        self.form.button_flexaid.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(2))
        self.form.button_surfaces.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(3))
        self.form.button_nrgten.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(4))
        self.form.button_modeller.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(5))
        self.form.button_ISOMIF.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(6))
        self.form.button_settings.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(7))
        self.form.button_home.clicked.connect(lambda: self.form.stackedWidget.setCurrentIndex(8))

        # Save/load session
        self.form.button_save.clicked.connect(lambda: general_functions.show_save_dialog(self.form, self.form.temp_line_edit.text()))
        self.form.button_load.clicked.connect(lambda: general_functions.show_save_dialog(self.form, self.form.temp_line_edit.text(), save=0))

        # GetCleft
        self.form.button_hide.clicked.connect(lambda: general_functions.pymol_hide_structures(self.form))
        self.form.cleft_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.cleft_select_object, self.form.output_box))
        self.form.button_start.clicked.connect(self.run_getcleft)

        # Partition Cleft
        self.form.cleft_partition_button_add.clicked.connect(
            lambda: spheres.display_sphere(self.form.cleft_partition_select_object.currentText(), self.form,
                                           self.form.cleft_partition_radius_slider, self.form.temp_line_edit.text()))
        self.form.cleft_partition_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.cleft_partition_select_object, self.form.output_box,filter_for='bd_site'))
        self.form.cleft_partition_button_move.clicked.connect(spheres.move_sphere)
        self.form.cleft_partition_radius_slider.valueChanged.connect( lambda: spheres.resize_sphere('SPHERE', self.form.cleft_partition_radius_slider.value()))
        self.form.cleft_partition_crop_button.clicked.connect(lambda: spheres.crop_cleft('SPHERE', self.form.cleft_partition_radius_slider.value() / 100, self.form.temp_line_edit.text(), self.form.cleft_partition_select_object.currentText(), self.form.output_box, self.form.cleft_partition_radius_slider))
        self.form.cleft_partition_button_delete.clicked.connect(lambda: spheres.delete_sphere('SPHERE', self.form.cleft_partition_radius_slider))

        # NRGRank:
        self.form.nrgrank_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown_target(self.form.nrgrank_select_target, self.form.output_box))
        self.form.nrgrank_select_target.currentIndexChanged.connect(lambda: general_functions.refresh_dropdown_bd_site(self.form.nrgrank_select_binding_site, self.form.nrgrank_select_target.currentText(), self.form.output_box, show_all_objects=self.form.show_all_obj_bd_checkbox.isChecked()))
        self.form.nrgrank_ligand_set_refresh.clicked.connect(lambda: general_functions.refresh_folder(self.ligand_set_folder_path, self.form.nrgrank_select_ligand))
        self.form.nrgrank_button_start.clicked.connect(self.run_nrgrank)
        self.form.nrgrank_button_cancel.clicked.connect(self.abort_nrgrank)
        self.form.nrgrank_result_browse_button.clicked.connect(lambda: general_functions.folder_browser(self.form.nrgrank_result_path, os.path.join(self.form.temp_line_edit.text(), 'NRGRank'), "CSV file (*.csv)"))
        self.form.nrgrank_result_table.selectionModel().selectionChanged.connect(lambda: nrgrank_on_target.show_ligand_from_table(self.form.nrgrank_result_table, self.form.nrgrank_select_binding_site.currentText(), self.form.nrgrank_select_ligand.currentText()))

        # NRGRank ligand manager:
        self.form.nrgrank_add_ligandset_button.clicked.connect(lambda: general_functions.folder_browser(self.form.nrgrank_add_ligand_file_path, self.ligand_set_folder_path, "Smiles Files (*.smi)"))
        self.form.nrgrank_delete_ligand_set_refresh.clicked.connect(lambda: general_functions.refresh_folder(self.ligand_set_folder_path, self.form.nrgrank_delete_ligand_set_dropdown, ignore_defaults=True))
        self.form.nrgrank_ligand_set_delete.clicked.connect(lambda: nrgrank_smiles_management.delete_ligand_set(self.form.nrgrank_delete_ligand_set_dropdown.currentText(), self.ligand_set_folder_path, self.form.output_box))
        self.form.nrgrank_button_ligandset_add.clicked.connect(self.run_generate_conformers)

        # FlexAID:
        self.form.flexaid_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown_target(self.form.flexaid_select_target,self.form.output_box))
        self.form.flexaid_select_target.currentIndexChanged.connect(lambda: general_functions.refresh_dropdown_bd_site(self.form.flexaid_select_binding_site, self.form.flexaid_select_target.currentText(), self.form.output_box, show_all_objects=self.form.show_all_obj_bd_checkbox.isChecked()))
        self.form.flexaid_ligand_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.flexaid_select_ligand, self.form.output_box))
        self.form.flexaid_button_start.clicked.connect(self.run_flexaid)
        self.form.flexaid_button_pause.clicked.connect(lambda: pause_resume_simulation(self.form, self.flexaid_manager.run_specific_simulate_folder_path))
        self.form.flexaid_button_stop.clicked.connect(lambda: stop_simulation(self.form, self.flexaid_manager.run_specific_simulate_folder_path))
        self.form.flexaid_button_abort.clicked.connect(lambda: abort_simulation(self.form, self.flexaid_manager.run_specific_simulate_folder_path, self.flexaid_manager))
        self.form.flexaid_result_table.selectionModel().selectionChanged.connect(lambda: flexaid_show_ligand_from_table(self.form))

        # Surfaces
        self.form.surfaces_refresh_object_1.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.surface_select_object_1, self.form.output_box))
        self.form.surfaces_refresh_object_1.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.surface_select_ligand_object_1, self.form.output_box, lig=1, add_none=1))
        self.form.surfaces_refresh_object_2.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.surface_select_object_2, self.form.output_box, add_none=1))
        self.form.surfaces_refresh_object_2.clicked.connect(lambda: general_functions.refresh_dropdown(self.form.surface_select_ligand_object_2, self.form.output_box, lig=1, add_none=1))
        self.form.surfaces_button_run.clicked.connect(lambda: run_Surfaces.load_surfaces(self.form, self.form.temp_line_edit.text(), install_dir, self.binary_folder_path, self.binary_suffix))
        self.form.surface_select_individual_result.currentIndexChanged.connect(lambda: run_Surfaces.load_csv_data(self.form, os.path.join(
            os.path.join(self.form.temp_line_edit.text(), 'Surfaces'), self.form.surface_select_individual_result.currentText() + '.txt')))
        self.form.surface_select_cf_comparison.currentIndexChanged.connect(lambda: run_Surfaces.load_csv_data(self.form, os.path.join(
            os.path.join(self.form.temp_line_edit.text(), 'Surfaces'), self.form.surface_select_cf_comparison.currentText() + '.csv')))
        self.form.surfaces_refresh_result.clicked.connect(lambda: run_Surfaces.refresh_res(self.form, os.path.join(self.form.temp_line_edit.text(), 'Surfaces')))
        self.form.surfaces_refresh_result.clicked.connect(lambda: run_Surfaces.load_csv_data(self.form, os.path.join(self.form.temp_line_edit.text(), 'Surfaces', self.form.surface_select_cf_comparison.currentText() + '.csv')))
        self.form.surfaces_button_interface.clicked.connect(lambda: run_Surfaces.read_and_select_residues(os.path.join(self.form.temp_line_edit.text(), 'Surfaces', self.form.surface_select_individual_result.currentText() + '.txt'),
            self.form.surface_select_individual_result.currentText()[5:-11], num_rows=self.form.TOPN_lineEdit_2.text()))

        # NRGTEN
        self.form.NRGten_target_refresh_object_1.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.NRGten_select_target_object_1, self.form.output_box))
        self.form.NRGten_target_refresh_object_1.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.NRGten_select_ligand_object_1, self.form.output_box, lig=1, add_none=1))
        self.form.NRGten_target_refresh_object_2.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.NRGten_select_target_object_2, self.form.output_box, add_none=1))
        self.form.NRGten_dynasig_run.clicked.connect(
            lambda: run_NRGTEN.dynamical_signature(self.form.NRGten_select_target_object_1.currentText(),
                                                   self.form.NRGten_select_ligand_object_1.currentText(),
                                                   self.form.NRGten_select_target_object_2.currentText(),
                                                   self.form.NRGten_dynasig_beta.text(), install_dir,
                                                   self.form.temp_line_edit.text()))
        self.form.NRGten_conf_ensem_run.clicked.connect(
            lambda: run_NRGTEN.conformational_ensemble(self.form.NRGten_select_target_object_1.currentText(),
                                                       self.form.NRGten_modes_lineEdit.text(),
                                                       self.form.NRGten_step_lineEdit.text(),
                                                       self.form.NRGten_max_conf_lineEdit.text(),
                                                       self.form.NRGten_max_dis_lineEdit.text(),
                                                       self.form.NRGten_optmizestates.isChecked(), install_dir,
                                                       self.form.temp_line_edit.text(), self.form))

        # Modeller
        self.form.Modeller_refresh_object.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.Modeller_select_object, self.form.output_box))
        self.form.Modeller_refresh_residue.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.Modeller_select_residue, self.form.output_box, lig=1))
        self.form.modeller_button_mutate.clicked.connect(lambda: run_modeller.model_mutations(self.form, self.form.temp_line_edit.text()))
        self.form.modeller_checkbox_all.clicked.connect(lambda: run_modeller.check_all(self.form))

        # isomif functions
        self.form.ISOMIF_target_refresh_object_1.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_target_object_1, self.form.output_box))
        self.form.ISOMIF_target_refresh_object_2.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_target_object_2, self.form.output_box, add_none=1))

        self.form.ISOMIF_cleft_refresh_object_1.clicked.connect(
            lambda: general_functions.refresh_dropdown_bd_site(self.form.ISOMIF_select_cleft_object_1, self.form.ISOMIF_select_target_object_1.currentText(), self.form.output_box, show_all_objects=self.form.show_all_obj_bd_checkbox.isChecked()))
        self.form.ISOMIF_cleft_refresh_object_1.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_ligand_object_1, self.form.output_box, lig=1, add_none=1))

        self.form.ISOMIF_cleft_refresh_object_2.clicked.connect(
            lambda: general_functions.refresh_dropdown_bd_site(self.form.ISOMIF_select_cleft_object_2, self.form.ISOMIF_select_target_object_2.currentText(), self.form.output_box, add_none=True, show_all_objects=self.form.show_all_obj_bd_checkbox.isChecked()))
        self.form.ISOMIF_cleft_refresh_object_2.clicked.connect(
            lambda: general_functions.refresh_dropdown(self.form.ISOMIF_select_ligand_object_2, self.form.output_box, lig=1, add_none=1))

        self.form.ISOMIF_run.clicked.connect(
            lambda: run_isomif.mif_plot(self.form, self.binary_folder_path, self.binary_suffix, install_dir))
    def run_getcleft(self):
        try:
            self.getcleftrunner = getcleft.GetCleftRunner(self.form, self.binary_folder_path, self.binary_suffix, install_dir, self.color_list)
        except ValueError as e:
            general_functions.output_message(self.form.output_box, e, 'warning')
        else:
            self.getcleftrunner.run_task()

    def run_nrgrank(self):
        self.nrgrankrunner = nrgrank_on_target.NRGRankManager(self.form, install_dir, self.ligand_set_folder_path, self.model)
        self.nrgrankrunner.run_nrgrank()

    def abort_nrgrank(self):
        self.nrgrankrunner.handle_thread_finished()
        self.nrgrankrunner = None

    def run_generate_conformers(self):
        self.conformer_generator = nrgrank_smiles_management.ConfGeneratorManager(self.form, install_dir, self.ligand_set_folder_path)
        self.conformer_generator.generate_conformer()

    def run_flexaid(self):
        self.flexaid_manager = FlexAIDManager(self.form, self.binary_folder_path, self.binary_suffix, install_dir, self.color_list, self.model)
        self.flexaid_manager.start_run()


class NRGSuitePlugin(QWidget):
    def __init__(self):
        super().__init__()
        self.form = loadUi(os.path.join(install_dir, 'plugin.ui'), self)
        self.binary_suffix = None
        self.operating_system = None
        self.get_os()
        self.binary_folder_path = os.path.join(install_dir, 'bin', self.operating_system)
        test_binary(self.binary_folder_path, self.operating_system)
        self.get_folders()
        self.manage_dirs()
        self.check_modeller()
        self.form.stackedWidget.setCurrentIndex(0)
        self.form.flexaid_tab.setTabEnabled(2, False)
        self.form.NRGRank_tabs.setTabEnabled(2, False)
        general_functions.refresh_dropdown(self.form.cleft_select_object, self.form.output_box, no_warning=True)
        general_functions.refresh_folder(self.ligand_set_folder_path, self.form.nrgrank_select_ligand)
        self.form.nrgrank_cpu_usage_target.setCurrentText("75%")
        self.color_list = general_functions.load_color_list(os.path.join(install_dir, 'deps', 'getcleft', 'color_list.txt'))
        self.form.nrgrank_progress_label.setText('')
        self.form.nrgrank_loading_gif.setText('')
        self.form.nrgrank_progress.hide()
        self.controller = Controller(self.form, self.binary_folder_path, self.binary_suffix, self.operating_system, self.ligand_set_folder_path, self.color_list)

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
        self.ligand_set_folder_path = os.path.join(install_dir, 'nrgrank_ligand_sets')
        self.plugin_tmp_output_path = os.path.join(os.path.expanduser('~'), 'Documents', 'NRGSuite_Qt')
        self.temp_path = os.path.join(self.plugin_tmp_output_path, 'temp')
        self.form.temp_line_edit.setText(self.temp_path)
        self.nrgrank_output_path = os.path.join(self.form.temp_line_edit.text(), 'NRGRank')
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
        os.mkdir(self.nrgrank_output_path)
        os.mkdir(self.modeller_save_path)
        os.mkdir(self.nrgten_save_path)
        os.mkdir(self.isomif_save_path)

    def check_modeller(self):
        if 'modeller' not in sys.modules:
            general_functions.output_message(self.form.output_box, 'Modeller install not detected. '
                                                                   'The modeller tab will be unavailable. '
                                                                   'It will not possible to optimise states in the NRGTEN tab. '
                                                                   '\nPlease install via conda.', 'warning')
            general_functions.output_message(self.form.output_box, '=====================', 'warning')
            self.form.NRGten_optmizestates.setEnabled(False)
            self.form.button_modeller.setEnabled(False)
            self.form.button_modeller.setStyleSheet("background-color: black; color: white;")