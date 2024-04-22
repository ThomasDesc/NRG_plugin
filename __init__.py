'''
PyMOL Demo Plugin

The plugin resembles the old "Rendering Plugin" from Michael Lerner, which
was written with Tkinter instead of PyQt.

(c) Schrodinger, Inc.

License: BSD-2-Clause
'''

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
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi

    import os
    import sys
    import shutil

    install_dir = os.path.dirname(__file__)
    sys.path.append(install_dir)
    import getcleft
    import general_functions
    import flexaid
    import platform
    import spheres
    import run_Surfaces
    import thread_test
    dialog = QtWidgets.QDialog()

    OS = platform.system().upper()
    folder = ''
    suffix = ''
    OSid = 'UNKNOWN'
    # SET the configuration location
    if OS == 'LINUX' or OS == 'BSD':
        OSid = 'LINUX'
        folder = 'linux'
    elif OS == 'DARWIN':
        OSid = 'MAC'
        folder = 'mac'
    elif OS == 'WINDOWS' or OS == 'MICROSOFT' or OS == 'WIN32':
        OSid = 'WIN'
        folder = 'win'
        suffix = '.exe'
    else:
        OSid = 'UNKNOWN'

    main_folder_path = os.path.dirname(__file__)
    uifile = os.path.join(os.path.dirname(__file__), 'demowidget.ui')
    getcleft_path = os.path.join(os.path.dirname(__file__), 'bin', folder, f'GetCleft{suffix}')
    flexaid_path = os.path.join(os.path.dirname(__file__), 'bin', folder, f'FlexAID{suffix}')
    vcon_path = os.path.join(os.path.dirname(__file__), 'bin', folder, f'vcon')
    process_ligand_path = os.path.join(os.path.dirname(__file__), 'bin', folder, f'Process_ligand{suffix}')
    plugin_tmp_output_path = os.path.join(os.path.expanduser('~'), 'Documents', 'NRGSuite_Qt')
    temp_path = os.path.join(plugin_tmp_output_path, 'temp')
    getcleft_output_path = os.path.join(temp_path, 'GetCleft')
    flexaid_output_path = os.path.join(temp_path, 'FlexAID')
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
    form = loadUi(uifile, dialog)
    form.stackedWidget.setCurrentIndex(0)
    form.flexaid_tab.setTabEnabled(2, False)


    # Refresh object dropdown menu
    general_functions.refresh_dropdown(form.cleft_select_object, form.output_box, no_warning=True)
    form.button_getcleft.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(0))
    form.button_partition_cleft.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(1))
    form.button_flexaid.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(2))
    form.button_nrgdock.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(3))
    form.button_Surfaces.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(4))

    form.cleft_partition_button_add.clicked.connect(lambda: spheres.display_sphere(form.cleft_partition_select_object.currentText(), form.cleft_partition_radius_slider, form.partition_sphere_select))

    form.button_hide.clicked.connect(lambda: general_functions.pymol_hide_structures(form))
    form.cleft_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_select_object, form.output_box))
    form.cleft_partition_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_partition_select_object, form.output_box, filter_for='_sph_'))
    form.cleft_partition_button_move.clicked.connect(lambda: spheres.move_sphere(form.partition_sphere_select.currentText()))

    form.button_start.clicked.connect(lambda: getcleft.run_getcleft(form, getcleft_path, getcleft_output_path,
                                                                    cleft_save_path, color_list))
    form.flexaid_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_target, form.output_box))
    form.flexaid_ligand_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_ligand, form.output_box))
    form.flexaid_binding_site_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_binding_site, form.output_box, filter_for='_sph_'))
    form.flexaid_button_start.clicked.connect(lambda: form.flexaid_tab.setTabEnabled(2, True))
    form.flexaid_button_start.clicked.connect(lambda: flexaid.run_flexaid(flexaid_output_path, form, cleft_save_path, process_ligand_path, flexaid_path, simulation_folder_path, hex_colour_list, plugin_tmp_output_path))

    form.flexaid_button_pause.clicked.connect(lambda: flexaid.pause_resume_simulation(form))
    form.flexaid_button_stop.clicked.connect(lambda: flexaid.stop_simulation(form))
    form.flexaid_button_abort.clicked.connect(lambda: flexaid.abort_simulation(form))

    form.cleft_partition_radius_slider.valueChanged.connect(lambda: spheres.resize_sphere(form.partition_sphere_select.currentText(), form.cleft_partition_radius_slider.value()))
    form.cleft_partition_crop_button.clicked.connect(lambda: spheres.crop_cleft(form.partition_sphere_select.currentText(), form.cleft_partition_radius_slider.value()/100, cleft_save_path, form.cleft_partition_select_object.currentText()))

    form.surfaces_refresh_button.clicked.connect(lambda: general_functions.refresh_dropdown(form.surface_select_result, form.output_box, filter_for='RESULT'))
    form.surfaces_run_button.clicked.connect(lambda: run_Surfaces.run_run_surfaces(form.surface_select_result.currentText(), surfaces_output_path, form.simulate_folder_path.text(), main_folder_path, vcon_path))
    # form.class_test.clicked.connect(lambda: getcleft.test_submit_command())
    return dialog
