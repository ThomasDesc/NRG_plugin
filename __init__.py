'''
PyMOL Demo Plugin

The plugin resembles the old "Rendering Plugin" from Michael Lerner, which
was written with Tkinter instead of PyQt.

(c) Schrodinger, Inc.

License: BSD-2-Clause
'''

from __future__ import absolute_import
from __future__ import print_function


def __init_plugin__(app=None):
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


def make_dialog():
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt
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


    uifile = os.path.join(os.path.dirname(__file__), 'demowidget2.ui')
    getcleft_path = os.path.join(os.path.dirname(__file__), 'bin', folder, f'GetCleft{suffix}')
    flexaid_path = os.path.join(os.path.dirname(__file__), 'bin', folder, f'FlexAID{suffix}')
    process_ligand_path = os.path.join(os.path.dirname(__file__), 'bin', folder, f'Process_ligand{suffix}')
    plugin_tmp_output_path = os.path.join(os.path.expanduser('~'), 'Documents', 'NRGSuite_Qt')
    temp_path = os.path.join(plugin_tmp_output_path, 'temp')
    getcleft_output_path = os.path.join(temp_path, 'GetCleft')
    flexaid_output_path = os.path.join(temp_path, 'FlexAID')
    cleft_save_path = os.path.join(getcleft_output_path, 'Clefts')
    simulation_folder_path = os.path.join(flexaid_output_path, 'Simulation')
    color_list = ['red', 'br8', 'tv_red', 'oxygen', 'iron', 'tv_orange', 'sulfur', 'gold', 'yelloworange', 'neodymium',
                  'limon', 'chartreuse', 'tv_green', 'limegreen', 'teal', 'rhodium', 'slate', 'tv_blue', 'blue',
                  'density']

    if os.path.isdir(plugin_tmp_output_path):
        shutil.rmtree(plugin_tmp_output_path)
    os.mkdir(plugin_tmp_output_path)
    os.mkdir(temp_path)
    os.mkdir(getcleft_output_path)
    os.mkdir(cleft_save_path)
    os.mkdir(flexaid_output_path)
    os.mkdir(simulation_folder_path)
    form = loadUi(uifile, dialog)

    # Refresh object dropdown menu
    general_functions.refresh_dropdown(form.cleft_select_object)
    form.button_getcleft.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(0))
    form.button_partition_cleft.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(1))
    form.button_flexaid.clicked.connect(lambda: form.stackedWidget.setCurrentIndex(2))

    form.cleft_partition_button_add.clicked.connect(lambda: spheres.display_sphere(form.cleft_partition_select_object.currentText(), form.cleft_partition_radius_slider))

    form.button_hide.clicked.connect(lambda: general_functions.pymol_hide_structures(form))
    form.cleft_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_select_object))
    form.cleft_partition_button_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.cleft_partition_select_object, filter_for='_sph_'))
    form.button_start.clicked.connect(lambda: getcleft.run_getcleft(form, getcleft_path, getcleft_output_path,
                                                                    cleft_save_path, color_list))
    form.flexaid_target_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_target))
    form.flexaid_ligand_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_ligand))
    form.flexaid_binding_site_refresh.clicked.connect(lambda: general_functions.refresh_dropdown(form.flexaid_select_binding_site, filter_for='_sph_'))
    form.flexaid_button_start.clicked.connect(lambda: flexaid.run_flexaid(flexaid_output_path, form, cleft_save_path, process_ligand_path, flexaid_path, simulation_folder_path))
    form.cleft_partition_radius_slider.valueChanged.connect(lambda: yo(form))


    def yo(form):
        print(form.cleft_partition_radius_slider.value())
    return dialog
