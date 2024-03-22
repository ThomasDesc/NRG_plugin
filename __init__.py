'''
PyMOL Demo Plugin

The plugin resembles the old "Rendering Plugin" from Michael Lerner, which
was written with Tkinter instead of PyQt.

(c) Schrodinger, Inc.

License: BSD-2-Clause
'''

from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

import os
import subprocess
import shutil

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('GetCleft', run_plugin_gui)


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
    # entry point to PyMOL's API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt

    dialog = QtWidgets.QDialog()

    uifile = os.path.join(os.path.dirname(__file__), 'demowidget.ui')
    getcleft_path = os.path.join(os.path.dirname(__file__), 'bin', 'GetCleft')
    plugin_tmp_output_path = os.path.join(os.path.expanduser('~'), 'Documents', 'GetCleft_tmp')
    getcleft_output_path = os.path.join(plugin_tmp_output_path, 'get_cleft_output')
    color_list = ['red', 'br8', 'tv_red', 'oxygen', 'iron', 'tv_orange', 'sulfur', 'gold', 'yelloworange', 'neodymium',
                  'limon', 'chartreuse', 'tv_green', 'limegreen', 'teal', 'rhodium', 'slate', 'tv_blue', 'blue',
                  'density']

    if os.path.isdir(plugin_tmp_output_path):
        shutil.rmtree(plugin_tmp_output_path)
    os.mkdir(plugin_tmp_output_path)
    os.mkdir(getcleft_output_path)
    form = loadUi(uifile, dialog)

    def select_object():
        list_pymol_objects = cmd.get_names('all')
        form.select_object.clear()
        form.select_object.addItems(list_pymol_objects)
        form.select_object.setCurrentText(list_pymol_objects[0])

    def get_arg_str(getcleft_path, object_path):
        min_radius = form.input_min_radii.text()
        max_radius = form.input_max_radii.text()
        resnumc = form.input_residue_in_contact.text()
        max_cleft_show = form.input_max_cleft_show.text()
        print('getcleft output path: ', getcleft_path)
        getcleft_output_name = os.path.join(getcleft_output_path, 'receptor')
        arg_string = f'{getcleft_path} -p "{object_path}" -l {min_radius} -u {max_radius} -t {max_cleft_show} -o "{getcleft_output_name}" -s'
        if resnumc != "":
            arg_string += f' -a {resnumc}'
        return arg_string

    def create_number_list(pTotColor, pTotalColorList):
        number_list = []
        modulo = (pTotalColorList - 1) % (pTotColor - 1)
        partition = (pTotalColorList - modulo - 1) / (pTotColor - 1)
        step_start = 0
        step_end = pTotalColorList - 1
        for i in range(0, pTotColor):

            if ((i % 2) == 0):
                number_list.append(step_start)
                step_start = step_start + partition
            else:
                number_list.append(step_end)
                step_end = step_end - partition
        number_list.sort()
        return [int(i) for i in number_list]

    def load_show_cleft(getcleft_output_path, color_list):
        auto_zoom = cmd.get("auto_zoom")
        cmd.set("auto_zoom", 0)
        all_files = os.listdir(getcleft_output_path)
        sph_file_list = []
        for filename in all_files:
            if filename.find('sph') != -1:
                sph_file_list.append({'path': os.path.join(getcleft_output_path, filename),
                                      'name': filename.split('.')[0]})
        # for Cleft in self.TempBindingSite.listClefts:
        sph_file_list = sorted(sph_file_list, key=lambda d: d['name'])
        number_color_list = create_number_list(len(sph_file_list), len(color_list))
        for cleft_counter, Cleft in enumerate(sph_file_list):
            try:
                cmd.load(Cleft['path'], Cleft['name'], state=1)
                cmd.hide('everything', Cleft['name'])
                if cleft_counter >= len(color_list):
                    cmd.color('grey50', Cleft['name'])
                else:
                    cmd.color(color_list[number_color_list[cleft_counter]], Cleft['name'])
                cmd.show('surface', Cleft['name'])
            except:
                print(f"ERROR: Failed to load cleft object  {Cleft['name']}")
                continue
        cmd.refresh()
        cmd.set("auto_zoom", auto_zoom)

    def run_getcleft(getcleft_path):
        pymol_object = form.select_object.currentText()
        if pymol_object != '':

            object_save_path = os.path.join(plugin_tmp_output_path, 'tmp.pdb')
            cmd.save(object_save_path, pymol_object)
        else:
            print('No object selected')
            return
        getcleft_command = f"{get_arg_str(getcleft_path, object_save_path)}"
        print(getcleft_command)
        subprocess.run(getcleft_command, shell=True)
        load_show_cleft(getcleft_output_path, color_list)

    # Refresh object dropdown menu
    select_object()

    form.button_refresh.clicked.connect(select_object)
    form.button_start.clicked.connect(lambda: run_getcleft(getcleft_path))

    return dialog
