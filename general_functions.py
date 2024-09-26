from pymol import cmd
import numpy as np
import os
from PyQt5.QtWidgets import QFileDialog,QMessageBox
import shutil


def output_message(output_box, text, message_type):
    # Type can be error, warning or valid
    out_color = None
    red = '<span style="color:red;">{}</span>'
    yellow = '<span style="color:orange;">{}</span>'
    green = '<span style="color:green;">{}</span>'
    if message_type == 'error':
        out_color = red
    elif message_type == 'warning':
        out_color = yellow
    elif message_type == ' valid':
        out_color = green
    output_box.append(out_color.format(text))

def show_popup(self,dir_path,temp_path,save_file):

        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        if not save_file:
            msg.setWindowTitle('No sessions found in this project')
            msg.setText('No sessions found in this project,\nDo you want to proceed?')
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

        if save_file:
            msg.setWindowTitle('NRGSuite_QT results folder already exists')
            msg.setText('NRGSuite_QT results folder already exists,\nDo you want to proceed?')
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

        response = msg.exec_()
        if response == QMessageBox.Ok:
            if not save_file:
                self.temp_line_edit.setText(dir_path)
            if save_file:
                sufix = 2
                base_dir=os.path.join(dir_path, 'NRGSuite_QT_results')
                while os.path.exists(f"{base_dir}_{sufix}"):
                    sufix += 1
                new_dir = f"{base_dir}_{sufix}"
                os.mkdir(new_dir)
                files = os.listdir(temp_path)
                for file_name in files:
                    source_file = os.path.join(temp_path, file_name)
                    target_file = os.path.join(new_dir, file_name)
                    if '.DS_Store' not in file_name:
                        if os.path.isdir(source_file):
                            shutil.copytree(source_file, target_file)
                        else:
                            shutil.copy2(source_file, target_file)
                self.temp_line_edit.setText(os.path.join(new_dir))
                cmd.save(os.path.join(new_dir,'load_project.pse'))
        else:
            show_save_dialog(self, temp_path, save=save_file)

def show_save_dialog(self,temp_path,save=1):

    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    dir_path = QFileDialog.getExistingDirectory(self, "Select Directory", "", options=options)
    if dir_path:
        if save:
            if temp_path!=dir_path:
                files=os.listdir(temp_path)
                if os.path.isdir(os.path.join(dir_path, 'NRGSuite_QT_results')):
                    show_popup(self, dir_path, temp_path, 1)
                else:
                    os.mkdir(os.path.join(dir_path, 'NRGSuite_QT_results'))
                    for file_name in files:
                        source_file = os.path.join(temp_path,file_name)
                        target_file = os.path.join(dir_path,'NRGSuite_QT_results',file_name)
                        if '.DS_Store' not in file_name:
                            if os.path.isdir(source_file):
                                shutil.copytree(source_file, target_file)
                            else:
                                shutil.copy2(source_file, target_file)
                    self.temp_line_edit.setText(os.path.join(dir_path,'NRGSuite_QT_results'))
                    cmd.save(os.path.join(dir_path,'NRGSuite_QT_results','load_project.pse'))
            else:
                cmd.save(os.path.join(dir_path,'load_project.pse'))
        if not save:
            files=os.listdir(dir_path)
            for file in files:
                if file=='load_project.pse':
                    cmd.load(os.path.join(dir_path,'load_project.pse'))
                    self.temp_line_edit.setText(dir_path)
                    break
            else:
                show_popup(self,dir_path,temp_path,0)



def refresh_dropdown(dropdown_to_refresh, output_box, filter_for='', no_warning=False,exclude=None, non_group=1, lig=0, add_none=0):
    list_pymol_objects = cmd.get_names('all')
    if non_group and not lig:
        list_pymol_objects_filtered = cmd.get_object_list('all')
        if 'surfaces_results' in list_pymol_objects:
            list_surfaces=cmd.get_object_list('surfaces_results')
            final_list=[]
            if list_surfaces:
                for obj in list_pymol_objects_filtered:
                    if obj not in list_surfaces:
                        final_list.append(obj)
                list_pymol_objects=final_list
        else:
            list_pymol_objects=list_pymol_objects_filtered
    if lig:
        list_pymol_objects = cmd.get_names('selections')
    if add_none:
        list_pymol_objects.append('None')
    if filter_for:
        list_pymol_objects = [x for x in list_pymol_objects if filter_for in x]
    if exclude:
        list_pymol_objects = [item for item in list_pymol_objects if exclude not in item]
    if len(list_pymol_objects) == 0 and no_warning is False:
        output_message(output_box, 'No objects found', 'warning')
    dropdown_to_refresh.clear()
    list_pymol_objects = sorted(list_pymol_objects)
    dropdown_to_refresh.addItems(list_pymol_objects)
    if len(list_pymol_objects) > 0:
        dropdown_to_refresh.setCurrentText(list_pymol_objects[0])




def refresh_folder(folder_path, dropdown_to_refresh):
    folders = next(os.walk(folder_path))[1]
    folders = [item.replace('_', ' ') for item in folders]
    dropdown_to_refresh.clear()
    dropdown_to_refresh.addItems(folders)


def folder_browser(text_window, ligand_set_path, file_extension):
    smile_file_path = QFileDialog.getOpenFileName(None, 'Select a File', ligand_set_path, file_extension)[0]
    if smile_file_path:
        print(smile_file_path)
        text_window.setText(smile_file_path)


def pymol_hide_structures(form):
    list_pymol_objects = cmd.get_names('all')
    list_pymol_objects = [x for x in list_pymol_objects if 'sph' in x]
    if form.button_hide.isChecked():
        if not list_pymol_objects:
            output_message(form.output_box, 'No clefts to hide', 'warning')
        else:
            form.button_hide.setText('Show')
            cmd.hide('everything', ','.join(list_pymol_objects))
    else:
        form.button_hide.setText('Hide')
        cmd.show('surface', ','.join(list_pymol_objects))


def get_mouse_config():
    config_mouse = ''
    try:
        name = cmd.get("button_mode_name")
        if name[0] == '1':
            config_mouse += 'one'
        elif name[0] == '2':
            config_mouse += 'two'
        elif name[0] == '3':
            config_mouse += 'three'
        config_mouse += '_button'
        if name[0] != '1':
            if name[9:] == 'Viewing':
                config_mouse += '_viewing'
            elif name[9:] == 'Editing':
                config_mouse += '_editing'
            elif name[9:] == 'Motions':
                config_mouse += '_motions'
        return config_mouse
    except:
        return 'three_button_viewing'


def read_coords_cleft(cleft_path):
    coords = []
    with open(cleft_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('AT'):
            line = line.split()
            temp_coords = (float(line[6]), float(line[7]), float(line[8]))
            coords.append(temp_coords)
    coords = np.array(coords)
    return lines, coords


def get_residue_info(selection):
    unique_residues = set()
    residue_info = []
    cmd.iterate(selection, 'unique_residues.add((resn, resi, chain))', space={'unique_residues': unique_residues})
    residue_info = [[resname, resn, chain] for resname, resn, chain in unique_residues]
    return residue_info


def create_number_list(length_TotColor, length_TotalColorList):
    if length_TotColor == 1:
        return [0]
    else:
        number_list = []
        modulo = (length_TotalColorList - 1) % (length_TotColor - 1)
        partition = (length_TotalColorList - modulo - 1) / (length_TotColor - 1)
        step_start = 0
        step_end = length_TotalColorList - 1
        for i in range(0, length_TotColor):
            if ((i % 2) == 0):
                number_list.append(step_start)
                step_start = step_start + partition
            else:
                number_list.append(step_end)
                step_end = step_end - partition
        number_list.sort()
        return [int(i) for i in number_list]


def surfaces_enable_buttons(form):
    form.flexaid_retrieve_nrgdock_ligands.setEnabled(True)
    form.surfaces_retreive_flexaid_result.setEnabled(True)

