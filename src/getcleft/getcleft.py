from PyQt5.QtWidgets import QWidget, QPushButton, QVBoxLayout, QApplication
from PyQt5.QtCore import pyqtSignal, QThread, QObject
import general_functions
import os
from pymol import cmd
# TODO: thread not working

class GetCleftRunner:
    def __init__(self, form):
        super().__init__()
        self.form = form

    def get_parameters(self):
        min_radius = self.form.input_min_radii.text()
        max_radius = self.form.input_max_radii.text()
        resnumc = self.form.input_residue_in_contact.text()
        max_cleft_show = self.form.input_max_cleft_show.text()
        receptor = os.path.basename(self.form.cleft_select_object.currentText())
        parameter_dictionary = {
            'min_radius': min_radius,
            'max_radius': max_radius,
            'resnumc': resnumc,
            'max_cleft_show': max_cleft_show,
            'receptor': receptor
        }
        return parameter_dictionary

    def run_task(self, binary_folder_path, binary_suffix, install_dir):
        print('Starting the task...')
        temp_path = self.form.temp_line_edit.text()
        pymol_object = self.form.cleft_select_object.currentText()
        if pymol_object == '':
            general_functions.output_message(self.form.output_box, 'No object selected', 'warning')
            return
        parameter_dictionary = self.get_parameters()
        self.worker = GetCleftWorker()
        self.worker.message_signal.connect(self.handle_message_signal)
        self.worker.run(binary_folder_path, binary_suffix, temp_path, install_dir, pymol_object, parameter_dictionary)

    def handle_message_signal(self, message, message_type):
        general_functions.output_message(self.form.output_box, message, message_type)


class GetCleftWorker(QThread):
    message_signal = pyqtSignal(str, str)

    def __init__(self):
        super().__init__()

    def get_arg_str(self, getcleft_path, receptor_pdb_path, cleft_save_path, parameter_dictionary):
        min_radius = parameter_dictionary['min_radius']
        max_radius = parameter_dictionary['max_radius']
        resnumc = parameter_dictionary['resnumc']
        max_cleft_show = parameter_dictionary['max_cleft_show']
        receptor = parameter_dictionary['receptor']
        getcleft_output_name = os.path.join(cleft_save_path, receptor)
        arg_string = f'{getcleft_path} -p "{receptor_pdb_path}" -l {min_radius} -u {max_radius} -t {max_cleft_show} -o "{getcleft_output_name}" -s'
        if resnumc != "":
            arg_string += f' -a {resnumc}'
        return arg_string


    def load_show_cleft(self, cleft_save_path, color_list, pymol_object):
        auto_zoom = cmd.get("auto_zoom")
        cmd.set("auto_zoom", 0)
        all_files = os.listdir(cleft_save_path)
        sph_file_list = []
        for filename in all_files:
            if filename.find('sph') != -1:
                sph_file_list.append({'path': os.path.join(cleft_save_path, filename),
                                      'name': filename.split('.')[0]})
        sph_file_list = sorted(sph_file_list, key=lambda d: d['name'])
        if len(sph_file_list) == 0:
            self.message_signal.emit('No clefts were found', 'warning')
        number_color_list = general_functions.create_number_list(len(sph_file_list), len(color_list))
        for cleft_counter, Cleft in enumerate(sph_file_list):
            try:
                cmd.load(Cleft['path'], Cleft['name'], state=1)
                cmd.group('Clefts',Cleft['name'])
                cmd.hide('everything', Cleft['name'])
                if cleft_counter >= len(color_list):
                    cmd.color('grey50', Cleft['name'])
                else:
                    cmd.color(color_list[number_color_list[cleft_counter]], Cleft['name'])
                cmd.show('surface', Cleft['name'])
            except:
                print(f"ERROR: Failed to load cleft object  {Cleft['name']}")
                continue
        cmd.zoom(pymol_object)
        cmd.refresh()
        cmd.set("auto_zoom", auto_zoom)

    def load_color_list(self, color_list_path):
        with open(color_list_path, 'r') as file:
            color_list = [line.strip() for line in file]
        return color_list


    def run(self, binary_folder_path, binary_suffix, temp_path, nrgsuite_base_path, pymol_object, parameter_dictionary):
        getcleft_binary_path = os.path.join(binary_folder_path, f'GetCleft{binary_suffix}')
        getcleft_output_path = os.path.join(temp_path, 'GetCleft')
        cleft_save_path = os.path.join(getcleft_output_path, 'Clefts')
        color_list_path = os.path.join(nrgsuite_base_path, 'deps', 'getcleft', 'color_list.txt')
        color_list = self.load_color_list(color_list_path)
        os.makedirs(getcleft_output_path, exist_ok=True)
        os.makedirs(cleft_save_path, exist_ok=True)
        object_save_path = os.path.join(getcleft_output_path, 'tmp.pdb')
        cmd.save(object_save_path, pymol_object)

        getcleft_command = self.get_arg_str(getcleft_binary_path, object_save_path, cleft_save_path, parameter_dictionary)
        print(getcleft_command)
        self.message_signal.emit('Running GetCleft...', 'valid')
        import time
        time.sleep(2)
        os.system(getcleft_command)
        self.load_show_cleft(cleft_save_path, color_list, pymol_object)
        self.message_signal.emit('Done GetClef' , 'valid')