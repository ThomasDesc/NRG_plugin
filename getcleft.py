import os
import subprocess
from pymol import cmd


def get_arg_str(form, getcleft_path, object_path, cleft_save_path):
    min_radius = form.input_min_radii.text()
    max_radius = form.input_max_radii.text()
    resnumc = form.input_residue_in_contact.text()
    max_cleft_show = form.input_max_cleft_show.text()
    print('getcleft output path: ', cleft_save_path)
    getcleft_output_name = os.path.join(cleft_save_path, 'receptor')
    arg_string = f'{getcleft_path} -p "{object_path}" -l {min_radius} -u {max_radius} -t {max_cleft_show} -o "{getcleft_output_name}" -s'
    # TODO: Check if correct
    if resnumc != "":
        arg_string += f' -a {resnumc}-'
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


def load_show_cleft(cleft_save_path, color_list):
    auto_zoom = cmd.get("auto_zoom")
    cmd.set("auto_zoom", 0)
    all_files = os.listdir(cleft_save_path)
    sph_file_list = []
    for filename in all_files:
        if filename.find('sph') != -1:
            sph_file_list.append({'path': os.path.join(cleft_save_path, filename),
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


def run_getcleft(form, getcleft_path, getcleft_output_path, cleft_save_path, color_list):
    pymol_object = form.cleft_select_object.currentText()
    if pymol_object != '':
        object_save_path = os.path.join(getcleft_output_path, 'tmp.pdb')
        cmd.save(object_save_path, pymol_object)
    else:
        print('No object selected')
        return
    getcleft_command = get_arg_str(form, getcleft_path, object_save_path, cleft_save_path)
    print(getcleft_command)
    subprocess.run(getcleft_command, shell=True)
    load_show_cleft(cleft_save_path, color_list)