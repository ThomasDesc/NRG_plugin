from pymol import cmd
import numpy as np
import wizard
from general_functions import read_coords_cleft
import os


def get_center(cleft_coordinates):
    length = cleft_coordinates.shape[0]
    sum_x = np.sum(cleft_coordinates[:, 0])
    sum_y = np.sum(cleft_coordinates[:, 1])
    sum_z = np.sum(cleft_coordinates[:, 2])
    return [sum_x / length, sum_y / length, sum_z / length]


def get_max_coords(cleft_coordinates, center):
    distances = np.linalg.norm(cleft_coordinates - center, axis=1)
    radius = np.max(distances)
    # max_diameter += 6
    return radius


def display_sphere(cleft_object_name, slider, partition_sphere_select):
    if cleft_object_name == '':
        print('No cleft object selected in step 1')
        return
    else:
        sphere_name = 'SPHERE_1'
        cleft_coordinates = np.array(cmd.get_model(cleft_object_name, 1).get_coord_list())
        center_coordinate = get_center(cleft_coordinates)
        max_vdw = get_max_coords(cleft_coordinates, center_coordinate)
        cmd.pseudoatom(sphere_name,
                       pos=center_coordinate,
                       vdw=max_vdw,
                       state=1)
        cmd.refresh()

        cmd.color('oxygen', sphere_name)
        cmd.hide('everything', sphere_name)
        cmd.show('spheres', sphere_name)
        cmd.refresh()
    slider_max = max_vdw*100
    slider.setEnabled(True)
    slider.setMaximum(slider_max)
    slider.setValue(slider_max)
    slider.setSingleStep(10)
    partition_sphere_select.addItem(sphere_name)
    partition_sphere_select.setCurrentText(sphere_name)


def resize_sphere(sphere_name, slider_value):
    cmd.alter(sphere_name, 'vdw=' + str(slider_value/100))
    cmd.rebuild(sphere_name)
    cmd.refresh()


def move_sphere(cleft_name):
    cmd.window('hide')
    cmd.window('show')
    wiz = wizard.Sphere()
    cmd.set_wizard(wiz)


def crop_cleft(object_name, sphere_vdw, cleft_save_path, cleft_name):
    sphere_coords = np.array(cmd.get_model(object_name).atom[0].coord)
    min_coords = np.array([sphere_coords-sphere_vdw])
    max_coords = np.array([sphere_coords+sphere_vdw])
    lines, partition_coords = read_coords_cleft(os.path.join(cleft_save_path, cleft_name + '.pdb'))
    indices = np.where((partition_coords > min_coords).all(axis=1) & (partition_coords < max_coords).all(axis=1))[0] + 1
    lines_to_output = [lines[0]] + [lines[i] for i in indices]
    partition_path = os.path.join(cleft_save_path, cleft_name + '_cropped.pdb')
    with open(partition_path, 'w') as f:
        f.writelines(lines_to_output)
    partition_name = os.path.basename(partition_path).split('.')[0]
    cmd.set("auto_zoom", 0)
    cmd.load(partition_path, format='pdb')
    cmd.hide('everything', partition_name)
    cmd.show('surface', partition_name)
    cmd.color('grey60', partition_name)
    cmd.disable(f'{object_name} {cleft_name}')
    cmd.enable(cleft_name)
    cmd.refresh()

