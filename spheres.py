from pymol import cmd
import numpy as np


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


def display_sphere(cleft_object_name, slider):
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
    slider.setMaximum(slider_max)
    slider.setValue(int(slider_max/2))
