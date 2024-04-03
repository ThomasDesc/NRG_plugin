from pymol import cmd
import numpy as np


def get_center(coordinate_array):
    length = coordinate_array.shape[0]
    sum_x = np.sum(coordinate_array[:, 0])
    sum_y = np.sum(coordinate_array[:, 1])
    sum_z = np.sum(coordinate_array[:, 2])
    return [sum_x / length, sum_y / length, sum_z / length]


def display_sphere(cleft_object_name):
    if cleft_object_name == '':
        print('No cleft object selected in step 1')
        return
    else:
        sphere_name = 'SPHERE_1'
        cleft_coordinates = cmd.get_coords(cleft_object_name, state=1)
        cmd.pseudoatom(sphere_name,
                       pos=get_center(cleft_coordinates),
                       vdw=self.SphereView.Radius,
                       state=self.State)
        cmd.refresh()

        cmd.color('oxygen', sphere_name)
        cmd.hide('everything', sphere_name)
        cmd.show('spheres', sphere_name)
        cmd.refresh()