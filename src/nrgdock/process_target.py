import os
from time import perf_counter

import numpy as np
from numba import njit
import sys
install_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(install_dir)
from src.nrgdock.main_processed_target import get_params_dict
import shutil
# from process_ligands import preprocess_ligands_one_target as preprocess_ligand
import timeit
import concurrent.futures
# from itertools import repeat
import argparse
from pathlib import Path
import json
import time

def load_rad_dict(filepath):
    with open(filepath, 'r') as file:
        loaded_atom_data = json.load(file)
    return loaded_atom_data


def get_radius_number(letter_type, rad_dict):
    if isinstance(letter_type, int):
        return list(rad_dict.values())[letter_type-1][1]
    else:
        letter_type = letter_type.upper().replace(' ', '')
        try:
            atm_info = [rad_dict[letter_type]['type_number'], rad_dict[letter_type]['radius']]
        except KeyError:
            atm_info = [39, 2.00]
        atm_type_num = atm_info[0]
        atm_rad = atm_info[1]
        return atm_type_num, atm_rad


def load_atoms_mol2(filename, rad_dict):
    coord_start = 0
    atoms_xyz = []
    atoms_numbers = []
    atoms_types = []
    atoms_radius = []
    with open(filename) as f:
        for line in f:
            if line.strip():
                line = line.split()
                if line[0] == '@<TRIPOS>ATOM':
                    coord_start = 1
                if line[0][0] == '@' and coord_start == 1 and line[0] != '@<TRIPOS>ATOM':
                    break
                if coord_start == 1 and line[0] != '@<TRIPOS>ATOM' and line[5].split(".")[0] != 'H':
                    atoms_xyz.append([float(line[2]), float(line[3]), float(line[4])])
                    atoms_numbers.append(int(line[0]))
                    atoms_type, atom_radius = get_radius_number(line[5], rad_dict)
                    atoms_types.append(atoms_type)
                    atoms_radius.append(atom_radius)
    atoms_xyz = np.array(atoms_xyz, dtype=np.float32)
    atoms_radius = np.array(atoms_radius, dtype=np.float32)
    atoms_numbers = np.array(atoms_numbers, dtype=np.int32)
    atoms_types = np.array(atoms_types, dtype=np.int32)
    atoms_types_sorted = atoms_types[atoms_numbers.argsort()]

    return atoms_xyz, atoms_types_sorted, atoms_radius


def find_cleft_file(target_folder, ligand_name=None):
    print('Running GetCleft')
    receptor_pdb_path = os.path.join(target_folder, "receptor.pdb")
    if not os.path.exists(receptor_pdb_path):
        exit('Could not find receptor.pdb file')
    get_cleft_folder = os.path.join(target_folder, "get_cleft")
    if not os.path.isdir(get_cleft_folder):
        os.mkdir(get_cleft_folder)
    get_cleft_path = ''
    output = os.path.join(get_cleft_folder, "bd_site")
    get_cleft_command = f'{get_cleft_path} -p "{receptor_pdb_path}" -o "{output}" -l 1.75 -u 3.6 -s -t 1'
    if ligand_name:
        get_cleft_command += f' -a {ligand_name}-'
    # print(get_cleft_command)
    os.system(get_cleft_command)
    for file in os.listdir(get_cleft_folder):
        if file.find("_sph_") != -1:
            return
    exit('Could not find binding site file')



def find_cleft_file_simple(target_folder):
    get_cleft_folder = os.path.join(target_folder, "get_cleft")
    for file in os.listdir(get_cleft_folder):
        if file.find("_sph_") != -1:
            return os.path.join(get_cleft_folder, file)


def load_binding_site_pdb(binding_site):
    atom_coord_type_list = []
    with open(binding_site) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('ATOM'):
                line = line.strip()
                temp_array = [float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[60:66])]  # sphere radius
                atom_coord_type_list.append(temp_array)
    return atom_coord_type_list


def load_binding_site_grid(dot_division, a, target_path, padding):
    # find min and max coords for box as well as remove or add sphere radius
    x = [np.min(a[:, 0]) - a[np.argmin(a[:, 0]), 3] - padding, np.max(a[:, 0]) + a[np.argmax(a[:, 0]), 3] + padding]
    y = [np.min(a[:, 1]) - a[np.argmin(a[:, 1]), 3] - padding, np.max(a[:, 1]) + a[np.argmax(a[:, 1]), 3] + padding]
    z = [np.min(a[:, 2]) - a[np.argmin(a[:, 2]), 3] - padding, np.max(a[:, 2]) + a[np.argmax(a[:, 2]), 3] + padding]
    bd_site_box = []
    for x_coord in x:
        for y_coord in y:
            for z_coord in z:
                bd_site_box.append([x_coord, y_coord, z_coord])
    # write_test(bd_site_box, f"{os.path.split(target_path)[-1]}_bd_site_box", f'./temp/ligand_poses', None, None)
    x_range = np.arange(x[0], x[1], dot_division)
    y_range = np.arange(y[0], y[1], dot_division)
    z_range = np.arange(z[0], z[1], dot_division)
    np.save(os.path.join(target_path, 'preprocessing_files', "range_array"), np.stack((x, y, z)))
    return x_range, y_range, z_range


def build_3d_cube_grid(params, target_atoms_xyz, atoms_radius, cw_factor=1):
    water_radius = params['WATER_RADIUS']
    grid_placeholder = params['GRID_PLACEHOLDER']
    max_rad = np.amax(atoms_radius, axis=0)
    cell_width = 2 * (max_rad + water_radius)
    #cell_width = 4
    # max_xyz = np.amax(target_atoms_xyz, axis=0)
    max_xyz = np.zeros(3)
    max_xyz[0] = np.max(target_atoms_xyz[:, 0]) + cell_width*cw_factor
    max_xyz[1] = np.max(target_atoms_xyz[:, 1]) + cell_width*cw_factor
    max_xyz[2] = np.max(target_atoms_xyz[:, 2]) + cell_width*cw_factor
    # min_xyz = np.amin(target_atoms_xyz, axis=0)
    min_xyz = np.zeros(3)
    min_xyz[0] = np.min(target_atoms_xyz[:, 0]) - cell_width*cw_factor
    min_xyz[1] = np.min(target_atoms_xyz[:, 1]) - cell_width*cw_factor
    min_xyz[2] = np.min(target_atoms_xyz[:, 2]) - cell_width*cw_factor
    lengths = ((max_xyz - min_xyz) / cell_width).astype(np.int32) + 1
    temp_grid = []
    for i in range(lengths[0]):
        temp_grid.append([])
        for j in range(lengths[1]):
            temp_grid[i].append([])
            for k in range(lengths[2]):
                temp_grid[i][j].append([])
    for i, row in enumerate(target_atoms_xyz):
        grid_indices = ((row[:3] - min_xyz) / cell_width).astype(np.int32)
        temp_grid[grid_indices[0]][grid_indices[1]][grid_indices[2]].append(i)
    max_cell_len = 0
    for row in temp_grid:
        for col in row:
            for cell in col:
                n = len(cell)
                if n > max_cell_len:
                    max_cell_len = n
    grid = np.full((lengths[0], lengths[1], lengths[2], max_cell_len), grid_placeholder, dtype=np.int32)
    for i in range(lengths[0]):
        for j in range(lengths[1]):
            for k in range(lengths[2]):
                for x, v in enumerate(temp_grid[i][j][k]):
                    grid[i][j][k][x] = v
    return grid, min_xyz, cell_width, max_xyz


def get_radius_list_from_nums(rad_dict, target_path, constant_radius):
    if constant_radius is not None:
        atom_number_list = np.arange(1, 40, 1, dtype=np.int32)
        num_rad_list = np.zeros((len(atom_number_list), 2), dtype=np.float32)
        for a, number in enumerate(atom_number_list):
            num_rad_list[a][0] = number
            num_rad_list[a][1] = constant_radius
    else:
        num_rad_list = [[details['type_number'], details['radius']] for details in rad_dict.values()]
        num_rad_list = np.array(sorted(num_rad_list, key=lambda x: x[0]))
    np.save(os.path.join(target_path, 'preprocessing_files', "minimum_radius_list"), num_rad_list)
    return num_rad_list


def prepare_preprocess_output(path_to_target, params_dict, config_file_path):
    numpy_output_path = os.path.join(path_to_target, 'preprocessing_files')
    if not os.path.isdir(numpy_output_path):
        os.mkdir(numpy_output_path)

    config_file_name = f"config_{params_dict['PRELOAD_GRID_DISTANCE']}"
    config_output = os.path.join(numpy_output_path, config_file_name + ".txt")
    if os.path.isfile(config_output):
        config_file_number = 1
        while os.path.isfile(os.path.join(numpy_output_path, config_file_name + f'_{config_file_number}.txt')):
            config_file_number += 1
        config_output = os.path.join(numpy_output_path, config_file_name + f'_{config_file_number}.txt')
    shutil.copyfile(config_file_path, config_output)
    return numpy_output_path


def load_ligand_test_dots(params, atom_coord_type_list):
    """ This will load the binding site spheres, find the extreme values for x, y, z, then build a grid of points
    and remove the ones that clash with the target atoms.
    """
    dot_division = params['DOT_DIVISION']
    grid_coords = []
    a = np.array(atom_coord_type_list)
    # find min and max coords for box as well as remove or add sphere radius
    x = [round(np.min(a[:, 0]) - a[np.argmin(a[:, 0]), 3], 3), round(np.max(a[:, 0]) + a[np.argmax(a[:, 0]), 3], 3)]
    y = [round(np.min(a[:, 1]) - a[np.argmin(a[:, 1]), 3], 3), round(np.max(a[:, 1]) + a[np.argmax(a[:, 1]), 3], 3)]
    z = [round(np.min(a[:, 2]) - a[np.argmin(a[:, 2]), 3], 3), round(np.max(a[:, 2]) + a[np.argmax(a[:, 2]), 3], 3)]

    for dot_x in np.arange(x[0], x[1], dot_division):
        for dot_y in np.arange(y[0], y[1], dot_division):
            for dot_z in np.arange(z[0], z[1], dot_division):
                coords = np.array([round(dot_x, 3), round(dot_y, 3), round(dot_z, 3)])
                for row in a:
                    distance = np.linalg.norm(coords - row[:3])
                    if distance < row[3]:
                        grid_coords.append(coords)
                        break
    return grid_coords


def clean_binding_site_grid(params, target_grid, binding_site_grid, min_xyz, cell_width, target_atoms_xyz, target_path):
    index = []
    dot_division = params['DOT_DIVISION']
    for a, point in enumerate(binding_site_grid):
        grid_index = ((point - min_xyz) / cell_width).astype(np.int32)
        for i_offset in [-1, 0, 1]:
            for j_offset in [-1, 0, 1]:
                for k_offset in [-1, 0, 1]:
                    i = i_offset + grid_index[0]
                    j = j_offset + grid_index[1]
                    k = k_offset + grid_index[2]
                    if i < len(target_grid) and j < len(target_grid[0]) and k < len(target_grid[0][0]):
                        for neighbour in target_grid[i][j][k]:
                            if neighbour == -1:
                                break
                            else:
                                dist = np.sqrt((target_atoms_xyz[neighbour][0] - point[0]) ** 2 +
                                               (target_atoms_xyz[neighbour][1] - point[1]) ** 2 +
                                               (target_atoms_xyz[neighbour][2] - point[2]) ** 2)
                                if dist <= 2.0:
                                    index.append(a)

                                    break
    cleaned_binding_site_grid = np.delete(binding_site_grid, index, 0)
    # write_test(cleaned_binding_site_grid, "cleaned grid", f'./temp/ligand_poses/', None, None)
    np.save(os.path.join(target_path, 'preprocessing_files', f"cleaned_grid_{dot_division}"), cleaned_binding_site_grid)


def write_dots_for_3d_cube(cell_width, min_xyz, max_xyz):
    coord_list = []
    name_list = []
    for x in np.arange(min_xyz[0], max_xyz[0], cell_width):
        for y in np.arange(min_xyz[1], max_xyz[1], cell_width):
            for z in np.arange(min_xyz[2], max_xyz[2], cell_width):
                coord_list.append([x, y, z])
                if not name_list:
                    name_list = ["O"]
                elif name_list == ["O"]:
                    name_list.append("N")
                else:
                    name_list.append("C")
    # write_test(coord_list, "binding_grid_test", r"C:\Users\thoma\Desktop\rotation_fix", name_list, None)


@njit
def get_cf_list(target_grid, atom_type_range, target_atom_types, energy_matrix, number_types):
    # test_dot = [[0, 0, 9]]
    # test_dot = np.array(test_dot)
    # test_dot = (test_dot * cell_width + min_xyz) - (cell_width/2)
    # write_test(test_dot, "test_dot", r"C:\Users\thoma\Desktop\rotation_fix", None, None)
    target_grid_x = len(target_grid)
    target_grid_y = len(target_grid[0])
    target_grid_z = len(target_grid[0][0])
    result_array = np.zeros((target_grid_x, target_grid_y, target_grid_z, number_types))
    for x in np.arange(0, target_grid_x):
        # print(f"{x}/{target_grid_x}")
        for y in np.arange(0, target_grid_y):
            for z in np.arange(0, target_grid_z):
                grid_index = [x, y, z]
                for counter, atom_type in enumerate(atom_type_range):
                    cf = 0.0
                    # if atom_type != 3:
                    for i_offset in [-1, 0, 1]:
                        for j_offset in [-1, 0, 1]:
                            for k_offset in [-1, 0, 1]:
                                i = i_offset + grid_index[0]
                                j = j_offset + grid_index[1]
                                k = k_offset + grid_index[2]
                                if 0 < i < len(target_grid) and 0 < j < len(target_grid[0]) and 0 < k < len(target_grid[0][0]):
                                    if target_grid[i][j][k][0] != -1:
                                        for neighbour in target_grid[i][j][k]:
                                            if neighbour == -1:
                                                break
                                            else:
                                                # Normal program:
                                                type_1 = atom_type
                                                type_2 = target_atom_types[neighbour]
                                                # Randomly chosen atom type
                                                # type_1 = np.random.randint(1, 41)
                                                # type_2 = np.random.randint(1, 41)
                                                energy_value = energy_matrix[type_1][type_2]
                                                cf += energy_value
                    result_array[x, y, z, counter] = cf
    return result_array


@njit
def get_clash_per_dot(x_range, y_range, z_range, target_grid, min_xyz, cell_width, target_atoms_xyz, max_size_array,
                      total_coords, use_clash):
    clash_list = np.zeros((max_size_array[0], max_size_array[1], max_size_array[2]), dtype=np.bool_)
    for a, x_value in enumerate(x_range):
        for b, y_value in enumerate(y_range):
            for c, z_value in enumerate(z_range):
                ligand_atom = np.array([x_value, y_value, z_value], dtype=np.float32)
                clash_list[a][b][c] = get_clash(ligand_atom, target_grid, min_xyz, cell_width, target_atoms_xyz,
                                                use_clash)
    return clash_list


def save_files(min_xyz, cell_width, preprocessed_file_path):
    np.save(os.path.join(preprocessed_file_path, f"min_xyz"), min_xyz)
    np.save(os.path.join(preprocessed_file_path, f"cell_width"), cell_width)


@njit
def get_clash(ligand_atom_coord, target_grid, min_xyz, cell_width, target_atoms_xyz, use_clash):
    # TODO: test if clashes per radius associated to each type is better
    # array_to_fill = np.full(type_to_test[:, 1].size, default_cf, dtype=np.float32)
    clash = False
    # for counter, radius in enumerate(type_to_test):
    grid_index = ((ligand_atom_coord - min_xyz) / cell_width).astype(np.int32)
    for i_offset in [-1, 0, 1]:
        for j_offset in [-1, 0, 1]:
            for k_offset in [-1, 0, 1]:
                i = i_offset + grid_index[0]
                j = j_offset + grid_index[1]
                k = k_offset + grid_index[2]
                if 0 <= i < len(target_grid) and 0 <= j < len(target_grid[0]) and 0 <= k < len(target_grid[0][0]):
                    if target_grid[i][j][k][0] != -1:
                        for neighbour in target_grid[i][j][k]:
                            if neighbour == -1:
                                break
                            else:
                                if use_clash is not None:
                                    dist = np.linalg.norm(target_atoms_xyz[neighbour] - ligand_atom_coord)
                                    if dist <= 2.0:
                                        clash = True
                                        return clash
    return clash


def preprocess_one_target(target, main_path, params_dict, energy_matrix, time_start, overwrite, run_getcleft,
                          config_file_path, deps_path, ligand_name=None, verbose=True):
    if verbose:
        print(f"Target: {target}")
    target_path = str(os.path.join(main_path, target))
    receptor_path = os.path.join(target_path, "receptor.mol2")
    print(receptor_path)
    if not os.path.isfile(receptor_path):
        exit(f"Could not find receptor.mol2 file at path: {receptor_path}")

    # ####################### GET BINDING SITE #######################
    if verbose:
        print("Finding Binding site")
    if run_getcleft:
        find_cleft_file(target_path, ligand_name)
    
    # ####################### BUILDING GRIDS #######################
    if verbose:
        print('Building grids')
    preprocessed_file_path = prepare_preprocess_output(target_path, params_dict, config_file_path)
    rad_dict = load_rad_dict(os.path.join(deps_path, "atom_type_radius.json"))
    use_clash = params_dict["USE_CLASH"]
    use_constant_radius = params_dict['CONSTANT_RADIUS']
    grid_distance = params_dict['PRELOAD_GRID_DISTANCE']
    precalculated_cf_box_padding = params_dict["CF_PRECALC_BOX_PADDING"]
    dot_division = params_dict['DOT_DIVISION']
    constant_radius = None
    if use_constant_radius:
        constant_radius = params_dict['ATOM_RADIUS']
    binding_site = find_cleft_file_simple(target_path)
    radius_list = get_radius_list_from_nums(rad_dict, target_path, constant_radius)
    target_atoms_xyz, target_atoms_types, atoms_radius = load_atoms_mol2(receptor_path, rad_dict)
    target_grid, min_xyz, cell_width, max_xyz = build_3d_cube_grid(params_dict, target_atoms_xyz, atoms_radius)
    # write_dots_for_3d_cube(cell_width, min_xyz, max_xyz)
    save_files(min_xyz, cell_width, preprocessed_file_path)
    binding_site_spheres = load_binding_site_pdb(binding_site)
    x_range, y_range, z_range = load_binding_site_grid(grid_distance, np.array(binding_site_spheres), target_path,
                                                       precalculated_cf_box_padding)

    # ####################### PRECALCULATE CF #######################

    if verbose:
        print('Precalculating CF')
    if not os.path.isfile(
            os.path.join(preprocessed_file_path, f"cf_list_{grid_distance}.npy")) or overwrite:
        max_size_array = np.array([len(x_range), len(y_range), len(z_range), len(radius_list)], dtype=np.int32)
        total_coords = max_size_array[0] * max_size_array[1] * max_size_array[2]
        atom_type_range = np.arange(1, len(radius_list)+1)

        if use_clash:
            print('Getting clashes')
            clash_list = get_clash_per_dot(x_range, y_range, z_range, target_grid, min_xyz, cell_width,
                                           target_atoms_xyz, max_size_array, total_coords, constant_radius)
            np.save(os.path.join(preprocessed_file_path, f"clash_list_{grid_distance}"), clash_list)
        cfs_list = get_cf_list(target_grid, atom_type_range, target_atoms_types, energy_matrix, len(radius_list))
        np.save(os.path.join(preprocessed_file_path, f"cf_list_{grid_distance}"), cfs_list)
    else:
        print(f"CF already precalculated for {grid_distance} A between atoms")

    # ####################### GENERATE AND CLEAN BINDING SITE GRID #######################

    if not os.path.isfile(os.path.join(preprocessed_file_path, f"cleaned_grid_{dot_division}.npy")) or overwrite:
        original_grid = load_ligand_test_dots(params_dict, binding_site_spheres)
        clean_binding_site_grid(params_dict, target_grid, original_grid, min_xyz, cell_width, target_atoms_xyz,
                                target_path)
    else:
        print(f"The file for binding site dots at {dot_division} A distance already exists")

    if verbose:
        total_run_time = timeit.default_timer() - time_start
        if total_run_time > 60.0:
            total_run_time /= 60
            print(f"{target}: {total_run_time} minutes to run")
        else:
            print(f"{target}: {total_run_time} seconds to run")


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", '--path_to_targets', required=True, type=str,
                        help='Path to folder containing target folders')
    parser.add_argument('-t', '--specific_target', type=str,
                        help='Specify target(s) to analyse. If multiple: separate with comma no space')
    parser.add_argument('-o', '--overwrite', action='store_true',
                        help='specify if you want to overwrite pre existing files')
    parser.add_argument("-c", '--run_getcleft', action='store_true',
                                   help='Specify if you want to run GetCleft')
    parser.add_argument("-l", '--ligand_name', type=str,
                        help='Ligand for GetCleft format (no spaces ex: RES999A): residue, number, chain')
    parser.add_argument("-d", '--deps_path', type=str,
                        help='path to deps')
    args = parser.parse_args()
    if args.ligand_name and args.run_getcleft is None:
        parser.error("-l (--ligand_name) requires -c (--run_getcleft).")

    path_to_targets = os.path.abspath(args.path_to_targets)
    if args.specific_target is None:
        target_list = next(os.walk(path_to_targets))[1]
    else:
        target_list = args.specific_target.split(',')
    if args.deps_path is None:
        deps_path=None
    else:
        deps_path = args.deps_path
    target_list = sorted(target_list)
    main(path_to_targets, target_list, overwrite=args.overwrite, run_getcleft=args.run_getcleft, ligand_name=args.ligand_name, deps_path=deps_path)


def main(path_to_targets, target_list, overwrite=False, run_getcleft=False, ligand_name=None, deps_path=os.path.join('.', 'deps')):
    # print('WARNING randomly assigning atom type when calculating cf')
    root_software_path = Path(__file__).resolve().parents[1]
    os.chdir(root_software_path)
    time_start = timeit.default_timer()
    print(deps_path)
    config_file = os.path.join(deps_path, "config.txt")
    params_dict = get_params_dict(config_file)
    matrix_name = params_dict['PRECALC_MATRIX_NAME']
    matrix_path = os.path.join(deps_path, fr"matrix/{matrix_name}.npy")
    energy_matrix = np.load(matrix_path)
    verbose = True
    if len(target_list) > 1:
        verbose = False
        print('Preprocessing targets: ', target_list)
    else:
        print('Preprocessing target: ', target_list[0])
    preprocess_one_target(target_list[0], path_to_targets, params_dict, energy_matrix, time_start, overwrite,
                          run_getcleft, config_file, deps_path, ligand_name, verbose)
    # with concurrent.futures.ProcessPoolExecutor(max_workers=32) as executor:
    #     executor.map(preprocess_one_target, target_list, repeat(path_to_targets), repeat(params_dict),
    #                  repeat(energy_matrix), repeat(time_start), repeat(overwrite), repeat(run_getcleft),
    #                  repeat(ligand_name), repeat(verbose))


if __name__ == "__main__":
    get_args()
