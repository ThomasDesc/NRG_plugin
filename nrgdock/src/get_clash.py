import numpy as np
from numba import njit


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
