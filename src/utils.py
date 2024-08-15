from typing import List, Tuple, Set, Dict
import os

import json
import numpy as np

JSON_DIR = "jsons"


def wrap_to_cell(w: int, x_width: int):
    """convert a 1d index to a 2d index in a square grid of width x_width"""
    return w % x_width, (w // x_width)


def unwrap_to_1d(x: int, y: int, x_width: int) -> int:
    """convert a 2d index to a 1d index in a square grid of width x_width"""
    return x + y*x_width


def get_neighbors(x_width: int, y_height: int) -> List[Tuple[int, int]]:
    """return a list of tuples of neighboring cells"""
    neighbors = []
    for i in range(x_width*y_height):
        if i % x_width != x_width-1:
            neighbors.append((i, i+1))
        if i // x_width != x_width-1:
            neighbors.append((i, i+x_width))
    return neighbors


def reduce_tile_list(
        tile_solutions: List[List[int]],
        equivalent_cell_indices: List[List[Tuple[int]]] = None,
        torso_length: int = 3,
        arm_length: int = 2,
        leg_length: int = 3
) -> Set[Tuple[int]]:
    """Remove duplicates from the tile solutions.

    If equivalent cell indices are given explicitely, length params are ignored.

    Args:
        tile_solutions: list of solutions, each solution is a list of integers (in terms of w coordinates).
        equivalent_cell_indices: list of lists of tuples, inside lists are sequences of equivalent cell indices.
        torso_length: int, the length of the torso. Can be used if equivalent_cell_indices are not explicitly given.
        arm_length: int, the length of the arm. Can be used if equivalent_cell_indices are not explicitly given.
        leg_length: int, the length of the leg. Can be used if equivalent_cell_indices are not explicitly given.
    """
    # two placements are equivalent if the w-coord of the cells on index (4,5) and (6,7) are swapped.
    # and/or the w-coord of the cells on index (8,9,10) and (11,12,13) are swapped
    # always sort them to make them unique. Put the smaller first position on the first limb.
    if equivalent_cell_indices is None:
        equivalent_cell_indices = get_equivalent_indices_man_tile(
            torso_length=torso_length,
            arm_length=arm_length,
            leg_length=leg_length
        )

    def _hash(pos_list: List[int]) -> Tuple[int]:
        list_copy = pos_list.copy()
        for equiv in equivalent_cell_indices:
            left = [pos_list[j] for j in equiv[0]]
            right = [pos_list[j] for j in equiv[1]]
            sorted_limbs = sorted([left,right])
            for j in [0,1]:
                for k, index in enumerate(equiv[j]):
                    list_copy[index] = sorted_limbs[j][k]
        return tuple(list_copy)

    return set(map(_hash, tile_solutions))


def get_equivalent_indices_man_tile(
        torso_length: int = 3,
        arm_length: int = 2,
        leg_length: int = 3
) -> List[List[Tuple[int]]]:
    """Return the equivalent indices (under mirror symmetry) of cells in the man tile"""
    indices = []
    pointer = torso_length + 1
    for limb_length in [arm_length, leg_length]:
        limb_indices = []
        for _ in [0, 1]:
            limb_indices.append(tuple(range(pointer, pointer + limb_length)))
            pointer += limb_length
        indices.append(limb_indices)
    return indices


def show_chunk(
        y_width: int,
        x_width: int,
        pts: Set[Tuple[int, int]],
        optional_pts: Set[Tuple[int, int]],
        chunk_history: List[Set[Tuple[int, int]]],
):
    """Visualize a chunk of the grid with points and optional points"""
    x_width = max([s[0] + 1 for s in pts.union(optional_pts)] + [x_width])

    chunk_array = [["|"]+[" " for __ in range(x_width)]+["|\n"] for _ in range(y_width)]
    for x,y in pts:
        chunk_array[y][x+1] = "■"
    for x,y in optional_pts:
        chunk_array[y][x+1] = "□"
    for h, old_chunk in enumerate(chunk_history):
        symbol = [" ", "⢀", "⣀", "⣠", "⣤", "⣴", "⣶", "⣾", "⣿"][h % 9]
        for x,y in old_chunk:
            chunk_array[y][x+1] = symbol
    top_bdy = "┌" + "─" * x_width + "┐\n"
    bottom_bdy = "└" + "─" * x_width + "┘"
    full_string = top_bdy + "".join(["".join(row) for row in chunk_array]) + bottom_bdy
    print(full_string)


def read_json_tiles(filename) -> List[np.array]:
    """read a json file with tiles and return a list of numpy arrays"""
    with open(f"{JSON_DIR}/{filename}") as infile:
        data = json.load(infile)
    return [np.array(tile, dtype=int) for tile in data["tiles"]]


def three_pts_collinear(pt1: tuple, pt2: tuple, pt3: tuple) -> bool:
    """check if three points are collinear"""
    return any(pt1[j] == pt2[j] == pt3[j] for j in [0,1])


def get_strands(juxt_list: List[Tuple[int]]) -> List[List[int]]:
    """return a list of strands of connected cells"""
    # identify all pts with more than 2 outgoing connections
    connex: Dict[int,List[int]] = {}
    for juxt in juxt_list:
        for h in [0,1]:
            connex.setdefault(juxt[h], []).append(juxt[1-h])

    start_pt = None
    for pt, connections in connex.items():
        if len(connections) == 1:
            start_pt = pt
            break
    if start_pt is None:
        raise ValueError("Not a tree")

    def strand_web(
            pt0: int,
            pt1: int,
    ) -> List[List[int]]:
        strand = [pt0, pt1]
        while len(connex[pt1]) == 2:
            a, b = connex[pt1]
            pt1 = a if b == strand[-2] else b
            strand.append(pt1)
        remaining_pts = [pt for pt in connex[pt1] if pt != strand[-2]]
        return [strand]+sum([strand_web(pt1, pt) for pt in remaining_pts],[])

    return strand_web(start_pt, connex[start_pt][0])


def _rot0(x, y):
    return x, y


def _rot90(x, y):
    return y, x


def _rot180(x, y):
    return x, -y


def _rot270(x, y):
    return -y, x


def rotate_eyes(sol: np.array, eye_index: int, pointing_towards: int, eye_coords: List[Tuple[float]]) -> List[Tuple[float,float]]:
    dxdy = sol[pointing_towards] - sol[eye_index]

    if dxdy[0] == 0:  # vertical
        if dxdy[1] < 0:
            f = _rot0
        else:
            f = _rot180
    elif dxdy[1] == 0:  # horizontal
        if dxdy[0] < 0:
            f = _rot90
        else:
            f = _rot270
    else:
        raise ValueError("Eyes not aligned with axis")
    return [tuple(sol[eye_index] + f(*eye_coord)) for eye_coord in eye_coords]


class RecurseConfig:
    """Just to store configuration. This object is shared by all recursion levels."""

    def __init__(
            self,
            dx_per_iter: int,
            y_width: int,
            min_x_bdy: int,
            max_x_bdy: int,
            optional_delta: int,
            solved: bool = False
    ):
        self.dx_per_iter = dx_per_iter
        self.y_width = y_width
        self.min_x_bdy = min_x_bdy
        self.max_x_bdy = max_x_bdy
        self.optional_delta = optional_delta
        self.solved = solved

def save_json(filename, results_dict):
    n_prev = sum([f.startswith(filename) for f in os.listdir(JSON_DIR)])
    with open(f"{JSON_DIR}/{filename}{n_prev}.json", "w") as outfile:
        json.dump(results_dict, outfile)

