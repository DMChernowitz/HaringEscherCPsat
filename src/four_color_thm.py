import numpy as np
from typing import List, Tuple, Dict

from ortools.sat.python import cp_model


def get_color_indices(sols: List[np.array]) -> List[int]:
    """Take a list of solutions that are a tiling, and return a list of numbers so no neighbors have the same number

    sols: list of np.arrays of shape [(x_1,y_1), ... (x_n,y_n)]

    We use the freedom of the objective function to return the solution that has the most equal number of each color.
    """
    l = len(sols)
    print(f"Attempting to color {l} tiles")
    x_width, y_width = (max([max(s[:,j]) for s in sols])+1 for j in [0,1])

    inverse_tile_hashmap = get_inverse_tile_hashmap(sols)

    # now get neighbor matrix
    neighbor_matrix = get_neighbor_matrix(inverse_tile_hashmap, l, x_width, y_width)

    for n_colors in [3,4,5]:
        model = cp_model.CpModel()
        color_assignment: List[List[cp_model.IntVar]] = []  # create variables
        for tile_nr in range(l):
            # we want a boolean per tile-color pair. It is unit if that tile has that color.
            color_assignment.append([model.NewBoolVar(f"color_{tile_nr}_{i}") for i in range(n_colors)])
            # a tile has exactly one color
            model.add(sum(color_assignment[tile_nr]) == 1)

        for i in range(l):
            for j in range(l):
                if neighbor_matrix[i][j] and i != j:
                    # by our method, a tile will neighbor itself. Still it should not have a different color to itself.
                    for c in range(n_colors):
                        model.add(color_assignment[i][c] + color_assignment[j][c] <= 1)  # neighbors must be distinct

        # total number of tiles per color
        # The boolean representation (not a color IntVar per tile) makes counting the number of tiles per color natural.
        total_colors = [sum(color_assignment[i][c] for i in range(l)) for c in range(n_colors)]
        minimum_color_count = model.NewIntVar(1, l, "minimum_color_count")
        model.add_min_equality(minimum_color_count, total_colors)
        # make the count of the rarest color as large as possible: most equal distribution
        model.maximize(minimum_color_count)

        solver = cp_model.CpSolver()
        status = solver.Solve(model)
        print(f"Status: {status}")
        if status == cp_model.OPTIMAL:
            print(f"Found a {n_colors} coloring")
            solution_numbers: List[int] = []
            for tile_nr in range(l):
                for color_nr in range(n_colors):
                    if solver.Value(color_assignment[tile_nr][color_nr]):
                        solution_numbers.append(color_nr)
                        break
            return solution_numbers
        else:
            print(f"No {n_colors} coloring found")
    return list(range(l))


def get_inverse_tile_hashmap(sols: List[np.array]) -> Dict[Tuple[int, int], int]:
    """Take a list of solutions (order is important) in terms of their (x,y) coordinates,
     And return a dict (hashmap) that, indexed per x,y-gridpoint, which tile is present there.
     """
    inverse_tile_hashmap: Dict[Tuple[int, int], int] = {}
    for h, tile_array in enumerate(sols):
        for x, y in tile_array:
            if (x, y) in inverse_tile_hashmap:
                print(f"Error: tile {h} overlaps with tile {inverse_tile_hashmap[(x, y)]}")
            inverse_tile_hashmap[(x, y)] = h
    return inverse_tile_hashmap


def get_neighbor_matrix(inverse_tile_hashmap, l, x_width, y_width) -> List[List[int]]:
    """Go cell by cell, and check its right and bottom neighbor (if present).
    Add the two tiles to the neighbor juxtaposition matrix.
    """
    neighbor_matrix = [[0]*l for _ in range(l)]
    for x in range(x_width):
        for y in range(y_width):
            if not (x,y) in inverse_tile_hashmap:
                continue
            if x + 1 < x_width and (x+1, y) in inverse_tile_hashmap:
                n1, n2 = inverse_tile_hashmap[(x, y)], inverse_tile_hashmap[(x + 1, y)]
                neighbor_matrix[n1][n2] = 1
            if y + 1 < y_width and (x, y+1) in inverse_tile_hashmap:
                n3, n4 = inverse_tile_hashmap[(x, y)], inverse_tile_hashmap[(x, y + 1)]
                neighbor_matrix[n3][n4] = 1
    return neighbor_matrix

