from typing import List, Tuple, Type, TypeVar, Dict

import numpy as np

T = TypeVar("T", bound="Tile")

from src.utils import get_strands, rotate_eyes, _rot0, _rot90, _rot180, _rot270, get_neighbors, unwrap_to_1d

from src.plotter import plot_tile, plot_bezier_style

class Tile:
    """base class for any shape on the grid"""

    def __init__(
            self,
            n_indices: int,
            juxt_index_list: List[Tuple[int, int]],
            double_apart_list: List[Tuple[int, int]],
            name: str,
            thickness: List[float],
            eyes: Dict[int, Tuple[List[Tuple[float]], int]]
    ):
        self.n_indices = n_indices
        self.juxt_index_list = juxt_index_list
        self.double_apart_index_list = double_apart_list
        self.name = name
        self.position_list: List[int] = None  # to be filled after optimizing: corresponds to index_list
        self.juxt_position_list: List[Tuple[int,int]] = None  # after optimizing: corresponds to juxt_index_list
        self.eye_position_list: List[Tuple[int,List[Tuple[float]]]] = None  # after optimizing: contains eye placement data
        self.thickness = thickness
        if eyes is None:
            self.eyes = {}
        else:
            self.eyes = eyes

    @staticmethod
    def grow_body(
            j_start: int,
            j_finish: int,
            _n_indices: int,
            juxt_index_list: List[Tuple[int,int]]
    ) -> (int, List[Tuple[int,int]]):
        """grow the body of the tile by adding connected cells from j_start to j_finish.

        Returns the updated number of indices and the updated juxt_index_list
        """
        started: bool = True
        for j in range(j_start, j_finish):
            _n_indices += 1
            if started:
                started = False
            else:
                juxt_index_list.append((j - 1, j))
        return _n_indices, juxt_index_list

    @classmethod
    def make_man(
        cls: Type[T],
        torso_length: int,
        arm_length: int,
        leg_length: int
            ) -> T:
        """Create a tile with a head, torso, arms, and legs"""

        head_length = 1

        n_indices = 0
        juxt_index_list = []

        thickness = [1.5,0.7] + [0.9]*(torso_length-2) + [1.2]
        for limb_length in [arm_length, leg_length]:
            for _ in [0, 1]:
                thickness.extend([0.7]*(limb_length-1)+[0.8])

        # add the head and torso cells
        n_indices, juxt_index_list = cls.grow_body(
            j_start=0,
            j_finish=head_length+torso_length,
            _n_indices=n_indices,
            juxt_index_list=juxt_index_list
        )

        limb_data = [(arm_length, head_length), (leg_length, torso_length+head_length-1)]

        for limb_length, joint in limb_data:
            for _left_right in [0,1]:
                juxt_index_list.append((joint,n_indices))
                n_indices, juxt_index_list = cls.grow_body(
                    j_start=n_indices,
                    j_finish=n_indices+limb_length,
                    _n_indices=n_indices,
                    juxt_index_list=juxt_index_list
                )

        double_apart_index_list = [(0,2)]

        eyes = {0: ([(-0.2,0.1),(0.2,0.1)], 1)}
        # mouth: , (-0.15,-0.2), (-0.05,-0.25), (0.05,-0.25), (0.15,-0.2)

        return cls(
            n_indices=n_indices,
            juxt_index_list=juxt_index_list,
            double_apart_list=double_apart_index_list,
            name="man",
            thickness=thickness,
            eyes=eyes
        )

    def add_position_list(self, position_list: List[int]) -> None:
        """Take a solution from the optimizer, and update the position list and juxt_position_list accordingly."""
        self.position_list = position_list
        self.juxt_position_list = [
            tuple(sorted([self.position_list[i], self.position_list[j]]))
            for i,j in self.juxt_index_list
        ]
        self.eye_position_list = self.map_eyes(position_list)

    def map_eyes(self, position_list: List[int]) -> List:
        """Take a solution from the optimizer, and process the eye_position_list accordingly."""
        eye_position_list = []
        for i, (eye_coords, pointing_towards) in self.eyes.items():
            eye_w = position_list[i]
            pointing_w = position_list[pointing_towards]
            d = pointing_w - eye_w
            if abs(d) == 1:  # horizontal
                if d == 1:
                    f = _rot270
                else:
                    f = _rot90
            else:  # vertical
                if d > 0:
                    f = _rot180
                else:
                    f = _rot0
            eye_position_list.append((eye_w, [f(*eye) for eye in eye_coords]))

        return eye_position_list

    @classmethod
    def make_worm(cls: Type[T], length: int) -> T:
        """Create a tile that is a worm"""
        n_indices = length
        juxt_index_list = [(i,i+1) for i in range(length-1)]
        double_apart_index_list = []
        thickness = [0.7] + [0.9]*(length-2) + [0.7]

        eyes = {0: ([(-0.2,0.1),(0.2,0.1)], 1)}

        return cls(
            n_indices=n_indices,
            juxt_index_list=juxt_index_list,
            double_apart_list=double_apart_index_list,
            name="worm",
            thickness=thickness,
            eyes=eyes
        )

    def generate_strands_and_plot_bezier(
            self,
            sols: List[np.array],
            color_list: List[int] = None,
            width_factor: float = 3.5,
            smoothing_param: float = 0.9
    ) -> None:
        """Take a list of solutions (x,y representation) and return a dictionary with the strands, widths, and color
        That the plot_tile function can use to plot the solutions. Then call the plot function.

        This will only work if all solutions are of the type (i.e. man / worm) that this instance of the Tile class is.

        The plotter will plot each tile as a collection of strands, think of strands roughly as appendages.
        So each arm, leg, and head will be a separate strand, but the torso will be a strand.
        A strand, is functionally a connected (but not straight) line segment that can be smoothed into a Bezier curve.
        Each strand is smoothed independently of others between its start and end points.
        Bezier curves are not well-defined for tree topologies.

        Args:
            sols: List of numpy arrays, each array is a solution to the set covering problem.
            color_list: List of integers, the color of each solution.
            width_factor: float, the factor to multiply the thickness of each strand by.
            smoothing_param: float, how far from the corners to start the bezier interpolation.
        """
        strands_indices_this_shape = get_strands(self.juxt_index_list)
        strands = []
        widths = []
        colors = []
        for sol, color in zip(sols, color_list):
            for strand_indices in strands_indices_this_shape:
                strands.append([tuple(sol[i]) for i in strand_indices])
                widths.append([self.thickness[i] for i in strand_indices])
                colors.append(color)

        # now create eyes for each man
        eyes = []
        for sol in sols:
            for eye_index, (eye_coords, pointing_towards) in self.eyes.items():
                eyes.extend(rotate_eyes(sol, eye_index, pointing_towards, eye_coords))

        plot_bezier_style(
            strands=strands,
            colors=colors,
            widths=widths,
            eyes=eyes,
            smoothing_param=smoothing_param,
            width_factor=width_factor
        )

    def plot_square_style(
            self,
            sols: List[np.array],
            color_list: List[int] = None,
    ) -> None:
        """Take a list of solutions (x,y representation) from the set covering formulation
        plot using the plot_tile function that represents each cell as a square on a grid.

        This will only work if all solutions are of the type (i.e. man / worm) that this instance of the Tile class is.

        Adds an aesthetic margin on the top and bottom of the grid.
        """
        x_width = max([max(s[:, 0]) for s in sols]) + 1
        max_y = max([max(s[:, 1]) for s in sols])
        y_margin = max_y // 2
        y_height = max_y + y_margin * 2 + 1
        if color_list is None:
            color_list = list(range(len(sols)))
        color_list_expanded = sum([[c] * len(tile_array) for c, tile_array in zip(color_list, sols)], [])
        all_cells: List[int] = []
        all_juxts: List[Tuple[int]] = []
        all_eyes_list: List[Tuple[int, List[Tuple[float]]]] = []
        for h, tile_array in enumerate(sols):
            these_cells = list(map(lambda x: unwrap_to_1d(*x, x_width), tile_array + [0, y_margin]))
            all_cells.extend(these_cells)
            all_juxts.extend([tuple(sorted([these_cells[i], these_cells[j]])) for i, j in self.juxt_index_list])
            all_eyes_list.extend(self.map_eyes(these_cells))
        disconnects = set(get_neighbors(x_width, y_height)) - set(all_juxts)

        plot_tile(
            cell_list=all_cells,
            separate_list=list(disconnects),
            eye_list=all_eyes_list,
            x_width=x_width,
            y_height=y_height,
            color_list=color_list_expanded
        )



