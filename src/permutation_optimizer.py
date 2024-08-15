# import a cp_sat solver
from typing import List, Tuple

from src.classes import Tile
from src.utils import reduce_tile_list, wrap_to_cell, save_json, get_neighbors

from src.plotter import plot_tile

from ortools.sat.python import cp_model


class TileSolver:
    """A class to create a CP model and solve the tile jigsaw problem."""

    def __init__(self, include_tiles: List[Tile], x_width: int):
        self.include_tiles = include_tiles
        self.x_width = x_width
        self.model = cp_model.CpModel()
        self.pos_vars = []
        self.allowed_diffs = [1, -1, x_width, -x_width]
        self.double_diffs = [2*d for d in self.allowed_diffs]
        self.pos_vars = []
        self.y_min = None
        self.y_max = None
        self.sols_pos: List[List[Tuple[int,int]]] = []

    def construct_model(self, w_max: int = None, w_start: int = 0) -> None:
        """
        solve the tile placement problem
        such that tiles that must be adjacent are adjacent
        and no cells overlap.

        Optionally, one can choose the wrapping number to start from and end at, so multiple solutions can be
        pasted together.

        Args:
            w_max: int, the maximum wrapping number to consider (tiles must be grouped below this number)
            w_start: int, the minimum wrapping number to consider (tiles must be grouped above this number)
        """

        total_n_cells = sum([tile.n_indices for tile in self.include_tiles])

        if w_max is None:
            w_max = total_n_cells + w_start - 1
        elif total_n_cells > w_max - w_start:
            raise ValueError("Tiles do not fit in the grid.")

        self.y_min, self.y_max = w_start // self.x_width, w_max // self.x_width

        allowed_diffs = [1, -1, self.x_width, -self.x_width]
        double_diffs = [2*d for d in allowed_diffs]

        print(
            f"Constructing model with {total_n_cells} cells, "
            f"wrapping between w={w_start} and w={w_max}, x_width={self.x_width}"
        )
        for h, tile in enumerate(self.include_tiles):
            self.pos_vars.append([
                self.model.NewIntVar(w_start, w_max, f"position_{tile.name}_{h}_{i}")
                for i in range(tile.n_indices)
            ])
            for i,j in tile.juxt_index_list:
                for k in range(self.y_min, self.y_max+1):
                    # juxt tiles are not wrapped around the boundary
                    # unfortunately, ORTools does not support modulo INequalities, only equalities
                    self.model.add((self.pos_vars[h][i] + self.pos_vars[h][j]) != 2*k*self.x_width-1)

                # For pairs of cells that are adjacent, the difference in their positions must be in this domain
                juxt_difference = self.model.NewIntVarFromDomain(cp_model.Domain.FromValues(allowed_diffs), 'diff')
                self.model.Add(self.pos_vars[h][i] - self.pos_vars[h][j] == juxt_difference)

            # For pairs of cells that are double apart (face and belly),
            # the difference in their positions must be in this domain
            for i,j in tile.double_apart_index_list:
                juxt_double = self.model.NewIntVarFromDomain(cp_model.Domain.FromValues(double_diffs), 'diff')
                self.model.Add(self.pos_vars[h][i] - self.pos_vars[h][j] == juxt_double)

        # break symmetry: impose an ordering of the tiles to prevent shuffling. This only works when using
        # one kind of tile, otherwise this isn't truly a symmetry of the problem.
        for g in range(len(self.include_tiles)-1):
            self.model.Add(self.pos_vars[g+1][0] < self.pos_vars[g][0])

        all_pos = [pos for tile in self.pos_vars for pos in tile]
        # all distinct constraint: no two cells can be in the same position
        self.model.AddAllDifferent(all_pos)

    def solve(self) -> None:
        print("Solving the permutation model.")
        solver = cp_model.CpSolver()

        # get the solutions
        status = solver.Solve(self.model)
        if status == cp_model.OPTIMAL:
            print("Optimal solution found")
        if status == cp_model.FEASIBLE:
            print("Feasible solution found")
        if status == cp_model.INFEASIBLE:
            print("No solution found")
            return None

        for h, tile in enumerate(self.include_tiles):
            tile.add_position_list([solver.Value(pos) for pos in self.pos_vars[h]])
            self.sols_pos.append(list(map(lambda x: wrap_to_cell(x, self.x_width), tile.position_list)))

        self.output_and_print()

    def output_and_print(self):
        results_dict = {"tiles": self.sols_pos,
                        "juxt_index_list": [tile.juxt_index_list for tile in self.include_tiles]}
        save_json(results_dict=results_dict, filename=f"KH_n{len(self.include_tiles)}_x{self.x_width}.json")
        color_list = list(range(len(self.sols_pos)))
        color_list_expanded = sum([[c] * len(tile_array) for c, tile_array in zip(color_list, self.sols_pos)], [])
        all_cells = sum([tile.position_list for tile in self.include_tiles], [])
        all_juxts = sum([tile.juxt_position_list for tile in self.include_tiles], [])
        all_eyes_list = sum([tile.eye_position_list for tile in self.include_tiles], [])
        disconnects = set(get_neighbors(self.x_width, self.y_max + 1)) - set(all_juxts)
        plot_tile(
            cell_list=all_cells,
            separate_list=list(disconnects),
            eye_list=all_eyes_list,
            x_width=self.x_width,
            y_height=self.y_max + 1,
            color_list=color_list_expanded
        )

    def all_solutions(self) -> List[List[int]]:
        """Return all possible foldings and rotations of a tile, given the constraints of the model.

        We add a constraint that one of the cells (magic elt) of the tile must touch the top boundary, and one must
        touch the left boundary, to uniquely position the tile translationally.
        """

        min_w = self.model.NewIntVar(0, self.x_width-1, "min_w")
        self.model.AddMinEquality(min_w, self.pos_vars[0])
        # stick to top of canvas

        # at least one tile must be 0 modulo x_width
        magic_index = self.model.NewIntVar(0, self.include_tiles[0].n_indices-1, "magic_index")
        # add a magic variable
        domain_used = self.pos_vars[0][0]._IntVar__var.domain  # duplicate the domain of the first variable

        magic_element = self.model.NewIntVar(domain_used[0], domain_used[1], "magic_element")
        self.model.AddElement(magic_index, self.pos_vars[0], magic_element)
        self.model.AddModuloEquality(0, magic_element, self.x_width)

        # Create a solver and solve.
        solver = cp_model.CpSolver()
        solution_collector = VarArraySolutionCollector(self.pos_vars[0])
        solver.SearchForAllSolutions(self.model, solution_collector)

        return solution_collector.solution_list


class VarArraySolutionCollector(cp_model.CpSolverSolutionCallback):

    def __init__(self, variables):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__variables = variables
        self.solution_list = []

    def on_solution_callback(self):
        self.solution_list.append([self.Value(v) for v in self.__variables])


def check_number_unique_foldings(x_width: int):
    """Create a grid, and place man tiles in it (all solutions) to check the number of unique foldings."""
    tiles = [Tile.make_man(torso_length=3, arm_length=2, leg_length=3) for _ in range(1)]

    tile_solver = TileSolver(tiles, x_width)
    tile_solver.construct_model(w_max=x_width * x_width)

    all_solutions = reduce_tile_list(
        tile_solver.all_solutions(),
        torso_length=3,
        arm_length=2,
        leg_length=3,
    )
    print(f"Found {len(all_solutions)} unique foldings of 3-2-3 man tiles")
    # For x_width 8 or 10 (8x8 or 10x10 grid) this gives 6124, which is the same number as we find generatively.
    # For odd grid numbers, for some reason the optimizer produces fewer solutions.
    # This should not depend on grid size, as long as the grid is large fit all solutions, which should from 7x7.

