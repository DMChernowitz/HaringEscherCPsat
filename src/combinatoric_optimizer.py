import numpy as np
import random

from typing import List, Tuple, Dict, Union, Set, TypeVar

from src.generative import GrouperOfMen, Outline
from src.classes import Tile

from src.utils import show_chunk, RecurseConfig, save_json

from ortools.sat.python import cp_model


SC = TypeVar("SC", bound="SolutionCollector")


class SolutionCollector(cp_model.CpSolverSolutionCallback):
    def __init__(self,
                 other_tiles: List[np.array],
                 rounding_up: int,
                 left_bdy: int,
                 skip_optional_pts: Set[Tuple[int, int]],
                 recurse_config: RecurseConfig,
                 grouper_of_men: GrouperOfMen = None,
                 target_x_bdy: int = None,
                 progress_tree_string: List[str] = None,
                 chunks: List[Set[Tuple[int, int]]] = None,
                 ) -> None:
        """Create a new SolutionCollector object that also constructs the local model.

        other_tiles: list of np.arrays of shape [(x_1,y_1), ... (x_n,y_n)] from previous iterations
        rounding_up: int, if non-zero, the buffer on the right until the straight edge attempt
        left_bdy: int, the left boundary of the grid
        grouper_of_men: GrouperOfMen object used to produce all the tile set coverings.
        skip_optional_pts: set of optional points already covered in previous iterations
        target_x_bdy: int, the right boundary of the grid, above which the solver first attempts a straight edge
            This value is increased by self.rc.dx_per_iter whenever the straight edge attempt fails, until self.rc.max_x_bdy.
        recurse_config: RecurseConfig object with the parameters for the recursion, and the solved boolean.
        progress_tree_string: list of strings that represent the progress of the search tree
        chunks: list of sets of points that are already covered in previous iterations
        """
        super().__init__()
        # grid parameters

        # global progress parameters
        self.other_tiles = other_tiles
        self.rounding_up = rounding_up
        self.left_bdy = max(0,left_bdy)
        self.grouper_of_men = grouper_of_men
        self.rc = recurse_config
        if target_x_bdy is None:
            target_x_bdy = self.rc.min_x_bdy
        self.target_x_bdy = target_x_bdy
        if grouper_of_men is None:
            grouper_of_men = GrouperOfMen()
        self.grouper_of_men = grouper_of_men
        if progress_tree_string is None:
            progress_tree_string = ["Search Tree: ----"]
        self.progress_tree_string = progress_tree_string
        if chunks is None:
            chunks = []
        self.chunks = chunks

        if rounding_up:
            self.right_bdy = left_bdy + rounding_up
            self.optional_pts = set()
        else:
            self.right_bdy = left_bdy + self.rc.dx_per_iter
            self.optional_pts = {
                (x, y) for x in range(self.right_bdy, self.right_bdy + self.rc.optional_delta)
                for y in range(self.rc.y_width)
            }
        pre_pts_in_grid = {(x, y) for x in range(self.left_bdy, self.right_bdy) for y in range(self.rc.y_width)}
        self.all_p, self.grouped_p = self.grouper_of_men.get_all_men_plus_translations(
            x_left=self.left_bdy,
            x_right=self.right_bdy + self.rc.optional_delta,
            y_bottom=0,
            y_top=self.rc.y_width
        )
        self.pts_in_grid = pre_pts_in_grid - skip_optional_pts
        show_chunk(
            y_width=self.rc.y_width,
            x_width=self.rc.max_x_bdy,
            pts=self.pts_in_grid,
            optional_pts=self.optional_pts,
            chunk_history=self.chunks
        )

        self.solution_list = []

        self.model = cp_model.CpModel()

        # maps a grid coordinate to the tile shape, dx and dy index that cover it
        self.set_cover_hash_map: Dict[Tuple[int, int], List[cp_model.IntVar]] = {}

        self.all_pts = self.pts_in_grid.union(self.optional_pts)
        self.n_men = (len(self.pts_in_grid) + len(self.optional_pts) // 2) // self.grouper_of_men.n_squares_per_tile
        self.all_bool_vars: Dict[Outline, List[List[Union[cp_model.IntVar, None]]]] = {k: [] for k in self.all_p.keys()}
        self.all_vars_list = []
        self.optional_cover_bools: Dict[Tuple[int, int], cp_model.IntVar] = {}
        self.print_tree_prog()

    @property
    def solved(self):
        return self.rc.solved

    @solved.setter
    def solved(self, value):
        self.rc.solved = value

    def construct_variables(self):
        """Creates booleans per tile (including position of the tile), whether to include it in the set covering.

        Variables corresponding to the optional points are created separately.
        Variables for tiles that don't fit in the grid are discarded.
        """
        print("Creating variables. ", end="")
        total_variables = 0
        discarded = 0
        for h, (outline, base_shape) in enumerate(self.all_p.items()):
            xw, yw, _, _ = base_shape.shape
            for dx in range(xw):
                self.all_bool_vars[outline].append([])
                for dy in range(yw):
                    if self.inside_grid(base_shape[dx, dy]):
                        self.all_bool_vars[outline][dx].append(self.model.NewBoolVar(f"shape_{h}_{dx}_{dy}"))
                        self.all_vars_list.append(self.all_bool_vars[outline][dx][dy])
                        for tile_nr in base_shape[dx, dy]:
                            x, y = tile_nr
                            self.set_cover_hash_map.setdefault((x, y), []).append(self.all_bool_vars[outline][dx][dy])
                        total_variables += 1
                    else:
                        self.all_bool_vars[outline][dx].append(None)
                        discarded += 1

        print(f"Number tiles: {self.n_men}. Boolean variables: {total_variables}. Discarded variables: {discarded}")

    def construct_constraints(self) -> bool:
        """
        Create the constraints for the set covering problem. This includes
        - each required grid point must be covered by exactly one tile
        - each optional grid point must be covered by at most one tile
        - if an optional point is covered, the point to its left must also be covered
        - no long peninsulas at the boundaries, because empirically the next iteration no covering is possible
        - the total number of tiles present should be the predicted amount based on a heuristic of the size of the grid.

        :return: Boolean whether the constraints could be created consistently
        """
        print("Creating constraints ", end="")
        # create the constraints
        for x, y in self.pts_in_grid:
            if (x, y) not in self.set_cover_hash_map:
                print(f"No covering possible with this edge. Missing tile at required pt {x},{y}")
                return False
            self.model.Add(sum(self.set_cover_hash_map[(x, y)]) == 1)
        for x, y in self.optional_pts:
            if (x, y) not in self.set_cover_hash_map:
                print(f"No covering possible with this edge. Missing tile at optional pt {x},{y}")
                return False
            # add an aux variable whether this optional point is used
            self.optional_cover_bools[(x, y)] = self.model.NewBoolVar(f"optional_{x}_{y}")
            self.model.Add(sum(self.set_cover_hash_map[(x, y)]) == self.optional_cover_bools[(x, y)])

        for x, y in self.optional_pts:
            if (x + 1, y) in self.set_cover_hash_map:
                # a variable can't be zero to its right is one
                self.model.Add(self.optional_cover_bools[(x, y)] >= self.optional_cover_bools[(x + 1, y)])

        if len(self.optional_pts) > 0:
            minx_opt, maxx_opt = [f(x for x, y in self.optional_pts) for f in [min, max]]
            # no long peninsulas at the boundaries
            for f, delt in [(min, 1), (max, -1)]:
                y_opt = f([y for x, y in self.optional_pts])
                self.model.Add(self.optional_cover_bools[(minx_opt, y_opt)] >= self.optional_cover_bools[(maxx_opt, y_opt + delt)])

        # create additional constraint
        self.model.Add(sum(self.all_vars_list) == self.n_men)
        print(f"for {len(self.all_pts)} grid points.")

        return True

    def solve_myself(self) -> bool:
        """Solve the set covering problem and collect the solutions.

        Will iterate over solutions and terminate when the end of a branch is successful at creating a straight edge.

        :return: Boolean whether the branch was successful
        """
        print("Solving...")
        # create the solver
        solver = cp_model.CpSolver()

        # tell the solver to get all solutions
        solver.parameters.enumerate_all_solutions = True

        # pass myself as the collector for callback, this class has the on_solution_callback method
        status = solver.Solve(model=self.model, solution_callback=self)

        prestr = f"At level {len(self.progress_tree_string)-1}: Solver time {str(solver.WallTime())[:6]} s, "

        print_dict = {
            cp_model.INFEASIBLE: "no solution found.",
            cp_model.FEASIBLE: f"{len(self.solution_list)} feasible solutions found.",
            cp_model.OPTIMAL: f"exhausted all {len(self.solution_list)} solution(s)."
        }
        print(f"{prestr}{print_dict.get(status, 'unknown status.')}")
        return bool(self.solved)

    def inside_grid(self,_tile):
        return {(int(_x), int(_y)) for _x, _y in _tile}.issubset(self.all_pts)

    def on_solution_callback(self):
        """Callback function for the solver, called when a solution is found. The recursion must take place here."""
        included_tiles = []

        for outline, base_shape in self.all_p.items():
            xw, yw, _, _ = base_shape.shape
            for dx in range(xw):
                for dy in range(yw):
                    this_var = self.all_bool_vars[outline][dx][dy]
                    if (this_var is not None) and self.Value(this_var):
                        n_options = len(self.grouped_p[outline])
                        minx, miny = min(base_shape[dx, dy, :, 0]), min(base_shape[dx, dy, :, 1])
                        # translate to correct place in grid
                        use_base_shape = self.grouped_p[outline][random.randint(0, n_options - 1)] + [minx, miny]
                        included_tiles.append(use_base_shape)

        used_optional_pts = set()

        for x, y in self.optional_pts:
            if self.Value(self.optional_cover_bools[(x, y)]):
                used_optional_pts.add((x, y))

        self.solution_list.append(
            {"included_tiles": included_tiles,
             "used_optional_pts": used_optional_pts
             }
        )
        self.print_tree_prog()
        self.continue_build()

        if self.solved:
            self.stop_search()

    def print_tree_prog(self):
        ns = str(len(self.solution_list))
        tree_str_list = self.progress_tree_string + [
            ("|"+ns+"|"*self.rounding_up if
             self.rounding_up else
             ns+ "----")
        ]
        print("".join(tree_str_list))

    def continue_build(self) -> bool:
        """Continue building the grid, or terminate the search.

        Returns a boolean whether this branch was successful.
        """

        if self.rounding_up:
            # there was a solution while in a rounding up phase; terminate without creating a deeper level
            self.output_solution()
            self.solved = True
            return True  # terminated successfully, don't create a deeper level
        elif self.right_bdy >= self.target_x_bdy:
            # start looking for a solution with a straight edge
            for rounding_up in range(self.rc.optional_delta, self.rc.dx_per_iter + 1):
                if self.create_deeper_level(rounding_up):
                    return True
            # no solution found with a straight edge: increase the allowed right boundary
            return self.create_deeper_level(0, increase_bdy=True)
        else:
            # continue building the grid
            return self.create_deeper_level(0)

    def create_deeper_level(self, rounding_up: int, increase_bdy: bool = False) -> bool:
        """Create a deeper level of the search tree."""

        new_target_x_bdy = self.target_x_bdy + (self.rc.dx_per_iter if increase_bdy else 0)
        if new_target_x_bdy > self.rc.max_x_bdy:
            # we have reached the maximum boundary: backtrack
            return False

        progress_tree_string = self.progress_tree_string + [f"({len(self.solution_list)})----"]
        chunks = self.chunks + [self.pts_in_grid.union(self.solution_list[-1]["used_optional_pts"])]

        deeper_solution_collector = SolutionCollector(
            other_tiles=self.other_tiles+self.solution_list[-1]["included_tiles"],
            rounding_up=rounding_up,
            left_bdy=self.right_bdy,
            grouper_of_men=self.grouper_of_men,
            skip_optional_pts=self.solution_list[-1]["used_optional_pts"],
            target_x_bdy=new_target_x_bdy,
            recurse_config=self.rc,
            progress_tree_string=progress_tree_string,
            chunks=chunks
        )

        deeper_solution_collector.construct_variables()
        if not deeper_solution_collector.construct_constraints():
            return False
        return deeper_solution_collector.solve_myself()  # success, but not terminated


    def output_solution(self):
        sols = self.other_tiles+self.solution_list[-1]["included_tiles"]

        man = Tile.make_man(torso_length=3, arm_length=2, leg_length=3)
        man.plot_square_style(sols)

        # save the results to a json
        results_dict = {
            "tiles": [sol.tolist() for sol in sols],
            "juxt_index_list": [man.juxt_index_list]*len(sols)
        }
        filename = f"KH_y{self.rc.y_width}_c{self.rc.min_x_bdy}_t{self.rc.dx_per_iter}"
        save_json(filename, results_dict)








