from src.utils import RecurseConfig
from src.combinatoric_optimizer import SolutionCollector
from src.optimizer import TileSolver
from src.classes import Tile
from src.four_color_thm import get_color_indices
from src.utils import read_json_tiles

PLOT_BEZIER = True

MODE = 0

if __name__ == "__main__" and MODE == 0:
    # use a recursive search over a local set covering approach to tesselate.

    # remember that the number of squares in a man tile is 14
    # so having a y_width of 7 or 14 means that every 2 or 1 x-widths, respectively, a solution might be possible.
    recurse_config = RecurseConfig(
        dx_per_iter=6,
        y_width=7,
        min_x_bdy=10,
        max_x_bdy=30,
        optional_delta=3,
    )

    # another configuration that works
    recurse_config2 = RecurseConfig(
        dx_per_iter=5,
        y_width=14,
        min_x_bdy=15,
        max_x_bdy=50,
        optional_delta=3,
    )

    # create the first solution collector
    solution_collector = SolutionCollector(
        other_tiles=[],
        rounding_up=0,
        left_bdy=-1,
        skip_optional_pts=set(),
        recurse_config=recurse_config
    )
    solution_collector.construct_variables()
    if not solution_collector.construct_constraints():
        raise ValueError("No initialization possible with this setup.")
    solution_collector.solve_myself()
    if solution_collector.solved:
        print("Solved and terminated.")
    else:
        print("No solution found.")

if __name__ == "__main__" and MODE == 1:
    # use a permutation of tile placement approach to jigsaw figures into a grid

    n_men = 4

    # collect the tiles we want to fit together into a list of Tile instances.
    tiles = [Tile.make_man(torso_length=3, arm_length=2, leg_length=3) for _ in range(n_men)]

    tile_solver = TileSolver(
        include_tiles=tiles,
        x_width=tiles[0].n_indices  # only works well if all tiles have the same number of indices
    )

    tile_solver.construct_model()
    tile_solver.solve()

if __name__ == "__main__" and MODE == 2:
    # load a solution from a JSON, apply the four color theorem, and plot the solution with Bezier curves.

    sols = read_json_tiles("KH_demo.json")

    coloring = get_color_indices(sols)

    man = Tile.make_man(3, 2, 3)

    if PLOT_BEZIER:
        man.generate_strands_and_plot_bezier(sols, coloring, width_factor=9.5)
    else:
        man.plot_square_style(sols, coloring)






