# HaringEscherCPsat

This is a hobby project that uses Constraint Programming and some recursion to generate artistic images. They are inspired by the work of Keith Haring (for their aesthetics) and MC Escher (for their pattern structure). The goal is to create large interlocking patterns of organic shapes. But unlike in Escher's work, and more like in Haring's work, we don't want translational symmetry (or repetition).

# Getting started

First, clone this repo or download a zip file and extract it, then open the project with your favorite IDE.
This project only uses three standalone packages, and their own prerequisites. These are `numpy`, `matplotlib`, and Google's `ortools`. I suspect Numpy is in fact a prerequisite of OR-Tools, but I'm mentioning it separately just in case.
Install them from the terminal with the command `pip install requirements.txt`, or navigate to them by hand if your IDE has integrated package management.

Once these packages are working, an example script utilizing the main functionalities can be found in `main.py`. Set the `MODE` parameter to 0,1 or 2, and run the script.
Take note that the script may take some minutes to terminate, if mode 0 or 1 are chosen for the demo.

### What this repo does

This repo contains scripts to generate tesselated figures without repetition, and to plot them. I will call this the Tiling Problem. See an example below.
In the example, faces were added later for clarity (and shits and giggles.)

![An example of a solution to the tiling problem](./readme_imgs/KH_crop.png)

The logic is contained in the `src` directory.

The starting point should be the `permutation_optimizer.py` and `classes.py`. 
These contain a set of classes and functions to define figures in terms of cells on a grid, and juxtaposition constraints on those cells.
Notably, there are no explicit constraints on the actual placement of each cell, though by construction it must lie inside a certain rectangular area.

These scripts can be used to create interlocking figures of little men with variable arm, torso, and leg length, or worms of varying length. 
They can also be modified to allow for figures of a different topology than lines or bipeds. However, figures of over 100 cells will not easily converge with a reasonable time (some hours).

The more advanced, but also more constrained in functionality, is `set_covering_optimizer`. it generates a ribbon shaped image of interlocking man-tiles. 
This contains a class that interacts with the google OR-tools CPsat solver in a recursive manner to perform a tree search on a subset of the solution space, only solving the tiling in one local chunk at a time.
This way, in principle the runtime should be approximately linear in the size of the area covered, assuming one is not too strict on the exact size of the image at the point of termination.

Finally, there is a script that applies a coloring to a solution to the tiling problem, such that no neighboring tiles have the same color. 
This is always possible with at most 4 colors, due to the four-color theorom. Hence the name `four_color_thm.py`.

Now let's look at the math / CS behind these two ways to describe the tiling problem, the _permutation_ and the _set covering_ formulations.

# Permutation Optimizer

The mathematical formulation of this approach is the following: 


Consider a finite, square grid of cells (pixels), on which we wish to tesselate figures. The width of the grid is `X`. The figures consist of neighboring cells. However, the pattern need not have repetition or any kind of translational symmetry.

A proposed representation of the grid is to have a single integer coordinate `w`, starting at 0, and wrapping around the grid in horizontal shelves or rows. For the case of width `X=7`, it would look like:

| 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| 7 | 8 | 9 | 10 | 11 | 12 | 13 |
| 14 | 15 | 16 | 17 | 18 | 19 | 20 |
| 21 | 22 | 23 | 24 | 25 | 26 | 27 |
| 28 | 29 | 30 | 31 | 32 | 33 | 34 |


Then the mapping to 2d is the obvious

    x(w)  =  w % X
    y(w)  =  floor(w/X)

Let us define a _tile_ as a tuple of an integer `n`, the number of cells in the tile, and a set of edges `E` of neighboring cells indices.

    n  =  14 
    E  =  (0,1), (1,2), (2,3), (3,4), (1,4), (4,5), (1,6), (6,7), (3,8), (8,9), (9,10), (3,11), (11, 12), (12,13)

Let us use this tile above as the motivating example: a man-shaped tile (man-tile) with torso length of 3 cells, arm length of 2 cells, and leg length of 3 cells. In one specific folding, it looks like the image below.
Remember that the numbers on this figure are not the `w` coordinate, but the index counting along the cells of this tile. So the head of the man has w-coordinate `w[0]`.

![The shape of man tile](./readme_imgs/mantile.png)

The constraint programming formulation is when we take the intuitive demands, valid independently for each tile.

    |w[i] - w[j]|  =  1 or X            for all (i,j) in E
    (w[i] + w[j] + 1) % X  !=  0        for all (i,j) in E

The first demand is simply the juxtaposition. If `w[i]-w[j]` is +1 or -1, the two cells are left and right. If it is +X or -X, they are vertically stacked.
The only caveat is that we don't want them to be in subsequent grid points, if those points are exactly on either side of the boundary (i.e. 6 and 7 in the example). That is the function of the second constraint.

Furthermore, we must create a set of `w[i]` integer variables for each tile, and all of the `w`'s must be distinct, so no tiles overlap.

    w[a][i]  !=  w[b][j]                for all tiles a,b, for all indices i,j

Google OR-Tools has an `AllDifferent` demand that the solver handles efficiently, immediately reducing the state space as to permutation problem. Hence the name of the module. If we create the variables to have minimum 0, and maximum the last square on the grid `w_max`,
and furthermore, we have `m` tiles, then `m*n=w_max+1` uses every cell on the grid exactly once. It turns out, that still not all grids that satisfy this have a solution. For instance, a 14x7 grid with 7 man-tiles has no solution, but 14x6 (6 tiles) and 14x8 (8 tiles) do.

As a small point of aesthetics, I have also added a _double juxt_ constraint between indices 0 and 2. This forces the arms to be on either side of the torso, with the head between them, and the hips as well, making it much easier to recognize as a human.

This is a constraint problem and there is no objective function. A satisfying solution is all that is required. In this representation, Google OR-tools seems to cap out (on my local machine) at a 14x8 grid: then there are 112 variables, and the state space cardinality is upper bounded by 112!, which is of course quite large. 

### Scripts

The central objects are the `TileSolver` from `permutation_optimizer.py` and `Tile` class from the `classes.py` file.

The `Tile` class is intended as a kind of superclass for all shapes of tiles that can be tesselated. However, instances are not obtained by inheritance,
but by specific class methods `make_man` and `make_worm`, that return an instance of the class with the correct parameters. The first creates a man-tile.

- `Tile.make_man()`
  - arguments:
    - `torso_length`: int, the number of cells connected to the head that form the torso.
    - `arm_length`: int, the number of cells per arm
    - `leg_length`: int, the number of cells per leg.
  - returns: An instance of the `Tile` class with properties:
    - `n_indices`: int, equalling `torso_length` + 2 * `arm_length` + 2 * `leg_length` + 1
    - `juxt_index_list`: list, with the appropriate edges connecting the body parts 
    - `double_apart_list`: list, equals [(0,2)], forcing the head opposite the hips.
    - `thickness`: list, of thickness (float) per cell for rendering the body parts organically by plotters.
    - `eyes`: list, containing information of where to plot markings that serve as eyes, on index 0, opposite index 1.
    
`Tile.make_worm()` is a similar function returning a tile without splitting in its topology.

- `Tile.make_worm()`
  - arguments:
    - `length`: int, the number of cells connected to the head that form the entire worm.
  - returns: An instance of the `Tile` class with properties:
    - `n_indices`: int, equalling `length`
    - `juxt_index_list`: list, with the appropriate edges connecting each index to the next.
    - `double_apart_list`: list, empty: no constraints on the shape.
    - `thickness`: list, of thickness (float) per cell
    - `eyes`: list, containing information of where to plot markings that serve as eyes, on index 0, opposite index 1.

The way to utilize this class is to make a list of `Tile` instances, and feed that list into the `TileSolver` class. This class is intended to solve the tiling problem.

- `TileSolver.__init__()`
  - arguments:
    - `include_tiles`: list, containing `Tile` instances you want to jigsaw into the grid.
    - `X_width`: int, the wrapping parameter for the grid in the 1-D representation
  - returns: An `TileSolver` instance initialized with a constraint solver. Notable properties:
    - `include_tiles`: list, containing the tile instances.
    - `model`: cp_model.CpModel()
    - `sols_pos`: list, now empty, but containing the tile positions after solving.

Then we must call

- `TileSolver.construct_model()`
  - arguments:
    - `w_start`: int, optional, the minimal `w` coordinate. Default 0.
    - `w_max`: int, optionalm the maximal `w` coordinate. Default exactly fitting the tiles used in initialization.
  - returns: 
    - Nothing, but constructs all variables for each cell of each tile, and adds the constraints as defined above to the model. 
  - notes:
    - It is possible to solve on a larger grid than necessary, then some cells will be empty. 
    - We can also jog the grid using the optional params, so the solutions start later. This way, one solution can be in domain `[0,w*)`, and the next in `[w*,w**)`. With some extra work, these can be pasted together, reducing the global complexity.

And finally, call

- `TileSolver.solve()`
  - arguments:
    - None
  - returns: Nothing, but solves the model, and if feasible, 
    - plots the result using the `plot_tile` function from the `plotter.py` file.
    - updates with the true `w` coordinates the `position_list`, `juxt_position_list` of each `Tile` instance in `include_tiles`. 
    - updates the `sols_pos` list with the (x,y) coordinates of each cell of each tile.
    - Outputs the solution to a json file and saves it in the `jsons` directory.

# Set Covering Optimizer

Intuitively, another way to approach this problem is like a jigsaw puzzle. Instead of folding a single tile in all allowed ways, let us find all possible foldings and rotations, and then attempt to fit them together to exactly cover the grid.
This moves complexity away from the solving, and into the formulation of the constraint program problem, and is a typical reframing. Now the question becomes a set-covering problem.

Let us consider the points on the grid. They form a set `S`.


|     | x=0 | x=1 | x=2 | x=3 | x=4 | x=5 |
|---|----|-----|---|---|---|---|
| y=0 | 0,0 | 1,0 | 2,0 | 3,0 | 4,0 | 5,0 |
| y=1 | 0,1 | 1,1 | 2,1 | 3,1 | 4,1 | 5,1 |
| y=2 | 0,2 | 1,2 | 2,2 | 3,2 | 4,2 | 5,2 |
| y=3 | 0,3 | 1,3 | 2,3 | 3,3 | 4,3 | 5,3 |

And additionally, we have at our disposal, a set of _subsets_ of `S` of points,

    S_0  =  {(x_0_0, y_0_0), (x_0_1, y_0_1), (x_0_2, y_0_2), ... }
    S_1  =  {(x_1_0, y_1_0), (x_1_1, y_1_1), (x_1_2, y_1_2), ... }
    ...  =  ...

etc. If these subsets `S_i` are well-chosen, then it will be possible to each grid point in `S` covered by 
exactly _one_ of the `S_i`'s. So let the `S_i` be the collection of points that would be covered by any possible way to fold, rotate, and translate the tiles.

Then the model we need to solve is the following: choose a Boolean variable for each set. `v_j`. If `v_j==1`, this means including the subset in the set cover. Otherwise, exclude it.
Now let `I_xy` be the collection of all indices `i` such that point `x,y` is in `S_i`.

    $\sum_{i\in\I_x,y} v_i == 1, \forall (x,y) \in S$