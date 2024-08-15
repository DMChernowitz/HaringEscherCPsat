import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib as mpl
import os

from typing import List, Tuple, Union

from src.utils import wrap_to_cell, three_pts_collinear
import numpy as np

IMAGE_DIR = "images"


def cart_to_map(
        x: Union[int,float],
        y: Union[int,float],
        square_size: float
) -> Tuple[float,float]:
    return x*square_size, (-y-1)*square_size


def color_spot(n):
    return ((2*n+1)%20)/20


def plot_tile(
        cell_list: List[int],
        separate_list: List[Tuple[int,int]],
        eye_list: List[Tuple[int, List[Tuple[float]]]],
        x_width: int,
        y_height: int,
        color_list: List[int],
        save_dpi: int = 2000,
) -> None:
    """plot a set of tiles as squares, in terms of the wrapping (w) coordinate.

    Args:
        cell_list: list of integers, the cells to plot.
        separate_list: list of tuples of integers, between each of these cells there is a line.
        eye_list: list of where to plot eyes and their orientation.
        x_width: int, the width of the grid.
        y_height: int, the height of the grid.
        color_list: list of integers, the color of each cell.
        save_dpi: int, the dpi of the saved image.
    """

    canvas_x_size = 100
    square_size = canvas_x_size / x_width
    canvas_y_size = y_height * square_size
    fig, ax = plt.subplots()

    rectangles = []
    colors = []

    colormap = mpl.colormaps.get_cmap('tab10')

    # create patch collection
    for cell, color_n in zip(cell_list,color_list):
        x, y = wrap_to_cell(cell, x_width)
        rectangles.append(Rectangle(cart_to_map(x,y,square_size), square_size, square_size))
        colors.append(colormap(color_spot(color_n)))
    collection = PatchCollection(rectangles, facecolors=colors)
    ax.add_collection(collection)

    # now add lines for each juxt
    for edge in separate_list:
        x1, y1 = wrap_to_cell(edge[0], x_width)
        x2, y2 = wrap_to_cell(edge[1], x_width)
        if x1 == x2:
            y_use = max(y1,y2)-1
            left_point = cart_to_map(x1, y_use, square_size)
            right_point = cart_to_map(x1+1, y_use, square_size)
            ax.plot([left_point[0],right_point[0]], [left_point[1],right_point[1]], 'k-')
        elif y1 == y2:
            x_use = max(x1,x2)
            top_point = cart_to_map(x_use, y1-1, square_size)
            bottom_point = cart_to_map(x_use, y2, square_size)
            ax.plot([top_point[0],bottom_point[0]], [top_point[1],bottom_point[1]], 'k-')

    scatter_x, scatter_y = [], []
    # plot eyes
    for eye_w, eye_coords in eye_list:
        base_x, base_y = wrap_to_cell(eye_w, x_width)
        for eye in eye_coords:
            x,y = cart_to_map(base_x+eye[0]+0.5, base_y+eye[1]-0.5, square_size)
            scatter_x.append(x)
            scatter_y.append(y)

    ax.scatter(scatter_x, scatter_y, s=20, c='k', marker='o')

    ax.set_aspect('equal')
    # plot outline
    ax.plot([0, canvas_x_size, canvas_x_size, 0, 0], [0, 0, -canvas_y_size, -canvas_y_size, 0], 'k-')
    # remove axis labels
    ax.set_xticks([])
    ax.set_yticks([])

    filename = f"KH_{len(cell_list)}_"
    n_files_before = sum([f.startswith(filename) for f in os.listdir(IMAGE_DIR)])
    # save the figure
    plt.savefig(f"{IMAGE_DIR}/{filename}{n_files_before}.png", dpi=save_dpi, bbox_inches='tight')

    plt.show()


def plot_bezier_style(
        strands: List[List[tuple[int,int]]],
        colors: List[int],
        widths: List[List[float]],
        eyes: List[Tuple[float, float]],
        smoothing_param: float = 0.5,
        width_factor: float = 10,
        save_dpi: int = 2000
) -> None:
    """plot the topologically separate strands in a bezier style.

    Args:
        strands: list of lists of tuples, each tuple is a point (x,y). Each strand will be connected
        start to finish by a series of bezier curves.
        colors: list of integers, the color of each strand.
        eyes: List of coordinates to plot eyes
        widths: list of lists of floats, the width of each strand at each point.
        smoothing_param: float, the smoothing parameter of the bezier curves. If 0: get straight lines and right angles.
            if ~1: smooth curves touch when going around corners.
        width_factor: float, the factor by which the width of the strands is multiplied. This needs to be chosen
            based on the total plot canvas size: for larger total x and y spans, the width_factor should be smaller.
        save_dpi: int, the dpi of the saved image.
    """
    if not 0 <= smoothing_param < 1:
        raise ValueError("smoothing_param must be in [0,1)")

    plot_pts_density = 40  # the number of points to plot per unit in the grid. Higher is smoother.

    prop = 0.5*smoothing_param
    oprop = prop/(1-prop)
    bezier_n_pts = int(1.62*prop*plot_pts_density)

    bg_width_thickness = 0.4 # the ratio of the background width to the width of the strands

    # 1.62 is the ratio of the arclength of bezier curve to one of its horizontal / vertical sides
    t = np.linspace(0, 1, bezier_n_pts)
    t0, t1, t2 = (1-t)**2, 2*(1-t)*t, t**2

    N_strands = len(strands)
    if not len(colors) == N_strands:
        raise ValueError("colors must be the same length as strands")
    if not len(widths) == N_strands:
        raise ValueError("widths must be the same length as strands")
    for n in range(N_strands):
        if not len(strands[n]) == len(widths[n]):
            raise ValueError(f"widths[{n}] must be the same length as strands[{n}]")

    xmin, ymin, xmax, ymax = (
        a+f([f([xy[j] for xy in strand]) for strand in strands])
        for j, f, a in zip([0, 1, 0, 1], [min, min, max, max], [-1, -1, 1, 1])
    )
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    # plot outline
    ax.plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], 'k-')
    #remove black box
    ax.set_frame_on(False)
    # remove axis labels
    ax.set_xticks([])
    ax.set_yticks([])
    x_pts = []
    y_pts = []
    w_pts = []
    bg_w_pts = []
    color_pts = []

    colormap = mpl.colormaps.get_cmap('tab10')

    def _extend_lists(xs, ys, wid0, wid1, col):
        x_pts.extend(xs)
        y_pts.extend(ys)
        w_now = (width_factor*np.linspace(wid0, wid1, len(xs)))**2
        bg_w_now = (width_factor*(bg_width_thickness+np.linspace(wid0, wid1, len(xs))))**2
        w_pts.extend(w_now)
        bg_w_pts.extend(bg_w_now)
        color_pts.extend([colormap(color_spot(col))]*len(xs))

    def _add_straight_line(pt0, pt1, wid0, wid1, col):
        # only for vertical or horizontal lines
        line_len = abs(pt0[0] - pt1[0]) + abs(pt0[1] - pt1[1])
        n_pts = int(line_len*plot_pts_density)
        xs = np.linspace(pt0[0], pt1[0], n_pts)
        ys = np.linspace(pt0[1], pt1[1], n_pts)
        _extend_lists(xs, ys, wid0, wid1, col)

    def _add_quadratic_bezier(pt0,pt1,pt2,wid0,wid2,col):
        # add a quadratic bezier curve with control point pt1, right angles
        xs = t0*pt0[0]+t1*pt1[0]+t2*pt2[0]
        ys = t0*pt0[1]+t1*pt1[1]+t2*pt2[1]
        _extend_lists(xs, ys, wid0, wid2, col)

    def _process_strand(_strand,_width,_color):
        # recursive function to plot a strand
        if len(_strand) == 1:
            return  # terminate recursion

        elif len(_strand) == 2 or three_pts_collinear(*_strand[:3]):
            # end of strand, or strand goes straight
            midway_pt_0 = _strand[1]
            midway_width_0 = _width[1]
            remaining_strand = _strand[1:]
            remaining_width = _width[1:]

        else:  # go around corner
            dxdy = [_strand[1][j]-_strand[0][j] for j in [0,1]]
            d = abs(sum(dxdy))
            # d is 1, or 1-prop

            prop1 = oprop if d < 1 else prop

            midway_pt_0 = tuple(_strand[1][j]-prop*round(dxdy[j]) for j in [0, 1])
            midway_width_0 = _width[0]*prop1+_width[1]*(1-prop1)

            # after corner
            midway_pt_1 = tuple(prop*_strand[2][j]+(1-prop)*_strand[1][j] for j in [0, 1])
            midway_width_1 = _width[2]*prop+_width[1]*(1-prop)

            _add_quadratic_bezier(midway_pt_0,_strand[1],midway_pt_1,midway_width_0,midway_width_1,_color)

            remaining_strand = [midway_pt_1]+_strand[2:]
            remaining_width = [midway_width_1]+_width[2:]

        _add_straight_line(_strand[0], midway_pt_0, _width[0], midway_width_0, _color)
        _process_strand(remaining_strand, remaining_width, _color)

    for strand, color, width in zip(strands, colors, widths):
        # the strand is assumed to be connected and neighboring points.
        _process_strand(strand, width, color)

    black_pts = [(0,0,0)]*len(x_pts)
    for _w,_c in [(bg_w_pts,black_pts),(w_pts,color_pts)]:
        ax.scatter(x_pts, y_pts, s=_w, c=_c)
        # scatter the points

    # plot eyes
    x_eye,y_eye = zip(*eyes)
    ax.scatter(x_eye, y_eye, s=(width_factor*bg_width_thickness/3)**2, c=[(0,0,0)]*len(eyes), marker='s')

    filename = f"KH_{len(strands)}_bezier_"
    n_files_before = sum([f.startswith(filename) for f in os.listdir(IMAGE_DIR)])
    # save the figure
    plt.savefig(f"{IMAGE_DIR}/{filename}{n_files_before}.png", dpi=save_dpi, bbox_inches='tight')

    plt.show()






