import numpy as np
from typing import List, Tuple, Dict


class Outline:
    """The outline of a particular instantiation of a tile is determined by the squares that are part of the tile.

    Thus two tiles with a different folding, but occupying the same squares, will have the same outline.

    In the set covering problem, we only care about the outline of the tiles, not the specific folding.

    The outlines will be used to index dictionaries, and create sets of unique outlines,
    so we need to define equality and hashing dunder methods.
    """

    def __init__(self, pts_array: np.array):
        """
        In a tile, the order of the squares determines the function of each square, and the juxtaposition (connections)

        For the outline, we only care about the squares that are part of the tile, not the order. So sorting the points
        will give us a unique identifier for the outline.

        Create the outline by just sorting the points in the array: then we don't care in what order they started."""
        self.pts = tuple(sorted([(int(a),int(b)) for a, b in pts_array]))

    def __hash__(self):
        """Tuples of tuples are hashable, and these are the unique identifiers for the outlines."""
        return hash(self.pts)

    def __repr__(self):
        return f"Outline {self.pts}"

    def __eq__(self, other):
        return self.pts == other.pts


def proposed_next_pieces(current_orientation):
    dxdy = [(0, 1), (1, 0), (0, -1), (-1, 0)]
    for dx, dy in dxdy:
        yield current_orientation[0]+dx, current_orientation[1]+dy


def grow(_ap, extra_check = None):
    """Expand on a list of possibilities by adding one square to the end of each possibility at each free spot.

    Prevent the square from being in the squares reserved by extra_check. This is to leave room for the other leg.
    """
    if extra_check is None:
        extra_check = []
    new_all_possibilities = []
    for possibility in _ap:
        for piece in proposed_next_pieces(possibility[-1]):
            if piece not in possibility + extra_check:
                new_all_possibilities.append(possibility + [piece])
    return new_all_possibilities


def check_duplicates(possibility):
    return len(possibility) != len(set(possibility))


def flip_x_and_y(possibility: np.array):
    return possibility[:, ::-1]


def invert_y(possibility: np.array):
    return possibility * [1,-1]


def get_bounding_box(possibility) -> Tuple[int, int, int, int]:
    xs, ys = zip(*possibility)
    return min(xs), min(ys), max(xs), max(ys)


def translate_to_top_left(possibility):
    """translate the collection of squares so that the top left square is at (0,0)"""
    min_x, min_y, _, _ = get_bounding_box(possibility)
    return possibility - [min_x, min_y]


def group_by_outline(all_possibilities: List[np.array]) -> Dict[Outline,List[np.array]]:
    """group possibilities by their outline."""
    result: Dict[Outline,List[np.array]] = {}
    for possibility in all_possibilities:
        outline = Outline(possibility)
        result.setdefault(outline, []).append(possibility)
    return result


def get_xy_translate_range(possibility, x_left, x_right, y_bottom, y_top):
    """return an array containing all feasible translations of the man tiling.

    first index: translate along x
    second index: translate along y
    third index: index (element) of the tiling (0-13)
    fourth index: x or y coordinate of the square (0-1)
    """
    min_x, min_y, max_x, max_y = get_bounding_box(possibility)
    translate_x = np.insert(np.arange(x_left-min_x, x_right-max_x)[:,np.newaxis,np.newaxis,np.newaxis], 1, 0, axis=3)
    translate_y = np.insert(np.arange(y_bottom-min_y, y_top-max_y)[np.newaxis,:,np.newaxis,np.newaxis], 0, 0, axis=3)
    return translate_x+translate_y+possibility[np.newaxis,np.newaxis,:,:]


class GrouperOfMen:
    """Object to hold in memory all possible foldings and rotations of a Tile.
    Currently only implemented for 3-2-3 men tiles.
    Can translate all these tiles to desired x,y coordinates."""

    def __init__(self):
        self.all_p = self.get_all_men()
        self.grouped_by_outline = group_by_outline(self.all_p)

    def get_all_men_plus_translations(self, x_left, x_right, y_bottom, y_top):
        """Return all possible foldings and rotations, and on top of that, all possible rotations
        of each of those foldings that lie inside the bounding box defined by x_left, x_right, y_bottom, y_top.

        The return value is a dictionary with the outlines as keys, and the values are lists of numpy arrays.
        The structure of the numpy arrays is as follows:
        first index: translate along x
        second index: translate along y
        third index: index (element) of the tiling (0-13 in the case of 3-2-3 man tiles)
        fourth index: x or y coordinate of the square (0-1)
        """

        all_translations = {
            k: get_xy_translate_range(v[0], x_left, x_right, y_bottom, y_top)
            for k,v in self.grouped_by_outline.items()
        }

        return all_translations, self.grouped_by_outline

    @classmethod
    def get_all_men(cls):
        """Return all possible foldings and rotations of a 3-2-3 man tile"""
        all_combs = [np.array(p, dtype=int) for p in cls.all_combinations_of_man()]
        all_possibilities_with_symmetries = []

        def id_fn(inp):
            return inp

        for possibility in all_combs:
            # there are four possible symmetry operations: identity, flip x and y, invert y, and both
            for f1 in [id_fn, flip_x_and_y]:
                for f2 in [id_fn, invert_y]:
                    all_possibilities_with_symmetries.append(translate_to_top_left(f1(f2(possibility))))

        print(f"Found {len(all_possibilities_with_symmetries)} possibilities for the man tiles")

        return all_possibilities_with_symmetries

    @staticmethod
    def all_combinations_of_man() -> List[List[Tuple[int]]]:
        """Return all possible foldings of a 3-2-3 man, standing up with his head on top."""
        start_orientation = [(0, 0), (0, 1), (0, 2)]
        # add final torso piece:
        all_possibilities = [start_orientation + hip for hip in [[(-1, 2)], [(0, 3)], [(1, 2)]]]
        # add left hand:
        all_possibilities = [p + [(1, 1)] for p in all_possibilities]
        all_possibilities = grow(all_possibilities)
        # add right arm piece:
        for possibility in all_possibilities:
            possibility.append((-1, 1))
        # add right hand:
        all_possibilities = grow(all_possibilities)
        new_all_possibilities = []
        for possibility in all_possibilities:
            # get squares eligible for start of legs
            leg_starts = [p for p in proposed_next_pieces(possibility[3]) if p not in possibility]
            for h in range(len(leg_starts) - 1):
                left_leg_start = leg_starts[h]
                new_possibility = possibility + [left_leg_start]
                for j in range(h + 1, len(leg_starts)):
                    right_leg_start = leg_starts[j]
                    subset_possibilities = grow([new_possibility], extra_check=[right_leg_start])
                    subset_possibilities = grow(subset_possibilities, extra_check=[right_leg_start])
                    subset_possibilities = [p + [right_leg_start] for p in subset_possibilities]
                    for _ in [0, 1]:
                        subset_possibilities = grow(subset_possibilities)
                    new_all_possibilities += subset_possibilities
        return new_all_possibilities

    @property
    def n_squares_per_tile(self):
        return len(self.all_p[0])


def test_grouper():
    """Just a sanity check to see if the translations reach but do not exceed the bounding box."""
    grouper = GrouperOfMen()

    bbox = (0,49,0,24)  # x_left, x_right, y_bottom, y_top

    all_translations, grouped_by_outline = grouper.get_all_men_plus_translations(*bbox)

    minx = min([min(p[:, :, :, 0].flatten()) for p in all_translations])
    maxx = max([max(p[:, :, :, 0].flatten()) for p in all_translations])
    miny = min([min(p[:, :, :, 1].flatten()) for p in all_translations])
    maxy = max([max(p[:, :, :, 1].flatten()) for p in all_translations])

    if not (minx, maxx, miny, maxy) == bbox:
        print(f"Bounding box: {minx, miny, maxx, maxy}")
        raise ValueError("Dimensions of output of translations are not as expected")



