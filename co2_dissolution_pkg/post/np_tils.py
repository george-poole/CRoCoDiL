from lucifex.fem import LUCiFExFunction as Function
from lucifex.utils import grid

import numpy as np


def spatial_average(
    f: Function  | np.ndarray,
    axis: str | int | None = None,
    slc: slice | tuple[slice, slice] | None = None, # TODO import as_slice
) -> np.ndarray:

    if isinstance(axis, str):
        axis = ('x', 'y', 'z').index(axis)

    if isinstance(f, Function):
        f = grid(f)

    if isinstance(slc, slice):
        f = f[slc]

    if isinstance(slc, tuple):
        f = f[slc[0], slc[1]]

    return np.mean(f, axis)


def grid_vertical_partition(
    f: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    contour: tuple[np.ndarray, np.ndarray],
    mask = np.nan,
) -> tuple[np.ndarray, np.ndarray]:
    """
    NOTE assumes single-valued contour spanning the grid's entire width
    """
        
    x_contour, y_contour = contour
    assert np.isclose(x_axis[0], min(x_contour))
    assert np.isclose(x_axis[-1], max(x_contour))
    # TODO assert single-valued contour

    grid_above = f.copy()
    grid_below = f.copy()
    
    for i, j in zip(range(f.shape[0]), range(f.shape[1])):
        x = x_axis[i]
        y = y_axis[i]
        contour_index = np.argmin(np.abs(x_contour - x))
        if y >= y_contour[contour_index]:
            grid_below[i, j] = mask
        else:
            grid_above[i, j] = mask

    return grid_above, grid_below
