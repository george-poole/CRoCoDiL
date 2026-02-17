import numpy as np
from typing import Iterable, overload

from lucifex.fem import Function
from lucifex.utils import grid, as_index


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
