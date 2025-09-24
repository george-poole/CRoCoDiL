from typing import Callable
import numpy as np

from dolfinx.fem import Function

from lucifex.utils import grid


def heaviside(
    fx: Callable[[np.ndarray], np.ndarray],
    f_plus: float = 1.0,
    f_minus: float = 0.0,
    eps: float | tuple[float, float] | None = None,
):
    """
    `H(f(x))` \\
    `= f₊` if `f(x) > 0` \\
    `= f₋` otherwise

    Optionally to smooth out discontinuities

    `H(f(x))` \\
    `= f₊ tanh(f(x) / ϵ)` if `f(x) > 0` \\
    `= f₋` otherwise

    `H(f(x))` \\
    `= (f₊ - f₋) tanh(f(x) / ϵ₊) / 2 + (f₊ + f₋) / 2` if `f(x) > 0` \\
    `= (f₊ - f₋) tanh(f(x) / ϵ₋) / 2 + (f₊ + f₋) / 2` otherwise
    """
    ind = lambda x: (fx(x) >= 0)

    if isinstance(eps, float):
        return lambda x: f_plus * np.tanh(fx(x) / eps) * ind(x) + f_minus
    elif isinstance(eps, tuple):
        eps_lt, eps_gt = eps
        return lambda x: (
            (0.5 * (f_plus - f_minus) * np.tanh(fx(x) / eps_gt) + 0.5 * (f_plus + f_minus)) * ind(x)
            + 0.5 * (f_plus - f_minus) * np.tanh(fx(x) / eps_lt) + 0.5 * (f_plus + f_minus)
        )
    else:
        return lambda x: f_plus * ind(x) + f_minus


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
