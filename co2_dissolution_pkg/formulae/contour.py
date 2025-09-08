from typing import overload, Callable
from operator import lt, gt

import numpy as np
from skimage import measure
from scipy.signal import find_peaks
from dolfinx.fem import Function
from scipy.interpolate import interp1d

from lucifex.utils import grid, optional_lru_cache


@optional_lru_cache
def contour_coordinates(
    fxy: Function | tuple[np.ndarray, np.ndarray, np.ndarray],
    alpha: float,
    adaptive: bool = False,
    primary: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns `([x₀, x₁, x₂, ...], [y₀, y₁, y₂, ...])`

    with orderered `x₀ < x₁ < x₂ < ...` for the contour defined by either an adaptive level

    `Γ = {(x,y) ∈ Ω : c = α min(c)  + (1 − α) max(c)}`

    or a fixed level

    `Γ = {(x,y) ∈ Ω : c = α}`
    """
    if isinstance(fxy, Function):
        f = grid(fxy)
        x_axis, y_axis = grid(use_cache=True)(fxy.function_space.mesh)
    else:
        f, x_axis, y_axis = fxy

    def _coordinates(contour: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        x_indices, y_indices = contour[:, 0], contour[:, 1]
        x_scale = (np.max(x_axis) - np.min(x_axis)) / (len(x_axis) - 1)
        y_scale = (np.max(y_axis) - np.min(y_axis)) / (len(y_axis) - 1)
        x_coordinates = x_scale * x_indices
        y_coordinates = y_scale * y_indices
        if x_coordinates[0] > x_coordinates[1]:
            x_coordinates = x_coordinates[::-1]
            y_coordinates = y_coordinates[::-1]
        return x_coordinates, y_coordinates
    
    if adaptive:    
        if alpha < 0 or alpha > 1:
            raise ValueError(f'`alpha` must be in interval [0, 1], not {alpha}.')     
        fmin = np.min(f)
        fmax = np.max(f)
        alpha = alpha * fmin + (1 - alpha) * fmax

    contours = measure.find_contours(f, alpha)

    if len(contours) == 0:
        raise ValueError(f'No contours found at level {alpha}.')

    if primary:
        sizes = [len(i) for i in contours]
        sizes_max = max(sizes)
        primary_contour = contours[sizes.index(sizes_max)]
        return _coordinates(primary_contour)
    else:
        x_coordinates, y_coordinates = [], []
        sorted_contours = sorted(contours, key=lambda c: _coordinates(c)[0][0], reverse=True)
        for contour in sorted_contours:
            x, y = _coordinates(contour)
            x_coordinates.extend(x)
            y_coordinates.extend(y)
        return np.array(x_coordinates), np.array(y_coordinates)


@overload
def contour_arclength(
    contour: tuple[np.ndarray, np.ndarray],
) -> float:
    ...


@overload
def contour_arclength(
    fxy: Function | tuple[np.ndarray, np.ndarray, np.ndarray], 
    alpha: float,
    adaptive: bool = False,
    primary: bool = True,
) -> float:
    ...
    

def contour_arclength(
    arg,
    alpha = None,
    adaptive = False,
    primary = True,      
):
    """
    `len(h) = `
    """
    if isinstance(arg, Function) or len(arg) == 3:
        assert alpha is not None
        use_cache = True if isinstance(arg, Function) else False
        x, y = contour_coordinates(use_cache=use_cache)(arg, alpha, adaptive, primary)
    else:
        x, y = arg

    return np.trapz(
            np.sqrt(1 + np.gradient(y, x) ** 2),
            x,
        )


@overload
def contour_peaks(
    contour: tuple[np.ndarray, np.ndarray],
    *,
    negative: bool = False,
) -> list[tuple[float, float]]:
    ...


@overload
def contour_peaks(
    fxy: Function | tuple[np.ndarray, np.ndarray, np.ndarray], 
    alpha: float,
    adaptive: bool = False,
    primary: bool = True,
    *,
    negative: bool = False,
) -> list[tuple[float, float]]:
    ...


def contour_peaks(
    arg,
    alpha = None,
    adaptive = False,
    primary = True,
    *,
    negative = False,       
):
    """
    Set of plume tip points defined by

    `P = {(xⱼ, yⱼ})ⱼ`
    """
    if isinstance(arg, Function) or len(arg) == 3:
        assert alpha is not None
        use_cache = True if isinstance(arg, Function) else False
        x, y = contour_coordinates(use_cache=use_cache)(arg, alpha, adaptive, primary)
    else:
        x, y = arg

    if negative:
        scale = -1.0
    else:
        scale = 1.0
    y_heights = scale * y
    peaks, _ = find_peaks(y_heights)
    return [(i, j) for (i, j) in zip(x[peaks], y[peaks], strict=True)]


@overload
def contour_peak_dimensions(
    contour: tuple[np.ndarray, np.ndarray],
    baseline: float,
    fraction: float,
    *,
    negative: bool = False,
) -> tuple[list[float], list[float], list[float], list[float]]:
    ...


@overload
def contour_peak_dimensions(
    peaks: list[tuple[float, float]],
    baseline: float,
    fraction: float,
) -> tuple[list[float], list[float], list[float], list[float]]:
    ...


@overload
def contour_peak_dimensions(
    fxy: Function | tuple[np.ndarray, np.ndarray, np.ndarray],
    baseline: float,
    fraction: float,
    alpha: float,
    adaptive: bool = False,
    primary: bool = True,
    *,
    negative: bool = False,
) -> tuple[list[float], list[float], list[float], list[float]]:
    ...


def contour_peak_dimensions(
    arg,
    baseline,
    fraction,
    alpha = None,
    adaptive: bool = False,
    primary: bool = True,
    *,
    negative = False,
):
    """
    Returns `widths, lengths, x_lefts, x_rights`
    """
    if isinstance(arg, Function) or (isinstance(arg, tuple) and len(arg) == 3):
        assert alpha is not None
        use_cache = True if isinstance(arg, Function) else False
        x, y = contour_coordinates(use_cache=use_cache)(arg, alpha, adaptive, primary)
        peaks = contour_peaks((x, y), negative=negative)
    elif isinstance(arg, tuple):
        x, y = arg
        peaks = contour_peaks((x, y), negative=negative)
    else:
        peaks = arg

    widths = []
    lengths = []
    lefts = []
    rights = []

    for (x_peak, y_peak) in peaks:
        if negative:
            assert y_peak < baseline
            y_target = baseline - fraction * (baseline - y_peak)
            comparison = lt
        else:
            assert y_peak > baseline
            y_target = baseline + fraction * (y_peak - baseline)
            comparison = gt

        x_peak_index = np.where(np.isclose(contour_coordinates[0], x_peak))[0][0]

        i = 0
        y_value_left = y_peak
        while comparison(y_value_left, y_target):
            try:
                y_value_left = contour_coordinates[1][x_peak_index - i]
                x_l = contour_coordinates[0][x_peak_index - i]
                i += 1
            except IndexError:
                break

        i = 0
        y_value_right = y_peak
        while comparison(y_value_right, y_target):
            try:
                y_value_right = contour_coordinates[1][x_peak_index + i]
                x_r = contour_coordinates[0][x_peak_index + i]
                i += 1
            except IndexError:
                break

        widths.append(x_r - x_l)
        lengths.append(0.5 * (y_value_left + y_value_right))
        lefts.append(x_l)
        rights.append(x_r)
        
    return widths, lengths, lefts, rights


@overload
def contour_peak_trajectories(
    fxy: list[Function | tuple[np.ndarray, np.ndarray, np.ndarray]],
    t: list[float],
    alpha: float,
    adaptive: bool = False,
    primary: bool = True,
    *,
    negative: bool = False,
    threshold: float = 0.1,
) -> tuple[
    list[list[tuple[float, float]]],
    list[list[float]]
]:
    ...


@overload
def contour_peak_trajectories(
    contours: list[tuple[np.ndarray, np.ndarray]],
    t: list[float],
    *,
    negative: bool = False,
    threshold: float = 0.1,
) -> tuple[
    list[list[tuple[float, float]]],
    list[list[float]]
]:
    ...


def contour_peak_trajectories(
    args: list[Function | tuple[np.ndarray, np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray]],
    t: list[float],
    alpha: float,
    adaptive: bool = False,
    primary: bool = True,
    *,
    negative: bool = False,
    threshold: float = 0.1,
) -> tuple[
    list[list[tuple[float, float]]],
    list[list[float]]
]:

    trajectories: list[list[tuple[float, float]]] = []
    times: list[list[float]] = []

    for i, (arg, t) in enumerate(zip(args, t, strict=True)):
        if isinstance(arg, Function) or len(arg) == 3:
            peaks = contour_peaks(arg, alpha, adaptive, primary, negative=negative)
        else:
            peaks = contour_peaks(arg, negative=negative)
        if i == 0:
            for p in peaks:
                trajectories.append([p])
                times.append([t])
        else:
            peaks_previous = [t[-1] for t in trajectories]
            for p in peaks:
                distances = [np.linalg.norm(np.array(p) - np.array(pp)) for pp in peaks_previous]
                distances_min = np.min(distances)
                if distances_min <= threshold:
                    traj_index = distances.index(distances_min)
                    trajectories[traj_index].append(p)
                    times[traj_index].append(t)
                else:
                    trajectories.append([p])
                    times.append([t])

    return trajectories, times


def filter_trajectories(
    trajs: list[list[tuple[float, float]]],
    times: list[list[float]],
    condition: int 
    | Callable[[list[tuple[float, float]]], bool] 
    | None = None,
) -> tuple[
    list[list[tuple[float, float]]],
    list[list[float]]
]:
    if condition is None:
        return trajs, times
    
    if isinstance(condition, int):
        _condition = lambda t: len(t) > condition
    else:
        _condition = condition
    
    trajs_filtered = []
    times_filtered = []
    
    for traj, time in zip(trajs, times, strict=True):
        if _condition(traj):
            trajs_filtered.append(traj)
            times_filtered.append(time)
    
    return trajs_filtered, times_filtered


def contour_peak_velocities(
    trajs: list[list[tuple[float, float]]],
    times: list[list[float]],
) -> list[tuple[np.ndarray, np.ndarray]]:
    """
    `[(ux, uy), ...]`
    """
    velocities = []
    for traj, time in zip(trajs, times, strict=True):
        x_traj = [xy[0] for xy in traj]
        y_traj = [xy[1] for xy in traj]
        ux = np.gradient(x_traj, time, edge_order=2)
        uy = np.gradient(y_traj, time, edge_order=2)
        velocities.append((ux, uy))
    return velocities


def contour_peak_mean_velocity(
    trajs: list[list[tuple[float, float]]],
    times: list[list[float]], 
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    `ux, uy, t`
    """
    velocities = contour_peak_velocities(trajs, times)
    n_peaks = len(velocities)

    ux_interp = []
    uy_interp = []
    for t, v in zip(times, velocities, strict=True):
        ux, uy = v
        ux_interp.append(interp1d(t, ux))
        uy_interp.append(interp1d(t, uy))

    ux_mean = lambda t: np.array(sum(u(t) for u in ux_interp)) / n_peaks
    uy_mean = lambda t: np.array(sum(u(t) for u in uy_interp)) / n_peaks

    time_mins = [min(t) for t in times]
    time_maxs = [max(t) for t in times]
    time_window = (max(time_mins), min(time_maxs))
    dt = np.min([np.min(np.diff(t)) for t in times])
    t = np.arange(time_window[0], time_window[1], dt)

    return ux_mean(t), uy_mean(t), t





    

