from typing import Callable

import numpy as np
from lucifex.fdm import FunctionSeries, ConstantSeries, NumericSeries, GridSeries
from lucifex.utils import as_index, as_indices
from lucifex.io.post import postprocess

from .contour import (
    contour_coordinates, 
    contour_peak_dimensions, contour_peaks,filter_trajectories,
    contour_arclength, contour_peak_mean_velocity, contour_peak_trajectories)
from ..pde.utils import spatial_average, grid_vertical_partition

when_lt = lambda s, t, val: np.array(t)[np.where(np.array(s) <= val)[0]]
when_gt = lambda s, t, val: np.array(t)[np.where(np.array(s) >= val)[0]]
when_first_gt = lambda s, t, val: when_gt(s, t, val)[0] if len(when_gt(s, t, val)) else None
when_last_gt = lambda s, t, val: when_gt(s, t, val)[-1] if len(when_gt(s, t, val)) else None
when_eq = lambda s, t, val: np.array(t)[np.isclose(s, val)]
when_first_eq = lambda s, t, val: when_eq(s, t, val)[0]


@postprocess
def compute_events_umax(
    umax: NumericSeries | ConstantSeries,
    onset: float,
    shut: float,
) -> tuple[float | None, float, float | None]:
    """
    `t_onset = argmin ₜ {umax(t) | umax > onset}`\\
    `t_max = argmax ₜ {umax(t)}`\\
    `t_shut = argmax ₜ {umax(t) | umax > shut}`

    `t_onset` < `t_max` < `t_shut`
    """
    if isinstance(umax, ConstantSeries):
        umax = NumericSeries.from_series(umax)
    t_onset = when_first_gt(umax.series, umax.time_series, onset)
    t_max = when_first_eq(umax.series, umax.time_series, np.max(umax.series))
    t_shut = when_last_gt(umax.series, umax.time_series, shut)
    return t_onset, t_max, t_shut


@postprocess
def compute_contour_series(
    f: GridSeries,
    window: tuple[float | int, float | int],
    alpha: float,
    adaptive: bool = False,
    name: str | None = None,
) -> NumericSeries:
    if isinstance(f, FunctionSeries):
        f = GridSeries.from_series(f)
    x, y = f.axes
    ti_min, ti_max = as_indices(f.time_series, window, window=True)

    contour_series = []
    time_series = []

    for si, ti in zip(f.series[ti_min: ti_max + 1], f.time_series[ti_min: ti_max + 1]):
        contour_xy = contour_coordinates((si, x, y), alpha, adaptive)
        contour_series.append(np.array(contour_xy))
        time_series.append(ti)

    return NumericSeries(contour_series, time_series, name)


def _compute_from_contour_series(
    func: Callable[[tuple[np.ndarray, np.ndarray]], int | float | np.ndarray],
    s: GridSeries,
    window: tuple[float | int, float | int],
    alpha: float,
    adaptive: bool = False,
    name: str | None = None,
) -> NumericSeries:
    contour = compute_contour_series(s, window, alpha, adaptive)
    return NumericSeries([func(i) for i in contour.series], contour.time_series, name)


@postprocess
def compute_front_average_height(
    s: GridSeries,
    s_value: float,
    window: tuple[float | int, float | int],
    name: str | None = None,
) -> NumericSeries:
    return _compute_from_contour_series(
        lambda i: np.mean(i[1]), 
        s, 
        window,
        s_value,
        name,
    )


@postprocess
def compute_front_sdev_height(
    s: GridSeries,
    window: tuple[float | int, float | int],
    s_value: float,
    name: str | None = None,
) -> NumericSeries:
    return _compute_from_contour_series(
        lambda i: np.std(i[1]), 
        s, 
        window,
        s_value,
        name,
    )


@postprocess
def compute_front_arclength(
    s: GridSeries,
    window: tuple[float | int, float | int],
    s_value: float,
    name: str | None = None,
):
    return _compute_from_contour_series(
        lambda i: contour_arclength(i),
        s,
        window,
        s_value,
        name,
    )


@postprocess
def compute_plume_number(
    c: GridSeries,
    window: tuple[float, float],
    alpha,
    adaptive,
    negative: bool = True,
    name: str | None = None,
) -> NumericSeries:
    contour = compute_contour_series(c, window, alpha, adaptive)
    # peaks = [contour_peaks((i, *c.axes), alpha, adaptive, negative) for i in c.series[ti_min: ti_max + 1]]
    peaks = [contour_peaks(i, negative=negative) for i in contour.series]
    number = [len(p) for p in peaks]
    return NumericSeries(number, contour.time_series, name)


@postprocess
def compute_plume_velocity_average(
    c: GridSeries,
    window: tuple[float, float],
    alpha: float,
    adaptive: bool,
    negative: bool = True,
    condition: int | None = None,
    name: str | None = None,
) -> NumericSeries:
    
    contour = compute_contour_series(c, window, alpha, adaptive)
    trajs, traj_times = contour_peak_trajectories(
        # [(i, x, y) for i in c.series[ti_min: 1 + ti_max]], 
        # c.time_series[ti_min: 1 + ti_max],
        contour.series,
        contour.time_series,
        negative=negative,
    )
    trajs, traj_times = filter_trajectories(trajs, traj_times, condition)
    _, uy_bar, t = contour_peak_mean_velocity(trajs, traj_times)
    return NumericSeries(uy_bar, t, name)


@postprocess
def compute_plume_width_average(
    c: GridSeries,
    window: tuple[float, float],
    alpha,
    adaptive,
    h0: float,
    fraction: float = 0.5,
    negative: bool = True,
    name: str | None = None,
) -> tuple[NumericSeries, NumericSeries, NumericSeries, NumericSeries, NumericSeries]:

    contour = compute_contour_series(c, window, alpha, adaptive)

    wbar_series = []
    for cntr in zip(contour.series):
        widths, *_ = contour_peak_dimensions(cntr, h0, fraction, negative=negative)
    wbar_series.append(np.mean(widths))

    return  NumericSeries(wbar_series, contour.time_series, name)


@postprocess
def compute_plume_wavenumber(
    c: GridSeries,
    window: tuple[float, float],
    alpha,
    adaptive,
    name: str | None = None,
) -> NumericSeries:
    contour = compute_contour_series(c, window, alpha, adaptive)
    fourier = lambda x: x
    NumericSeries([fourier(i) for i in contour.series], contour.time_series, name)
    raise NotImplementedError


@postprocess
def compute_horizontal_average(
    c: GridSeries,
    name: str | None = None,
) -> GridSeries:
    """
    `c(y,t) = ⟨c(x,y,t)⟩ₓ`
    """
    return GridSeries(
        [spatial_average(i, 'x') for i in c.series], 
        c.time_series, 
        (c.axes[0], ), 
        name,
    )


# FIXME
@postprocess
def compute_subdomain_averages(
    c: GridSeries,
    h: float | NumericSeries,
):
    """
    `c⁺(t)`, `c⁻(t)` averages in upper and lower subdomains partitioned by either `h₀`, `h(t)` or `h(x,t)`
    """
    if isinstance(h, float):
        y_index = as_index(c.axes[1], h)
        cplus = [np.mean(i[:, y_index:]) for i in c.series]
        cminus = [np.mean(i[:, :y_index]) for i in c.series]
    else:
        if h.shape == ():
            y_indices = [as_index(c.axes[1], i) for i in h.series]
            cplus = [np.mean(i[:, y:]) for i, y in zip(c.series, y_indices)]
            cminus = [np.mean(i[:, :y]) for i, y in zip(c.series, y_indices)]
        else:
            x, y = c.axes
            partitions = [grid_vertical_partition(i, x, y, cntr) for i, cntr in zip(c.series, h.series)]
            cplus = [np.mean(i[0]) for i in partitions]
            cminus = [np.mean(i[1]) for i in partitions]

    cplus = NumericSeries([], c.time_series, f'{c.name}plus')
    cminus = NumericSeries([], c.time_series, f'{c.name}minus')

    return cplus, cminus






# def compute_c_lower_tstar(
#     c: FunctionSeries,
#     onset: float = 0.01,
#     saturation: float = 0.8,
# )-> dict[str, float]:
#     """
#     `c(x,y=0, t)`
#     """
#     return _compute_cross_section_tstar(c, 0, onset, saturation)


# def compute_c_mid_timestamps(
#     c: FunctionSeries,
#     onset: float = 0.01,
#     saturation: float = 0.8,
# ) -> dict[str, float]:
#     """
#     `c(x,y=0.5, t)`
#     """
#     return _compute_cross_section_tstar(c, 0.5, onset, saturation)


# def compute_spurious_timestamps(
#     cminmax: ConstantSeries,
#     Sminmax: ConstantSeries,
#     Sr: float,
# ) -> dict[str, np.ndarray]:
#     cmin = [i[0] for i in cminmax.value_series]
#     cmax = [i[1] for i in cminmax.value_series]
#     Smin = [i[0] for i in Sminmax.value_series]
#     Smax = [i[1] for i in Sminmax.value_series]
#     return {
#         'cmin_spurious': when_lt(cmin, cminmax.time_series, 0),
#         'cmax_spurious': when_gt(cmax, cminmax.time_series, 1),
#         'Smin_spurious': when_lt(Smin, cminmax.time_series, 0),
#         'Smax_spurious': when_gt(Smax, cminmax.time_series, Sr),
#     }


# # def _compute_cross_section_tstar(
#     c: FunctionSeries,
#     y_value: float,
#     onset: float,
#     saturation: float,
# )-> dict[str, float]:
#     """
#     `c(x,y=yᵢ, t)`
#     """

#     axes = grid(use_cache=True)(c.mesh)
#     y_index = np.argmin(np.abs(axes[1] - y_value))
#     y_value = axes[1][y_index]

#     grid_series = [grid(use_cache=True)(_c)[0] for _c in c.series]
#     max_series = np.array([np.max(_g[:, y_index]) for _g in grid_series])
#     min_series = np.array([np.min(_g[:, y_index]) for _g in grid_series])

#     timestamps = {
#         f't_{c.name}(y={y_value:.2f})_onset': when_first_gt(min_series, c.time_series, onset),
#         f't_{c.name}(y={y_value:.2f})_saturation': when_first_gt(max_series, c.time_series, saturation)
#     }

#     return timestamps