from collections.abc import Iterable
from typing import Any, Callable

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.cm import ScalarMappable
from matplotlib.animation import FuncAnimation
from lucifex.fdm import FunctionSeries, ConstantSeries, GridSeries, NumericSeries
from lucifex.viz import (plot_line, plot_colormap, plot_twins, 
                         plot_contours, animate_colormap, create_mosaic_figure)
from lucifex.io.post import postprocess
from lucifex.utils import StrSlice, as_index, as_indices

from .compute import compute_horizontal_average
from ..formulae.contour import contour_peak_trajectories, filter_trajectories, contour_coordinates


class TeX:
    X = 'x'
    """
    x-coordinate
    """
    Y = 'y'
    """
    y-coordinate
    """
    T = 't'
    """
    time
    """
    MIN_X = staticmethod(lambda u: f'\\min_{{\mathbf{{x}}}}({u})')
    """
    `minₓ(u)`
    """
    MAX_X = staticmethod(lambda u: f'\\max_{{\mathbf{{x}}}}({u})')
    """
    `maxₓ(u)`
    """
    ABS_MIN_X = staticmethod(lambda u: f'\\min_{{\mathbf{{x}}}}|\mathbf{{{u}}}|')
    """
    `minₓ|𝐮|`
    """
    ABS_MAX_X = staticmethod(lambda u: f'\\max_{{\mathbf{{x}}}}|\mathbf{{{u}}}|')
    """
    `maxₓ|𝐮|`
    """
    BRAKET = lambda a: f'\langle {a}\\rangle'
    """
    `⟨u⟩`
    """


@postprocess
def plot_colormaps(
    u: GridSeries | FunctionSeries,
    t: range | list[int | float] | int | StrSlice,
    contours: GridSeries | FunctionSeries | list[tuple[np.ndarray, np.ndarray]] | None = None,
    u_label: str | None = None,
    suptitle: str | None = None,
    **kwargs,
) -> list[tuple[Figure, Axes]]:
    """
    Plots `u(x, y, t)` for all `t ∈ {tᵢ}`
    """
    if isinstance(u, FunctionSeries):
        u = GridSeries.from_series(u)
    if u_label is None:
        u_label = u.name

    time_indices = as_indices(u.time_series, t)

    fig_axs: list[tuple[Figure, Axes]] = []
    for i in time_indices:
        title = f'{u_label}({TeX.X}, {TeX.Y}, {TeX.t}={u.time_series[i]:.6f})'
        fig, ax = plot_colormap((u.series[i], *u.axes), title=title, **kwargs)
        if suptitle:
            fig.suptitle(suptitle)
        fig_axs.append((fig, ax))
        
    if isinstance(contours, (GridSeries, FunctionSeries)):
        for (fig, ax), i in zip(fig_axs, time_indices):
            t = u.time_series[i]
            t_cntr_index = as_index(contours.time_series, t)
            t_cntr = contours.time_series[t_cntr_index]
            if np.isclose(t, t_cntr):
                plot_contours(fig, ax, contours.series[i], **kwargs)
            else:
                print(f'Contour timeseries has no data at time t={t}.')
    elif isinstance(contours, list):
        for (fig, ax), i in zip(fig_axs, time_indices):
            x, y = contours[i]
            plot_line(fig, ax, (x, y), **kwargs)

    return fig_axs


@postprocess
def plot_colormaps_mosaic(
    u: GridSeries | FunctionSeries,
    t: range | list[int | float] | int | StrSlice,
    # contours: GridSeries | FunctionSeries | list[tuple[np.ndarray, np.ndarray]] | None = None,   #TODO
    colorbar: bool | tuple[float, float] = True,
    n_cols: int = 2,
    suptitle: str | None = None,
    figscale: float = 1.0,
    width_ratio: float = 0.025,
    u_label: str | None = None,
    **kwargs,
) -> Figure:
    
    if isinstance(u, FunctionSeries):
        u = GridSeries.from_series(u)
    if u_label is None:
        u_label = u.name

    time_indices = as_indices(u.time_series, t)
    assert len(time_indices) % n_cols == 0
    n_rows = len(time_indices) // n_cols

    if colorbar is True:
        subplot_kws = {'width_ratios': np.array([(1, width_ratio)] * n_cols).flatten()}
        n_cols = 2 * n_cols
    else:
        subplot_kws = {}

    fig, axs = create_mosaic_figure(
        n_rows,
        n_cols,
        suptitle,
        figscale,
        tex=True,
        **subplot_kws,
    )

    if colorbar is True:
        cmap_axs: list[Axes] = list(axs.flat[0::2])
    else:
        cmap_axs = list(axs.flat)

    for cmap_ax, i in zip(cmap_axs, time_indices):
        title = f'{u_label}(x, y, t={u.time_series[i]:.6f})'
        plot_colormap(fig, cmap_ax, (u.series[i], *u.axes), title=title, colorbar=False, **kwargs)
        if isinstance(colorbar, tuple):
            cmap: ScalarMappable = cmap_ax.collections[0]
            cmap.set_clim(*colorbar)
            if cmap_ax is cmap_axs[n_cols - 1]:
                fig.colorbar(cmap, ax=cmap_ax)

    if colorbar is True:
        cbar_axs: list[Axes] = list(axs.flat[1::2])
        [fig.colorbar(cmap_ax.collections[0], cbar_ax) for cmap_ax, cbar_ax in zip(cmap_axs, cbar_axs)]

    return fig


@postprocess
def plot_timeseries(
    u: NumericSeries | ConstantSeries,
    label: str | None = None,
    t_events: Iterable[float] = (),
    t_labels: Iterable[str] = (),
    suptitle: str | None = None,
    **kwargs,
) -> tuple[Figure, Axes]:
    if isinstance(u, ConstantSeries):
        u = NumericSeries.from_series(u)
    assert u.shape == ()
    
    if label is None:
        label = u.name

    fig, ax = plot_line((u.time_series, u.series), x_label=TeX.T, y_label=label, **kwargs)

    if t_events:
        ylims = ax.get_ylim()
        ax.vlines(t_events, [ylims[0]], [ylims[1]], color='gray', linestyles='dashed')
        ax.set_ylim(*ylims)

    if t_labels:
        raise NotImplementedError
    
    if suptitle:
        fig.suptitle(suptitle)
    
    return fig, ax


@postprocess
def plot_twinned_timeseries(
    u: Iterable[NumericSeries | ConstantSeries] | NumericSeries, 
    labels: Iterable[str],
    individual: bool = True,
    together: bool = True,
    twins: bool = False,
    suptitle: str | None = None,
) -> list[tuple[Figure, Axes]]:

    if isinstance(u, ConstantSeries):
        u = NumericSeries.from_series(u)
    if isinstance(u, NumericSeries):
        assert len(u.shape) == 1
        u = [NumericSeries([i[j] for i in u.series], u.time_series) for j in range(u.shape[0])]
    u = [i if isinstance(i, NumericSeries) else NumericSeries.from_series(i) for i in u]

    fig_axs: list[tuple[Figure, Axes]] = []
    if individual:
        for i, l in zip(u, labels, strict=True):
            fig_ax = plot_line((i.time_series, i.series), x_label=TeX.T, y_label=l)
            fig_axs.append(fig_ax)

    if together:
        x_lims = (min(i.time_series[0] for i in u),
                  min(i.time_series[-1] for i in u))
        fig_ax = plot_line(
            [(i.time_series, i.series) for i in u],  x_lims=x_lims, x_label=TeX.T, legend_labels=labels)
        fig_axs.append(fig_ax)

    if twins:
        assert len(u) == 2
        assert len(labels) == 2
        fig_ax = plot_twins(
            (u[0].time_series, u[1].time_series), (u[0].series, u[1].series), labels,  x_label=TeX.T)
        fig_axs.append(fig_ax)

    if suptitle:
        [fig.suptitle(fig) for fig, _ in fig_axs]

    return fig_axs


@postprocess
def plot_horizontal_average(
    u: GridSeries,
    t: range | list[int | float] | int | StrSlice,
    x_lims: tuple[float, float] | None = None,
    y_lims: tuple[float, float] | None = None,
    u_label: str | None = None,
    title: str | None = None,
    legend: bool = True,
    fig_ax: tuple[Figure, Axes] | None = None,
    **kwargs,
):
    if u_label is None:
        u_label = u.name
    if fig_ax is None:
        fig_ax = ()
    y = u.axes[1]
    if y_lims is None:
        y_lims = y[0], y[-1]
    
    time_indices = as_indices(u.time_series, t)
    time_values = [u.time_series[i] for i in time_indices]
    u_braket = compute_horizontal_average(u)

    legend_kwargs = {}
    if legend:
        legend_kwargs = dict(
            legend_labels=(min(time_values), max(time_values)),
            legend_title=TeX.T,
        )

    return plot_line(
        *fig_ax,
        [(u_braket.series[i], y) for i in time_indices], 
        x_label=f'{TeX.BRAKET(u_label)}_x', 
        y_label='y', 
        x_lims=x_lims,
        y_lims=y_lims,
        cycler='jet',
        title=title,
        **legend_kwargs,
        **kwargs,
    )


@postprocess
def plot_finger_trajectories(
    c: GridSeries, 
    window: tuple[float, float],
    alpha: float,
    adaptive: bool = False,
    negative: bool = True,
    condition: int | Callable[[list[tuple[float, float]]], bool] | None = None,
    c_label: str | None = None,
    line_kwargs: dict[str, Any] | tuple[tuple[str, Any], ...] = (('color', 'cyan'), ('lw', '0.75')),
    scatter_kwargs: dict[str, Any] | tuple[tuple[str, Any], ...] = (('color', 'cyan'), ),
    **kwargs,
):
    if isinstance(c, FunctionSeries):
        c = GridSeries.from_series(c)
    if c_label is None:
        c_label = c.name
    if isinstance(line_kwargs, tuple):
        line_kwargs = dict(line_kwargs)
    if isinstance(scatter_kwargs, tuple):
        scatter_kwargs = dict(scatter_kwargs)

    x, y = c.axes

    ti_min, ti_max = as_indices(c.time_series, window, window=True)

    trajs, traj_times = contour_peak_trajectories(
        [(i, x, y) for i in c.series[ti_min: 1 + ti_max]], 
        c.time_series[ti_min: 1 + ti_max],
        alpha,
        adaptive,
        negative=negative,
    )
    trajs, traj_times = filter_trajectories(trajs, traj_times, condition)
    xy_contour = contour_coordinates(
        (c.series[ti_max], x, y), 
        alpha, 
        adaptive,
    )
    
    title = f'{c_label}({TeX.X}, {TeX.Y}, {TeX.T}={c.time_series[ti_max]:.3f})'
    fig, ax = plot_colormap((c.series[ti_max], x, y), title=title, **kwargs)
    ax.plot(*xy_contour, **line_kwargs)

    for traj in trajs:
        xy_final = traj[-1]
        ax.scatter(*xy_final, **scatter_kwargs)
        ax.plot([xy[0] for xy in traj], [xy[1] for xy in traj], **line_kwargs)

    return fig, ax


@postprocess
def animate_colormaps(
    u: GridSeries,
    colorbar: bool | tuple[float, float] = True,
    interval: int = 100,
    u_label: str | None = None,
    **anim_kwargs: Any,
) -> FuncAnimation:   
    if u_label is None:
        u_label = u.name
    x, y = u.axes
    cmap_series = [(i, x, y) for i in u.series]
    title_series = [f'{u_label}(x, y, t={t:.6f})' for t in u.time_series]
    return animate_colormap(cmap_series, title_series, colorbar=colorbar, interval=interval, **anim_kwargs)
