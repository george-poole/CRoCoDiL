from collections.abc import Iterable
from typing import Any

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.cm import ScalarMappable
from lucifex.fdm import ConstantSeries, NumericSeries, GridSeries, FunctionSeries
from lucifex.viz import (
    plot_line,
    plot_colormap, create_mosaic_figure, set_axes)
from lucifex.io.post import co_postprocess
from lucifex.utils import StrSlice, as_index

from .plot import plot_horizontal_average
from .tex import TeX


@co_postprocess
def co_plot_colormaps_mosaic(
    u: Iterable[GridSeries | FunctionSeries],
    t: int | float,
    colorbar: bool | tuple[float, float] = True,
    n_cols: int = 2,
    titles: Iterable[str | None] | None = None,
    suptitle: str | None = None,
    figscale: float = 1.0,
    width_ratio: float = 0.025,
    **kwargs,
) -> Figure:
    
    u = [GridSeries.from_series(i) if isinstance(i, FunctionSeries) else i for i in u]
    assert len(u) % n_cols == 0
    n_rows = len(u) // n_cols

    if titles is None:
        titles = [None] * len(u)

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
    
    for cmap_ax, ui, title in zip(cmap_axs, u, titles, strict=True):
        time_index = as_index(ui.time_series, t)
        plot_colormap(fig, cmap_ax, (ui.series[time_index], *ui.axes), title=title, colorbar=False, **kwargs) 
        if isinstance(colorbar, tuple):
            cmap: ScalarMappable = cmap_ax.collections[0]
            cmap.set_clim(*colorbar)
            if cmap_ax is cmap_axs[n_cols - 1]:
                fig.colorbar(cmap, ax=cmap_ax)

    if colorbar is True:
        cbar_axs: list[Axes] = list(axs.flat[1::2])
        for cmap_ax, cbar_ax in zip(cmap_axs, cbar_axs):
            fig.colorbar(cmap_ax.collections[0], cbar_ax)
    
    return fig


@co_postprocess
def co_plot_horizontal_average_mosaic(
    u: Iterable[GridSeries],
    t: Iterable[range | list[int | float] | int | StrSlice],
    n_cols: int = 2,
    titles: Iterable[str | None] | None = None,
    x_lims: tuple[float, float] | None = None,
    y_lims: tuple[float, float] | None = None,
    u_label: str | None = None,
    figscale: float = 1.0,
    width_ratio: float = 0.1,
    **kwargs,
) -> Figure:
    assert len(u) == len(t)
    u = [GridSeries.from_series(i) if isinstance(i, FunctionSeries) else i for i in u]  
    assert len(u) % n_cols == 0
    n_rows = len(u) // n_cols

    if titles is None:
        titles = [None] * len(u)

    subplot_kws = {'width_ratios': np.array([(1, width_ratio)] * n_cols).flatten()}
    n_cols = 2 * n_cols

    fig, axs = create_mosaic_figure(
        n_rows,
        n_cols,
        figscale=figscale,
        tex=True,
        **subplot_kws,

    )
    cmap_axs: list[Axes] = list(axs.flat[0::2])

    for cmap_ax, ui, ti, title in zip(cmap_axs, u, t, titles):
        plot_horizontal_average(
            ui, ti, x_lims, y_lims, u_label, title=title,
            fig_ax=(fig, cmap_ax), **kwargs)
        
    cpad_axs: list[Axes] = list(axs.flat[1::2])
    for cpad_ax in cpad_axs:
        cpad_ax.set_visible(False)

    return fig


@co_postprocess
def co_plot_timeseries(
    u: Iterable[NumericSeries | ConstantSeries],
    legend_labels: Iterable[Any] = (),
    legend_title: str | Iterable[str] = (),
    **kwargs,
) -> tuple[Figure, Axes]:
    u = [NumericSeries.from_series(i) if isinstance(i, ConstantSeries) else i for i in u]
    if not isinstance(legend_title, str):
        legend_title = ', '.join(legend_title)
    return plot_line(
        [(i.time_series, i.series) for i in u], 
        legend_labels, 
        legend_title, 
        x_label=TeX.T, 
        **kwargs,
    )


@co_postprocess
def co_plot_timeseries_mosaic(
    u: Iterable[NumericSeries | ConstantSeries],
    n_cols: int = 2,
    titles: Iterable[str | None] | None = None,
    suptitle: str | None = None,
    figscale: float = 1.0,
    **kwargs,
) -> Figure:
    u = [NumericSeries.from_series(i) if isinstance(i, ConstantSeries) else i for i in u]
    assert len(u) % n_cols == 0
    n_rows = len(u) // n_cols

    if titles is None:
        titles = [None] * len(u)

    fig, axs = create_mosaic_figure(
        n_rows,
        n_cols,
        suptitle,
        figscale,
        tex=True,
    )

    for ax, ui, title in zip(axs.flat, u, titles, strict=True):
        plot_line(fig, ax, (ui.time_series, ui.series), title=title, **kwargs) 

    return fig


@co_postprocess
def co_plot_scalars(
    scalars: list[float],
    xy_params: list[tuple[float, float]],
    cbar: str | None = None,
    cmap: str = 'viridis',
    suptitle: str | None = None,
    **kwargs,
) -> tuple[Figure, Axes]:
    
    x_axis = np.sort(np.unique([x for x, _ in xy_params]))
    y_axis = np.sort(np.unique([y for _, y in xy_params]))
    mat = np.zeros((len(x_axis), len(y_axis)))

    _filled = []
    for s, (x, y) in zip(scalars, xy_params, strict=True):
        i = np.argmin(np.abs(x_axis - x)) 
        j = np.argmin(np.abs(y_axis - y))
        assert (i, j) not in _filled
        _filled.append((i, j))
        mat[i, j] = s

    assert len(_filled) == mat.shape[0] * mat.shape[1]

    fig, ax = plt.subplots()
    set_axes(ax, **kwargs, tex=True)
    im = ax.imshow(mat.T, cmap=cmap, origin='lower')
    colorbar = fig.colorbar(im, ax=ax)
    if cbar:
        colorbar.set_label(cbar, rotation=360, ha='left')
    ax.set_xticks(np.arange(len(x_axis)))
    ax.set_yticks(np.arange(len(y_axis)))
    ax.set_xticklabels(x_axis)
    ax.set_yticklabels(y_axis)

    if suptitle:
        fig.suptitle(suptitle)

    return fig, ax

     



# @co_postprocess
# def co_plot_xy_scatter(
#     x: Iterable[float],
#     y: Iterable[float],
#     **kwargs,
# ) -> tuple[Figure, Axes]:
#     return plot_xy_scatter((x, y), **kwargs)


# @co_postprocess
# def co_plot_xyz_scatter(
#     x: Iterable[float],
#     y: Iterable[float],
#     z: Iterable[float],
#     **kwargs,
# ) -> tuple[Figure, Axes]:
#     return plot_xyz_scatter((x, y, z), **kwargs)

# tmax of umax

# make_tmax = lambda *umax: [compute_events_umax(i)[1] for i in umax]

# co_plot_xy_scatter(directories)(
#     ld(('Ra', PARAMETERS)),
#     ld(make_tmax, (GridSeries, 'umax', NUMERIC_SERIES)),
#     x_label='Ra',
#     y_label='$t_{\mathrm{max}}$',
# )

# co_plot_xyz_scatter(directories)(
#     ld(('Ra', PARAMETERS)),
#     ld(('Da', PARAMETERS)),
#     ld(make_tmax, (GridSeries, 'umax', NUMERIC_SERIES)),
#     x_label='Ra',
#     y_label='Da',
#     y_label='$t_{\mathrm{max}}$',
# )




# UMAX_TEX = '\max|\mathbf{u}|'

# @postprocess_many_to_one
# @optional_save
# def plot_umax_Nx_Ny(
#     umax: Iterable[SimulationData[ConstantSeries]],
#     Nx: Iterable[ParameterData[int]],
#     Ny: Iterable[ParameterData[int]],
#     title: str | None = None,
# ) -> tuple[Figure, Axes]: 
#     return _plot_timeseries(umax, {'N_x': Nx, 'N_y': Ny}, y_label=UMAX_TEX, title=title)


# @postprocess_many_to_one
# @optional_save
# def plot_umax_supg(
#     umax: Iterable[SimulationData[ConstantSeries]],
#     supg: Iterable[ParameterData[str]],
#     title: str | None = None,
# ) -> tuple[Figure, Axes]:
#     return _plot_timeseries(umax, {'/SUPG/': supg}, y_label=UMAX_TEX, title=title)


# @postprocess_many_to_one
# @optional_save
# def plot_cmax_Nx_Ny(
#     cmax: Iterable[SimulationData[ConstantSeries]],
#     Nx: Iterable[ParameterData[int]],
#     Ny: Iterable[ParameterData[int]],
#     title: str | None = None,
# ) -> tuple[Figure, Axes]:
#     return _plot_timeseries(cmax, {'N_x': Nx, 'N_y': Ny}, y_index=1, y_label='\max(c)', title=title)

