from functools import partial

import numpy as np

from lucifex.io import create_dir_path, find_dir_paths
from lucifex.fdm import GridFunctionSeries, NPyConstantSeries
from lucifex.mesh import GridMesh
from lucifex.fem import (
    average_grid, 
    resample_grid, 
    crop_grid, 
    mirror_grid, 
    cross_section_grid,
    copy_grid,
    GridFunction,
)
from lucifex.plt import (
    plot_colormap, plot_line, save_figure, create_multifigure,
    plot_colormap_multifigure, plot_line_multifigure, 
    plot_scatter, scatter_size, plot_contours, plot_quiver,
    plot_twin_lines, plot_stacked_lines,
    configure_matplotlib,
)
from lucifex.utils.npy_utils import as_index, derivative, resample
from lucifex.utils.py_utils import FrozenDict, as_int_if_close
from lucifex.sim import GridSimulationFromNPZ
from crocodil.dns.system_a import SYSTEM_A_REFERENCE
from crocodil.theory.system_a import (
    mass_dissolved_asymptote, mass_capillary_asymptote,
    mass_dissolved_initial, mass_capillary_initial,
)

from ipynb_utils import *


PARAMS_NUMERICAL = FrozenDict(
    c_stabilization=None,
    c_limits=True,
)
DIR_ROOT = create_dir_path(
    PARAMS_NUMERICAL, 
    dir_root='./',
    dir_prefix='data', 
    dir_params=PARAMS_NUMERICAL.keys(), 
)
DIR_FIGS = f'{DIR_ROOT}/figures'

T_STOP = 120.0

SIM_DIR_PATHS = find_dir_paths(
    DIR_ROOT, 
    include=f't_stop={T_STOP}_*',
    contains=('CHECKPOINT.h5', 'c.npz'),
)