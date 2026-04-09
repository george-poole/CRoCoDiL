import os
from typing import Callable

import numpy as np
from mpi4py import MPI
from dolfinx.mesh import Mesh

from lucifex.mesh import MeshBoundary, mesh_boundary, rectangle_mesh
from lucifex.solver import OptionsPETSc
from lucifex.pde.scaling import ScalingOptions
from lucifex.utils.py_utils import FloatEnum


CROCODIL_JIT_DIR = os.path.abspath(
    os.path.join(
        __file__,
        '../../..',
        '__jit__',
    )
)


SCALINGS = ScalingOptions(
    ('Ad', 'Di', 'Ki', 'Bu', 'X'),
    lambda Ra, Da=0: {
        'advective': (1, 1/Ra, Da, 1, 1),
        'diffusive': (1, 1, Ra * Da, Ra, 1),
        'advective_diffusive': (1, 1, Da/Ra, 1, Ra),
        'reactive': (1, 1, 1, np.sqrt(Ra / Da) if Da else np.inf, np.sqrt(Ra * Da)),
    }
)
"""
Choice of length scale `ℒ`, velocity scale `𝒰`
and time scale `𝒯` in the non-dimensionalization.

`'advective'` \\
`ℒ` = domain size \\
`𝒰` = advective speed

`'diffusive'` \\
`ℒ` = domain size \\
`𝒰` = diffusive speed

`'advective_diffusive'` \\
`ℒ` = diffusive length \\
`𝒰` = advective speed

`'reactive'` \\
`ℒ` = diffusive length \\
`𝒯` = reactive time
"""


class CONVECTION_CONSTANTS(FloatEnum):
    """
    Constants from the theory of convection in porous media
    """

    RA_CRITICAL = 4 * np.pi **2
    """
    Critical Rayleigh for the onset of Rayleigh-Benard convection
    """
    FLUX_FACTOR = 0.008
    """
    https://doi.org/10.1038/s43247-025-02554-0
    """



def heaviside(
    fx: Callable[[np.ndarray], np.ndarray],
    f_plus: float = 1.0,
    f_minus: float = 0.0,
    eps: float | tuple[float, float] | None = None,
):
    """
    `H(f(x))` \\
    `= f₊` if `f(x) >= 0` \\
    `= f₋` otherwise

    Optionally to smooth out discontinuities

    `H(f(x))` \\
    `= f₊ tanh(f(x) / ϵ)` if `f(x) >= 0` \\
    `= f₋` otherwise

    `H(f(x))` \\
    `= (f₊ - f₋) tanh(f(x) / ϵ₊) / 2 + (f₊ + f₋) / 2` if `f(x) >= 0` \\
    `= (f₊ - f₋) tanh(f(x) / ϵ₋) / 2 + (f₊ + f₋) / 2` otherwise
    """
    geq = lambda x: (fx(x) >= 0)
    le = lambda x: (fx(x) < 0)

    if isinstance(eps, float):
        return lambda x: 0.0 + f_plus * np.tanh(fx(x) / eps) * geq(x) + f_minus * le(x)
    elif isinstance(eps, tuple):
        eps_lt, eps_gt = eps
        return lambda x: (
            0.0
            + (0.5 * (f_plus - f_minus) * np.tanh(fx(x) / eps_gt) + 0.5 * (f_plus + f_minus)) * geq(x)
            + (0.5 * (f_plus - f_minus) * np.tanh(fx(x) / eps_lt) + 0.5 * (f_plus + f_minus)) * le(x)
        )
    else:
        return lambda x: 0.0 + f_plus * geq(x) + f_minus * le(x)


def rectangle_mesh_closure(
    Lx: float,
    Ly: float,
    Nx: int,
    Ny: int,
    cell: str,
    name: str = 'Omega',
    clockwise_names: tuple[str, str, str, str] = ('upper', 'right', 'lower', 'left'),
    comm: MPI.Comm | str = MPI.COMM_WORLD,
) -> tuple[Mesh, MeshBoundary]:
    """
    `Ω ∪ ∂Ω`
    """
    mesh = rectangle_mesh(Lx, Ly, Nx, Ny, cell, name, comm)
    boundary = mesh_boundary(
        mesh,
        {
            clockwise_names[0]: lambda x: x[1] - Ly,
            clockwise_names[1]: lambda x: x[0] - Lx,
            clockwise_names[2]: lambda x: x[1],
            clockwise_names[3]: lambda x: x[0],
        },
    )
    return mesh, boundary


def use_streamfunction(
    flow_petsc: OptionsPETSc
    | tuple[OptionsPETSc, OptionsPETSc | None],   
) -> bool:
    """
    `flow_petsc = (psi_petsc, u_petsc)` to use the streamfunction formulation.
    If `u_petsc=None` then velocity is obtained by interpolation rather than projection.
    
    `flow_petsc = up_petsc` to use the velocity-pressure formulation.
    """
    return isinstance(flow_petsc, tuple)


def use_continuous_galerkin(
    stabilization: str | float | tuple[float, float] | None,
) -> bool:
    return isinstance(stabilization, (str, type(None)))
