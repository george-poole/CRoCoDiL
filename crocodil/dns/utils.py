import os
from typing import Callable

import numpy as np
from mpi4py import MPI
from dolfinx.mesh import Mesh

from lucifex.mesh import MeshBoundary, mesh_boundary, rectangle_mesh


DEFAULT_JIT_DIR = os.path.abspath(
    os.path.join(
        __file__,
        '../../..',
        '__jit__',
    )
)


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
