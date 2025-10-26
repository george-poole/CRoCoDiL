from typing import Callable

import numpy as np
from dolfinx.mesh import Mesh
from dolfinx.fem import Function, Constant
from ufl.core.expr import Expr

from lucifex.mesh import MeshBoundary, mesh_boundary, rectangle_mesh


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
    

def mass_capillary_trapped(
    s: Function, 
    epsilon: Constant | float,
) -> Expr:
    """
    Mass capillary-trapped per unit per unit area (2D) or volume (3D)
    
    `ρᶜ = s / ε`

    for the mass integral

    `mᶜ = ∫ ρᶜ dx` 
    """
    return s / epsilon


def mass_dissolved(
    c: Function, 
    s: Function,
) -> Expr:
    """
    Mass dissolved per unit per unit area (2D) or volume (3D)

    `ρᴰ = ∫ (1 - s)·c dx` 

    for the mass integral

    `mᴰ = ∫ ρᴰ dx` 
    """
    return (1 - s) * c


def rectangle_domain(
    Lx: float,
    Ly: float,
    Nx: int,
    Ny: int,
    cell: str,
    name: str = 'LxLy',
    clockwise_names: tuple[str, str, str, str] = ('upper', 'right', 'lower', 'left'),
) -> tuple[Mesh, MeshBoundary]:
    
    mesh = rectangle_mesh(Lx, Ly, Nx, Ny, name, cell)
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
