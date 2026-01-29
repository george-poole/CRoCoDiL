from typing import Callable

import numpy as np
from dolfinx.mesh import Mesh
from ufl.core.expr import Expr

from lucifex.fem import Function, Constant
from lucifex.mesh import MeshBoundary, mesh_boundary, rectangle_mesh
from lucifex.utils import mesh_integral, as_index, mesh_axes
from lucifex.pde.advection_diffusion import flux
from lucifex.pde.scaling import ScalingOptions


CONVECTION_REACTION_SCALINGS = ScalingOptions(
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


def critical_wavelength(
    scaling: str,
    Ra: float,
) -> float:
    """
    `λ = 90 / Ra`
    """
    # TODO effects of scaling? other versions?
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra)
    Di = scaling_map['Di']
    return 90.0 * Di


def critical_dx(
    scaling: str,
    Ra: float,
    n_per_cell: int,
):
    """
    `Δx ≤ λ / N`
    """
    return critical_wavelength(scaling, Ra) / n_per_cell


def critical_Nx(
    scaling: str,
    Ra: float,
    n_per_cell: int,
    Lx: float
):
    """
    `Nₓ ≥ n Lₓ / λ`
    """
    return np.ceil(Lx / critical_dx(scaling, Ra, n_per_cell))


def vertical_flux(
    u: Function,
    a: Function,
    d: Function,
    y_target: float | int,
    Lx: float,
    tol: float | None = 1e-6,
) -> np.ndarray:
    """
    Evaluates the vertical advective and diffusive fluxes per unit length
     
    `Fᵁ = 1/Lₓ ∫ (𝐧·𝐚)u ds` \\
    `Fᴰ = 1/Lₓ ∫ 𝐧·(-D·∇u) ds`

    at heights `y ≃ y₀, y₀⁺, y₀⁻`. Returns `np.ndarray` of shape `(3, 2)`.
    """
    mesh = u.function_space.mesh
    y_axis = mesh_axes(mesh)[1]
    h0_index = as_index(
        y_axis, 
        y_target, 
        lambda aprx, trgt: aprx <= trgt and np.abs(aprx - trgt) < tol, 
        'Mesh resolution must be chosen such that `h0` is aligned with cell facets.',
    )
    h0_approx = y_axis[h0_index]
    h0_plus = y_axis[h0_index + 1]
    h0_minus = y_axis[h0_index - 1]
    return (1 / Lx) * flux(
        'dS', 
        lambda x: x[1] - h0_approx, 
        lambda x: x[1] - h0_plus, 
        lambda x: x[1] - h0_minus, 
        facet_side="+",
    )(u, a, d)


@mesh_integral
def mass_capillary_trapped(
    s: Function, 
    epsilon: Constant | float,
) -> Expr:
    """
    `mᶜ = ∫ s / ε dx` 
    """
    return s / epsilon


@mesh_integral
def mass_dissolved(
    c: Function, 
    phi: Function,
) -> Expr:
    """
    `mᴰ = ∫ ϕc dx` 
    """
    return phi * c


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
) -> tuple[Mesh, MeshBoundary]:
    """
    `Ω ∪ ∂Ω`
    """
    mesh = rectangle_mesh(Lx, Ly, Nx, Ny, cell, name)
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
