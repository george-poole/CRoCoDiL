import numpy as np
from ufl.core.expr import Expr

from lucifex.fem import Function, Constant
from lucifex.utils.fenicsx_utils import mesh_integral, mesh_axes
from lucifex.utils.array_utils import as_index
from lucifex.pde.advection_diffusion import flux


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
    y_index = as_index(
        y_axis, 
        y_target, 
        condition=lambda aprx, trgt: aprx <= trgt and np.abs(aprx - trgt) < tol, 
        msg='Mesh resolution must be chosen such that `y_target` is aligned with cell facets',
    )
    y_approx = y_axis[y_index]
    y_approx_plus = y_axis[y_index + 1]
    y_approx_minus = y_axis[y_index - 1]
    return (1 / Lx) * flux(
        'dS', 
        lambda x: x[1] - y_approx, 
        lambda x: x[1] - y_approx_plus, 
        lambda x: x[1] - y_approx_minus, 
        facet_side="+",
    )(u, a, d)


@mesh_integral
def mass_capillary(
    s: Function, 
    varphi: Function,
    epsilon: Constant | float,
) -> Expr:
    """
    `mᶜ = ∫ 𝜑s / ε dx` 
    """
    return varphi * s / epsilon


@mesh_integral
def mass_dissolved(
    c: Function, 
    phi: Function,
) -> Expr:
    """
    `mᴰ = ∫ ϕc dx` where `ϕ = 𝜑(1 - s)`
    """
    return phi * c

