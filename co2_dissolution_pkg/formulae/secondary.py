from typing import Iterable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from dolfinx.fem import Function, Constant
from ufl import FacetNormal
from ufl.core.expr import Expr

from lucifex.fdm import inner, grad


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


def advective_flux(
    c: Function,
    u: Function | Constant,
) -> Expr:
    """
    `fᵁ = 𝐮·∇c`
    """
    return inner(u, grad(c))


def diffusive_flux(
    c: Function,
    d: Function,
    Ra: Constant
) -> Expr:
    """
    `fᴰ = 1/Ra 𝐧·(D·∇c)`
    """
    n = FacetNormal(c.function_space.mesh)
    return (1/Ra) * inner(n, d * grad(c))


def flux(
    c: Function,
    u: Function | Constant, 
    d: Function,
    Ra: Constant
) -> tuple[Expr, Expr]:
    """
    Advective flux 
    `fᵁ = 𝐮·∇c`

    and diffusive flux
    `fᴰ = 1/Ra 𝐧·(D·∇c)`

    per unit length (2D) or area (3D) for the flux integrals

    `Fᵁ = ∫ fᵁ ds` \\
    `Fᴰ = ∫ fᴰ ds`
    """
    return advective_flux(c, u), diffusive_flux(c, d, Ra)


# def interfacial_flux(
#     mC: Iterable[float],
#     t: Iterable[float],
#     s: Function  | np.ndarray, 
#     s_contour: float = 0.0,
# ) -> np.ndarray:
#     """
#     F = -1/|Γ| dmᶜ/dt 
#     """
#     return -np.gradient(mC, t) / contour_length(s, level=s_contour)

