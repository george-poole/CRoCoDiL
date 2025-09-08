from dolfinx.fem import Function, Constant
from ufl import exp, as_tensor, cos, sin
from ufl.core.expr import Expr

from lucifex.fdm import FunctionSeries, ExprSeries


def rock_porosity_layered(
):
    """
    `𝜑(𝐱) = ...`
    """
    ...

def effective_porosity(
    varphi: Function | Constant | float,
    s: Function | FunctionSeries,
) -> Function | FunctionSeries:
    """
    `ϕ(𝐱,t) = 𝜑(𝐱)(1 - S(𝐱,t))`
    """
    return varphi * (1 - s)


def permeability_power_law(
    phi: Function,
    n: float,
) -> Expr:
    """
    `K(ϕ) = ϕⁿ`
    """
    return phi ** n


def permeability_cross_bedded(
    phi,
    kappa,
    vartheta,
    n,
):
    """
    `𝖪(ϕ) = K(ϕ) (
        (cos²ϑ + κsin²ϑ , (1 - κ)cosϑsinϑ), 
        ((1 - κ)cosϑsinϑ , κcos²ϑ + sin²ϑ), 
    )`
    """
    k = permeability_power_law(phi, n)
    cs = cos(vartheta)
    sn = sin(vartheta)  
    tensor = as_tensor(
        (
            (cs**2 + kappa*sn**2, (1 - kappa)*cs*sn),
            ((1 - kappa)*cs*sn, kappa*cs**2 + sn**2), 
        ),
    )
    return k * tensor


def permeability_exponential(
    phi,
    a: float,
):
    """
    `K(ϕ) = exp(aϕ)`
    """
    return exp(a * phi)
    


def density_power_law(
    c: Function,
    n: float,
) -> Expr:
    """
    `ρ(c) = cⁿ`
    """
    return c ** n


def viscosity_power_law(
    c: Function,
    beta: float,
    n: float,
) -> Expr:
    """
    `μ(c) = 1 + βcⁿ`
    """
    return 1 + beta * c ** n


def reaction_power_law(
    s: Function | FunctionSeries,
    c: Function,
    a: float,
    b: float,
    n: float,
    ce: float | Function | Expr,
) -> Expr | ExprSeries:
    """
    `r(s,c) = sᵃ(1 - s)ᵇ(ce - c)ⁿ`
    """
    return (s**a) * (1 - s)**b * (ce - c)**n