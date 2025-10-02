from dolfinx.fem import Function, Constant
from ufl import exp, as_tensor, cos, sin
from ufl.core.expr import Expr

from lucifex.fdm import FunctionSeries, ExprSeries


def effective_porosity(
    varphi: Function | Constant | float,
    s: Function | FunctionSeries,
) -> Function | FunctionSeries:
    """
    `ϕ(𝐱,t) = 𝜑(𝐱)(1 - S(𝐱,t))`
    """
    return varphi * (1 - s)


def permeability_cross_bedded(
    Kphi,
    kappa,
    vartheta,
):
    """
    `𝖪(ϕ) = K(ϕ) (
        (cos²ϑ + κsin²ϑ , (1 - κ)cosϑsinϑ), 
        ((1 - κ)cosϑsinϑ , κcos²ϑ + sin²ϑ), 
    )`
    """
    cs = cos(vartheta)
    sn = sin(vartheta)  
    tensor = as_tensor(
        (
            (cs**2 + kappa*sn**2, (1 - kappa)*cs*sn),
            ((1 - kappa)*cs*sn, kappa*cs**2 + sn**2), 
        ),
    )
    return Kphi * tensor