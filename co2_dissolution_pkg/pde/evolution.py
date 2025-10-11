from typing import Callable

from dolfinx.fem import Function, Constant
from ufl import dx, TestFunction, Form
from ufl.core.expr import Expr

from lucifex.fdm import (
    DT, FiniteDifference, FunctionSeries, ConstantSeries, Series, 
    apply_finite_difference,
)
from lucifex.solver import (
    interpolation_solver, ivp_solver, 
    OptionsPETSc,
)

from .utils import ExplicitDiscretizationError



def evolution_forms(
    s: FunctionSeries,
    dt: Constant,
    varphi: Function | Constant | float,
    epsilon: Constant,
    Da: Constant,
    r: Function | Expr | Series | tuple[Callable, tuple],
    D_reac: FiniteDifference | tuple[FiniteDifference, ...],
) -> tuple[Form, Form]:
    """
    `𝜑 ∂s/∂t = -ε Da R(s, c)`
    """
    v = TestFunction(s.function_space)

    F_dsdt = v * DT(s, dt) * dx
    r = apply_finite_difference(D_reac, r, s)
    F_reac = v * (epsilon * Da / varphi) * r * dx

    return F_dsdt, F_reac


def evolution_expression(
    s: FunctionSeries,
    dt: Constant | ConstantSeries,
    varphi: Function | Constant | float,
    epsilon: Constant,
    Da: Constant,
    r: Series | Expr | Function,
    D_reac: FiniteDifference | tuple[FiniteDifference, ...],
) -> Expr:
    """
    `𝜑 ∂s/∂t = -ε Da R`

    rearranged after finite difference discretization into the algebraic expression

    `s¹ = s⁰ - Δt ε Da 𝒟(R) / 𝜑`.

    under the assumption that 𝒟(R) is explicit in `s`.
    """
    if isinstance(dt, ConstantSeries):
        dt = dt[0]
        
    if isinstance(D_reac, FiniteDifference):
        if D_reac.is_implicit:
            raise ExplicitDiscretizationError(D_reac, 'Reaction must be explicit w.r.t. saturation')
    else:
        if D_reac[0].is_implicit:
            raise ExplicitDiscretizationError(D_reac[0], 'Reaction must be explicit w.r.t. saturation')

    r = apply_finite_difference(D_reac, r, s)
    return s[0] - dt * (epsilon * Da / varphi) * r


def evolution_solver(
   s,
   dt,
   varphi,
   epsilon,
   Da,
   r,
   D_reac,    
   s_limits,
   s_petsc,
):
    if s_petsc is Ellipsis:
        return interpolation_solver(s, evolution_expression, dofs_corrector=s_limits, future=True)(
            s, dt, varphi, epsilon, Da, r, D_reac,
        )
    else:
        s_petsc = OptionsPETSc("gmres", "none") if s_petsc is None else s_petsc
        return ivp_solver(evolution_forms, petsc=s_petsc, dofs_corrector=s_limits)(
            s, dt, varphi, epsilon, Da, r, D_reac,
        )