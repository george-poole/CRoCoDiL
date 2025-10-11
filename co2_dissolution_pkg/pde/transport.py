from typing import Callable

import numpy as np
from dolfinx.fem import Function, Constant
from ufl import (dx, dS, Form, CellDiameter, FacetNormal,
                 div, TestFunction,
                 jump, avg, dot, conditional, gt)
from ufl.geometry import CellDiameter
from ufl.core.expr import Expr

from lucifex.solver import IBVP, ibvp_solver, OptionsPETSc, BoundaryConditions
from lucifex.fdm import (DT, AB1, FiniteDifference, FunctionSeries, ConstantSeries, Series, 
                        apply_finite_difference)
from lucifex.fdm.ufl_operators import inner, grad
from lucifex.utils import is_tensor, extract_integrand

from .stabilization import supg_diffusivity, supg_reaction, supg_velocity, supg_tau
from .utils import ExplicitDiscretizationError


def advection_diffusion_cg(
    c: FunctionSeries,
    dt: Constant,
    phi: Series | Function | Expr,
    u: FunctionSeries,
    Ra: Constant,
    d: Series | Function | Expr,
    D_adv: FiniteDifference | tuple[FiniteDifference, FiniteDifference],
    D_diff: FiniteDifference,
    D_phi: FiniteDifference = AB1,
    supg: str | None = None,
    bcs: BoundaryConditions | None = None,
    expand: bool = False,
) -> list[Form]:
    """
    `ϕ∂c/∂t + 𝐮·∇c = 1/Ra ∇·(D·∇c)`
    """
    v = TestFunction(c.function_space)

    if isinstance(phi, Series):
        phi = D_phi(phi)
    if isinstance(d, Series):
        d = D_phi(d)

    dcdt = DT(c, dt)
    F_dcdt = v * dcdt * dx

    match D_adv:
        case D_adv_u, D_adv_c:
            adv = (1 / phi) * inner(D_adv_u(u, False), grad(D_adv_c(c)))
        case D_adv:
            adv = (1 / phi) * D_adv(inner(u, grad(c)))
    # NOTE equivalent to 
    # adv = (1 / phi) * apply_finite_difference(
    #     D_adv,
    #     (lambda u, c: inner(u, grad(c)), (u, c)),
    #     c,
    # )
    F_adv = v * adv * dx

    diff = -(1/Ra) * (1/phi) * div(d * grad(D_diff(c)))
    if expand:
        assert not is_tensor(d)
        diff_expand = (
            -(1/Ra) * (d / phi) * div(grad(D_diff(c))), 
            -(1/Ra) * (1/phi) * inner(grad(d), grad(D_diff(c))), 
        )
        F_diff = (1/Ra) * inner(grad(v), grad(D_diff(c))) * dx + v * diff_expand[1] * dx
    else:
        F_diff = (1/Ra) * inner(grad(v / phi), d * grad(D_diff(c))) * dx

    forms = [F_dcdt, F_adv, F_diff]

    if supg is not None:
        u_eff = supg_velocity(phi, u, Ra, D_adv, D_diff)
        d_eff = supg_diffusivity(Ra, D_diff)
        tau = supg_tau(supg, c.function_space.mesh, u_eff, d_eff)        
        res = dcdt + adv + diff
        F_res = tau * inner(grad(v), u_eff) * res * dx
        forms.append(F_res)

    if bcs is not None:
        ds, c_neumann = bcs.boundary_data(c.function_space, 'neumann')
        F_neumann = sum([-(1 / Ra) * v * cN * ds(i) for i, cN in c_neumann])
        forms.append(F_neumann)

    return forms


# TODO debug and test
def advection_diffusion_dg(
    c: FunctionSeries,
    dt: Constant,
    phi: Function | Expr | Series,
    u: FunctionSeries,
    Ra: Constant,
    d: Function | Expr | Series, 
    alpha: float,
    gamma: float,
    D_adv: tuple[FiniteDifference, FiniteDifference],
    D_diff: FiniteDifference,
    D_phi: FiniteDifference = AB1,
    bcs: BoundaryConditions | None = None,
) -> list[Form]:
    if bcs is None:
        bcs = BoundaryConditions()
    ds, c_dirichlet, c_neumann = bcs.boundary_data(c.function_space, 'dirichlet', 'neumann')

    if isinstance(phi, Series):
        phi = D_phi(phi)
    if isinstance(d, Series):
        d = D_phi(d)

    v = TestFunction(c.function_space)
    h = CellDiameter(c.function_space.mesh)
    n = FacetNormal(c.function_space.mesh)

    F_dcdt = v * DT(c, dt) * dx

    D_adv_u, D_adv_c = D_adv
    uEff = (1/phi) * (D_adv_u(u) - (1/Ra) * grad(phi))
    cAdv = D_adv_c(c)
    outflow = conditional(gt(dot(uEff, n), 0), 1, 0)

    F_adv_dx = -inner(grad(v), uEff * cAdv) * dx
    F_adv_dS = 2 * inner(jump(v, n), avg(outflow * uEff * cAdv)) * dS
    F_adv_ds = inner(v, outflow * inner(uEff, n) * cAdv) * ds 
    F_adv_ds += sum([inner(v, (1 - outflow) * inner(uEff, n) * cD) * ds(i) for cD, i in c_dirichlet])
    #### alternatice c.f. wells, burkardt
    # un = 0.5 * (inner(u, n) + abs(inner(u, n)))
    # un = outflow * inner(u, n)
    # F_adv_dS = inner(jump(v), un('+') * c('+') - un('-') * c('-')) * dS
    # F_adv_dS = inner(jump(v), jump(un * c)) * dS
    # F_adv_ds = v * (un * c + ud * cD) * ds(i)   # NOTE this applies DiricletBC on inflow only?
    ####
    F_adv = F_adv_dx + F_adv_dS + F_adv_ds

    cDiff = D_diff(c)
    # + ∫ ∇v⋅∇c dx
    F_diff_dx = inner(grad(v), grad(cDiff)) * dx
    # - ∫ [vn]⋅{∇c} dS
    F_diff_dS = -inner(jump(v, n), avg(grad(cDiff))) * dS
    # - ∫ [vn]⋅{∇c} dS
    F_diff_dS += -inner(avg(grad(v)), jump(cDiff, n)) * dS
    # + ∫ (α / h)[vn]⋅[cn] dS
    F_diff_dS += (alpha / avg(h)) * inner(jump(v, n), jump(cDiff, n)) * dS # TODO h('+') or avg(h) ?
    # ...
    F_diff_ds = sum([-(inner(grad(v), (cDiff - cD) * n) + inner(v * n, grad(cDiff))) * ds(i) for cD, i in c_dirichlet])
    F_diff_ds += sum([(gamma / h) * v * (cDiff - cD) * ds(i) for cD, i in c_dirichlet])
    F_diff_ds += sum([-v * cN * ds(i) for cN, i in c_neumann])
    F_diff = (1/Ra) * (F_diff_dx + F_diff_dS + F_diff_ds)

    return [F_dcdt, F_adv, F_diff]


def advection_diffusion_reaction_cg(
    c: FunctionSeries,
    dt: Constant,
    phi: Function | Expr | Series,
    u: FunctionSeries,
    Ra: Constant,
    d: Function | Expr | Series, 
    Da: Constant,
    r: Function | Expr | Series | tuple[Callable, tuple],
    D_adv: FiniteDifference | tuple[FiniteDifference, FiniteDifference],
    D_diff: FiniteDifference,
    D_reac: FiniteDifference | tuple[FiniteDifference, ...],
    D_phi: FiniteDifference = AB1,
    supg: str | None = None,
    bcs: BoundaryConditions | None = None,
    expand: bool = False,
) -> list[Form]:
    """
    `ϕ∂c/∂t + 𝐮·∇c = 1/Ra ∇·(D·∇c) + Da R`
    """
    if np.isclose(float(Da), 0):
        return advection_diffusion_cg(c, dt, phi, u, Ra, d, D_adv, D_diff, D_phi, supg, bcs, expand)
    else:
        forms = advection_diffusion_cg(c, dt, phi, u, Ra, d, D_adv, D_diff, D_phi, None, bcs, expand)

    if isinstance(phi, Series):
        phi = D_phi(phi)

    v = TestFunction(c.function_space)
    r = apply_finite_difference(D_reac, r, c)
    reac = -Da * r / phi 
    F_reac = v * reac * dx

    forms.append(F_reac)

    if supg is not None:
        u_eff = supg_velocity(phi, u, Ra, D_adv, D_diff)
        d_eff = supg_diffusivity(Ra, D_diff)
        r_eff = supg_reaction(dt, phi, Da, D_reac)
        tau = supg_tau(supg, c.function_space.mesh, u_eff, d_eff, r_eff)   
        dcdt, adv, diff = (extract_integrand(f) for f in forms[:3])
        res = dcdt + adv + diff + reac
        F_res = tau * inner(grad(v), u_eff) * res * dx
        forms.append(F_res)

    return forms


def advection_diffusion_reaction_dg(
    c: FunctionSeries,
    dt: Constant,
    phi: Function | Expr | Series,
    u,
    Ra: Constant,
    d,
    Da: Constant,
    r,
    alpha: float,
    gamma: float,
    D_adv: tuple[FiniteDifference, FiniteDifference],
    D_diff: FiniteDifference,
    D_reac: FiniteDifference | tuple[FiniteDifference, FiniteDifference],
    D_phi: FiniteDifference = AB1,
    bcs: BoundaryConditions | None = None,
) -> list[Form]:
    
    forms = advection_diffusion_dg(c, dt, phi, u, Ra, d, D_adv, D_diff, D_phi, alpha, gamma, bcs)

    if np.isclose(float(Da), 0):
        return forms
    
    if isinstance(phi, Series):
        phi = D_phi(phi)

    r = apply_finite_difference(D_reac, r, c)
    v = TestFunction(c.function_space)
    F_reac = -v * Da * (r / phi) * dx
    forms.append(F_reac)

    return forms


def advection_diffusion_solver(
    c,
    dt,
    phi,
    u,
    Ra,
    d,
    D_adv,
    D_diff,
    bcs,
    petsc,
    limits,
    stabilization,
) -> IBVP:
    petsc = OptionsPETSc("gmres", "none") if petsc is None else petsc
    if use_cts(stabilization):
        concentration_solver = ibvp_solver(advection_diffusion_cg, bcs=bcs, petsc=petsc, dofs_corrector=limits)
        return concentration_solver(
            c, dt, phi, u, Ra, d, D_adv, D_diff, AB1, stabilization)
    else:
        concentration_solver = ibvp_solver(advection_diffusion_dg, petsc=petsc)
        alpha_dg, gamma_dg = dg_penalty_parameters(stabilization)
        return concentration_solver(
            c, dt, phi, u, Ra, d, D_adv, D_diff, AB1, alpha_dg, gamma_dg, bcs)


def advection_diffusion_reaction_solver(
    c,
    dt,
    phi,
    u,
    Ra,
    d,
    Da, 
    r,
    D_adv,
    D_diff,
    D_reac,
    bcs,
    petsc,
    limits,
    stabilization,
):
    petsc = OptionsPETSc("gmres", "none") if petsc is None else petsc
    if use_cts(stabilization):
        concentration_solver = ibvp_solver(advection_diffusion_reaction_cg, bcs=bcs, petsc=petsc, dofs_corrector=limits)
        return concentration_solver(
            c, dt, phi, u, Ra, d, Da, r, D_adv, D_diff, D_reac, AB1, stabilization)
    else:
        concentration_solver = ibvp_solver(advection_diffusion_reaction_dg, petsc=petsc, dofs_corrector=limits)
        alpha_dg, gamma_dg = dg_penalty_parameters(stabilization)
        return concentration_solver(
            c, dt, phi, u, Ra, d, Da, r, alpha_dg, gamma_dg, D_adv, D_diff, D_reac, AB1, bcs)
    

def dg_penalty_parameters(
    stabilization: str | tuple[float, float] | None,
) -> tuple[float, float]:
    """
    `α` and `γ`
    """
    if stabilization == 'dg':
        alpha, gamma = 10.0, 10.0
    else:
        alpha, gamma = stabilization
    return alpha, gamma


def use_cts(
    stabilization: str | tuple[float, float] | None,
) -> bool:
    if stabilization is None:
        return True
    if isinstance(stabilization, str) and stabilization != 'dg':
        return True
    return False

    
