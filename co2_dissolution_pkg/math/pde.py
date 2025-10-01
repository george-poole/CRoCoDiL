from typing import Callable

from dolfinx.fem import Function, Constant
from ufl import (dx, dS, Form, CellDiameter, FacetNormal,
                 as_matrix, Dx, div, sin, cos, TrialFunction, TestFunction,
                 jump, avg, det, transpose, TestFunctions, TrialFunctions, 
                 inv, as_vector, dot, conditional, gt, as_vector)
from ufl.geometry import CellDiameter
from ufl.core.expr import Expr

from lucifex.solver import BoundaryConditions
from lucifex.fdm import (DT, AB1, FiniteDifference, FunctionSeries, ConstantSeries, Series, 
                        finite_difference_argwise)
from lucifex.fdm.ufl_operators import inner, grad
from lucifex.utils import is_tensor, extract_integrand

from .stabilization import supg_diffusivity, supg_reaction, supg_velocity, supg_tau


def darcy_streamfunction(
    psi: FunctionSeries,
    rho: Expr | Function,
    k: Expr | Function | Constant | float,
    mu: Expr | Function | Constant | float,
    egx: Expr | Function | Constant | float | None = None,
    egy: Expr | Function | Constant | float | None = None,
) -> list[Form]:
    """
    `∇·(μKᵀ·∇ψ / det(K)) = ∂(ρ·e₉ʸ)/∂x - ∂(ρ·e₉ˣ)/∂y`

    for tensor-valued `K` or

    `∇·(μ·∇ψ / K) = ∂(ρ·e₉ʸ)/∂x - ∂(ρ·e₉ˣ)/∂y`

    for scalar-valued `K`.
    """
    v = TestFunction(psi.function_space)
    psi_trial = TrialFunction(psi.function_space)
    if is_tensor(k):
        F_lhs = -(mu / det(k)) * inner(grad(v), transpose(k) * grad(psi_trial)) * dx 
    else:
        F_lhs = -(mu / k) * inner(grad(v), grad(psi_trial)) * dx
    forms = [F_lhs]
    if egx is not None:
        F_egx = v * Dx(egx * rho, 1) * dx
        forms.append(F_egx)
    if egy is not None:
        F_egy = -v * Dx(egy * rho, 0) * dx
        forms.append(F_egy)
    return forms


def streamfunction_velocity(psi: Function) -> Expr:
    """
    `𝐮 = ((0, 1), (-1, 0))·∇ψ`
    """
    return as_matrix([[0, 1], [-1, 0]]) * grad(psi)


def darcy_incompressible(
    u_p: FunctionSeries,
    rho,
    k,
    mu,
    egx: Expr | Function | Constant | float,
    egy: Expr | Function | Constant | float,
    egz: Expr | Function | Constant | float | None = None,
    p_bcs: BoundaryConditions | None = None,
) -> list[Form]:
    """
    `∇⋅𝐮 = 0` \\
    `𝐮 = -(K/μ)⋅(∇p + ρĝ)`
    
    `F(𝐮,p;𝐯,q) = ∫ q(∇·𝐮) dx ` \\
    `+ ∫ 𝐯·(μ K⁻¹⋅𝐮) dx - ∫ p(∇·𝐯) dx - ∫ 𝐯·ρĝ dx + ∫ p(𝐯·n) ds`
    """
    v, q = TestFunctions(u_p.function_space)
    u, p = TrialFunctions(u_p.function_space)
    n = FacetNormal(u_p.function_space.mesh)

    dim = u_p.function_space.mesh.geometry.dim
    if dim == 2:
        g = as_vector([egx, egy])
    if dim == 3:
        assert egz is not None
        g = as_vector([egx, egy, egz])

    if is_tensor(k):
        F_velocity = inner(v, mu * inv(k) * u) * dx
    else:
        F_velocity = inner(v, mu * u / k) * dx
    F_pressure = -p * div(v) * dx
    F_buoyancy = -inner(v, g) * rho * dx
    F_div = q * div(u) * dx

    forms = [F_velocity, F_pressure, F_buoyancy, F_div]

    if p_bcs is not None:
        ds, p_natural = p_bcs.boundary_data(u_p.function_space, 'natural')
        F_bcs = sum([inner(v, n) * pN * ds(i) for i, pN in p_natural])
        forms.append(F_bcs)

    return forms


def darcy_compressible(
    u_p: FunctionSeries,
    rho,
    k,
    mu,
    c: Function,
    s: Function,
    reaction: Callable,
    epsilon: Constant,
    eta: Constant,
    Da: Constant,
    egx: Expr | Function | Constant | float,
    egy: Expr | Function | Constant | float,
    egz: Expr | Function | Constant | float | None = None,
    p_bcs: BoundaryConditions | None = None,
) -> list[Form]:
    """
    `∇⋅𝐮 = -ε(1 - η)Da R(s,c)` \\
    `𝐮 = -(K/μ)⋅(∇p + ρĝ)`
    """
    forms = darcy_incompressible(u_p, rho, k, mu, egx, egy, egz, p_bcs)
    q = TestFunctions(u_p.function_space)[1]
    F_reac = q * epsilon * (1 - eta) * Da * reaction(s, c) * dx
    forms.append(F_reac)
    return forms


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
            adv = (1 / phi) * inner(D_adv_u[False](u), grad(D_adv_c(c)))
        case D_adv:
            adv = (1 / phi) * D_adv(inner(u, grad(c)))
    # NOTE equivalent to 
    # adv = (1 / phi) * finite_difference_argwise(
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


# FIXME # TODO tensor valued d
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
    if not Da:
        return advection_diffusion_cg(c, dt, phi, u, Ra, d, D_adv, D_diff, D_phi, supg, bcs, expand)
    else:
        forms = advection_diffusion_cg(c, dt, phi, u, Ra, d, D_adv, D_diff, D_phi, None, bcs, expand)

    if isinstance(phi, Series):
        phi = D_phi(phi)

    v = TestFunction(c.function_space)
    r = finite_difference_argwise(D_reac, r, c)
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

    if not Da:
        return forms
    
    if isinstance(phi, Series):
        phi = D_phi(phi)

    r = finite_difference_argwise(D_reac, r, c)
    v = TestFunction(c.function_space)
    F_reac = -v * Da * (r / phi) * dx
    forms.append(F_reac)

    return forms


def evolution(
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
    r = finite_difference_argwise(D_reac, r, s)
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
        
    DiscretizationError =  RuntimeError('Expression requires reaction term to be explicit in saturation')
    if isinstance(D_reac, FiniteDifference):
        if D_reac.is_implicit:
            raise DiscretizationError
    else:
        if D_reac[0].is_implicit:
            raise DiscretizationError

    r = finite_difference_argwise(D_reac, r, s)
    return s[0] - dt * (epsilon * Da / varphi) * r