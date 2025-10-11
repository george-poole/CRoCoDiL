from typing import Callable

from dolfinx.fem import Function, Constant

from ufl import (dx, Form, FacetNormal,
                 as_matrix, Dx, div, TrialFunction, TestFunction,
                 det, transpose, TestFunctions, TrialFunctions, 
                 inv, as_vector, as_vector)
from ufl.core.expr import Expr

from lucifex.fdm import FunctionSeries
from lucifex.solver import (
    projection_solver, 
    interpolation_solver, bvp_solver,
    OptionsPETSc, BoundaryValueProblem, ProjectionProblem, InterpolationProblem)
from lucifex.solver import BoundaryConditions
from lucifex.fdm.ufl_operators import inner, grad
from lucifex.utils import is_tensor


def streamfunction_velocity(psi: Function) -> Expr:
    """
    `𝐮 = ((0, 1), (-1, 0))·∇ψ`
    """
    return as_matrix([[0, 1], [-1, 0]]) * grad(psi)


def darcy_streamfunction(
    psi: FunctionSeries,
    k: Expr | Function | Constant | float,
    mu: Expr | Function | Constant | float,
    rho: Expr | Function,
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
        F_egx = -v * Dx(egx * rho, 1) * dx
        forms.append(F_egx)
    if egy is not None:
        F_egy = v * Dx(egy * rho, 0) * dx
        forms.append(F_egy)
    return forms


def streamfunction_solvers(
    psi,
    u,
    k,
    mu,
    rho,
    egx,
    egy,
    psi_bcs = None,
    petsc: tuple = (None, ...),
) -> tuple[BoundaryValueProblem, ProjectionProblem | InterpolationProblem]:
    psi_petsc, u_petsc = petsc
    psi_petsc = OptionsPETSc("gmres", "none") if psi_petsc is None else psi_petsc
    psi_solver = bvp_solver(darcy_streamfunction, psi_bcs, petsc=psi_petsc)(
        psi, k, mu, rho, egx, egy,
    ) 
    if u_petsc is Ellipsis:
        u_solver = interpolation_solver(u, streamfunction_velocity)(psi[0])  # FIXME u_bcs
    else:
        u_petsc = OptionsPETSc("gmres", "none") if u_petsc is None else u_petsc
        u_solver = projection_solver(u, streamfunction_velocity, petsc=u_petsc)(psi[0])  # FIXME u_bcs, petsc=
        
    return psi_solver, u_solver


def darcy_incompressible(
    up: FunctionSeries,
    rho,
    k,
    mu,
    egx: Expr | Function | Constant | float,
    egy: Expr | Function | Constant | float,
    egz: Expr | Function | Constant | float | None = None,
    p_bcs: BoundaryConditions | None = None,
) -> list[Form]:
    """
    `𝐮 = -(K/μ)⋅(∇p + ρe₉)` \\
    `∇⋅𝐮 = 0`
    
    `F(𝐮,p;𝐯,q) = ∫ q(∇·𝐮) dx ` \\
    `+ ∫ 𝐯·(μ K⁻¹⋅𝐮) dx - ∫ p(∇·𝐯) dx - ∫ 𝐯·ρe₉ dx + ∫ p(𝐯·n) ds`
    """
    v, q = TestFunctions(up.function_space)
    u, p = TrialFunctions(up.function_space)
    n = FacetNormal(up.function_space.mesh)

    dim = up.function_space.mesh.geometry.dim
    if dim == 2:
        eg = as_vector([egx, egy])
    if dim == 3:
        assert egz is not None
        eg = as_vector([egx, egy, egz])

    if is_tensor(k):
        F_velocity = inner(v, mu * inv(k) * u) * dx
    else:
        F_velocity = inner(v, mu * u / k) * dx
    F_pressure = -p * div(v) * dx
    F_buoyancy = inner(v, eg) * rho * dx
    F_div = q * div(u) * dx

    forms = [F_velocity, F_pressure, F_buoyancy, F_div]

    if p_bcs is not None:
        ds, p_natural = p_bcs.boundary_data(up.function_space, 'natural')
        F_bcs = sum([inner(v, n) * pN * ds(i) for i, pN in p_natural])
        forms.append(F_bcs)

    return forms


def darcy_solver(
    up,
    rho,
    k,
    mu,
    egx,
    egy,
    egz,
    u_bcs,
    p_bcs,
    petsc: OptionsPETSc | None = None,
):
    petsc = OptionsPETSc("gmres", "lu") if petsc is None else petsc
    petsc['pc_factor_mat_solver_type'] = 'mumps'
    return bvp_solver(darcy_incompressible, u_bcs, petsc=petsc)(up, rho, k, mu, egx, egy, egz, p_bcs)


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
    `𝐮 = -(K/μ)⋅(∇p + ρĝ)` \\
    `∇⋅𝐮 = -ε(1 - η)Da R(s,c)` 
    """
    forms = darcy_incompressible(u_p, rho, k, mu, egx, egy, egz, p_bcs)
    q = TestFunctions(u_p.function_space)[1]
    F_reac = q * epsilon * (1 - eta) * Da * reaction(s, c) * dx
    forms.append(F_reac)
    return forms


def use_streamf(
    petsc: OptionsPETSc | tuple,
) -> bool:
    if isinstance(petsc, tuple):
        return True
    return False
