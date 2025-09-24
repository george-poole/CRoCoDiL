from lucifex.fdm import AB1
from lucifex.solver import (
    projection_solver, 
    interpolation_solver, bvp_solver, ivp_solver, ibvp_solver, 
    OptionsPETSc, BoundaryValueProblem, ProjectionProblem, InterpolationProblem)


from .pde import (darcy_streamfunction, streamfunction_velocity, darcy_incompressible, 
                  advection_diffusion_cg, advection_diffusion_dg, 
                  advection_diffusion_reaction_dg, advection_diffusion_reaction_cg,
                  evolution_expression, evolution)


def streamfunction_solvers(
    psi,
    u,
    rho,
    k,
    mu,
    beta,
    psi_bcs = None,
    petsc: tuple = (None, ...),
) -> tuple[BoundaryValueProblem, ProjectionProblem | InterpolationProblem]:
    psi_petsc, u_petsc = petsc
    psi_petsc = OptionsPETSc("gmres", "none") if psi_petsc is None else psi_petsc
    psi_solver = bvp_solver(darcy_streamfunction, psi_bcs, petsc=psi_petsc)(
        psi, rho, k, mu, beta,
    ) 
    if u_petsc is Ellipsis:
        u_solver = interpolation_solver(u, streamfunction_velocity)(psi[0])  # FIXME u_bcs
    else:
        u_petsc = OptionsPETSc("gmres", "none") if u_petsc is None else u_petsc
        u_solver = projection_solver(u, streamfunction_velocity, petsc=u_petsc)(psi[0])  # FIXME u_bcs, petsc=
        
    return psi_solver, u_solver


def darcy_solver(
    up,
    rho,
    k,
    mu,
    beta,
    u_bcs,
    p_bcs,
    petsc: OptionsPETSc | None = None,
):
    petsc = OptionsPETSc("gmres", "none") if petsc is None else petsc
    petsc['pc_factor_mat_solver_type'] = 'mumps'
    return bvp_solver(darcy_incompressible, u_bcs, petsc=petsc)(up, rho, k, mu, beta, p_bcs)


def streamfunction_method(
    petsc: OptionsPETSc | tuple,
) -> bool:
    if isinstance(petsc, tuple):
        return True
    return False


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
):
    petsc = OptionsPETSc("gmres", "none") if petsc is None else petsc
    if continuous(stabilization):
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
    if continuous(stabilization):
        concentration_solver = ibvp_solver(advection_diffusion_reaction_cg, bcs=bcs, petsc=petsc, dofs_corrector=limits)
        return concentration_solver(
            c, dt, phi, u, Ra, d, Da, r, D_adv, D_diff, D_reac, AB1, stabilization)
    else:
        concentration_solver = ibvp_solver(advection_diffusion_reaction_dg, petsc=petsc, dofs_corrector=limits)
        alpha_dg, gamma_dg = dg_penalty_parameters(stabilization)
        return concentration_solver(
            c, dt, phi, u, Ra, d, Da, r, alpha_dg, gamma_dg, D_adv, D_diff, D_reac, AB1, bcs)
    

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
        return ivp_solver(evolution, petsc=s_petsc, dofs_corrector=s_limits)(
            s, dt, varphi, epsilon, Da, r, D_reac,
        )


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


def continuous(
    stabilization: str | tuple[float, float] | None,
) -> bool:
    if stabilization is None:
        return True
    if isinstance(stabilization, str) and stabilization != 'dg':
        return True
    return False
