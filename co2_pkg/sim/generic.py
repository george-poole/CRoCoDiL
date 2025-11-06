from typing import Callable, Iterable, TypeAlias
from types import EllipsisType

import numpy as np
from ufl import as_vector
from ufl.core.expr import Expr
from dolfinx.mesh import Mesh

from lucifex.fem import Function, Constant
from lucifex.mesh import MeshBoundary
from lucifex.fdm import (
    FunctionSeries, ConstantSeries, FiniteDifference, AB1, Series, 
    ExprSeries, FiniteDifferenceArgwise, finite_difference_order,
    cflr_timestep, cfl_timestep, reactive_timestep,
)
from lucifex.solver import (
    BoundaryConditions, OptionsPETSc, bvp, ibvp, ivp, 
    interpolation, projection, evaluation, integration,
)
from lucifex.utils import CellType, extremum, div_norm
from lucifex.sim import Simulation
from lucifex.pde.streamfunction import streamfunction_velocity
from lucifex.pde.darcy import darcy_streamfunction, darcy_incompressible
from lucifex.pde.transport import advection_diffusion_reaction, advection_diffusion, flux
from lucifex.pde.transport_dg import advection_diffusion_reaction_dg, advection_diffusion_dg
from lucifex.pde.evolution import evolution_forms, evolution_expression

from .utils import mass_dissolved, mass_capillary_trapped


Phi: TypeAlias = FunctionSeries
C: TypeAlias = FunctionSeries
Theta: TypeAlias = FunctionSeries
S: TypeAlias = FunctionSeries
U: TypeAlias = FunctionSeries
def thermosolutal_transport_generic(
    # domain
    Omega: Mesh,
    dOmega: MeshBoundary,
    # gravity
    eg: tuple[Expr | Function | Constant | float, ...] | None = None,
    # physical 
    epsilon: float | None = None,
    # initial conditions
    c_ics = None,
    theta_ics = None,
    s_ics = None,
    # boundary conditions
    flow_bcs: BoundaryConditions | EllipsisType | None | tuple[BoundaryConditions | None, BoundaryConditions | None] = ...,
    c_bcs: BoundaryConditions | EllipsisType | None = None,
    theta_bcs: BoundaryConditions | EllipsisType | None = None,
    # constitutive relations
    rock_porosity: Callable[[np.ndarray], np.ndarray] | float = 1,
    permeability: Callable[[Phi], Series] = lambda phi: phi**2,
    dispersion_solutal: Callable[[Phi, U], Series] | None = lambda phi, _: phi,
    dispersion_thermal: Callable[[Phi, U], Series] | None = lambda phi, _: phi,
    density: Callable[[C | Theta], Series] 
    | Callable[[C, Theta], Series] = lambda c, theta: c - theta,
    viscosity: Callable[[C | Theta], Series]
    | Callable[[C, Theta], Series] = lambda *_: 1 + 0 * _[0], 
    reaction: Callable[[S, C], Series]
    | Callable[[S, C, Theta], Series] = lambda s, c, theta: s * (1 - c) * theta,
    # time step
    dt_min: float = 0.0,
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float | None = 0.75,
    k_courant: float | None = None,
    # time discretization
    D_adv_solutal: FiniteDifference | FiniteDifferenceArgwise = AB1,
    D_diff_solutal: FiniteDifference = AB1,
    D_reac_solutal: FiniteDifference 
    | FiniteDifferenceArgwise = AB1,
    D_adv_thermal: FiniteDifference | FiniteDifferenceArgwise = AB1,
    D_diff_thermal: FiniteDifference = AB1,
    D_reac_evol: FiniteDifference 
    | FiniteDifferenceArgwise = AB1,
    # stabilization
    c_stabilization: str | tuple[float, float] | None = None,
    c_limits: tuple[float, float] | EllipsisType | None = None,
    theta_stabilization: str | tuple[float, float] = None,
    theta_limits: tuple[float, float] | EllipsisType | None = None,
    s_limits: tuple[float, float] | EllipsisType | None = None,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc | None, OptionsPETSc | EllipsisType | None] | OptionsPETSc | None = (None, ...),
    c_petsc: OptionsPETSc | None = None,
    theta_petsc: OptionsPETSc | None = None,
    s_petsc: OptionsPETSc | EllipsisType | None = ...,
    # optional solvers
    secondary: bool = False,    
    # solvers
    secondary_extras: Iterable = (),
    namespace_extras: Iterable = (),
) -> Simulation:    
    """
    `ϕ∂c/∂t + 𝐮·∇c =  ∇·(D(ϕ,𝐮)·∇c) + R(s,c,θ)` \\
    `ϕ∂θ/∂t + 𝐮·∇θ = ∇·(G(ϕ,𝐮)·∇θ)`\\
    `∇⋅𝐮 = 0` \\
    `𝐮 = -(∇p + ρ(c,θ)e₉)` \\
    `𝜑∂s/∂t = -εR(s,c,θ)`

    `ϕ = 𝜑(1 - s)` is the effective porosity.
    
    Default boundary conditions are no flux of fluid, solute and heat everywhere on `∂Ω`. 

    Default gravity unit vector is `e₉ = -eʸ` in 2D or `e₉ = -eᶻ` in 3D.
    
    Default constitutive relations are uniform rock porosity `𝜑 = 1`, 
    isotropic quadratic permeability `K(ϕ) = ϕ²`, isotropic linear solutal
    dispersion `D(ϕ) = ϕ`, isotropic linear thermal dispersion `G(ϕ) = ϕ`, 
    uniform viscosity `μ = 1`.

    In general: rock porosity is a spatial function `𝜑(𝐱)`; permeability `K(ϕ)` is a 
    function of porosity; solutal and thermal dispersions are functions `D(ϕ, 𝐮)`, `G(ϕ, 𝐮)` of porosity
    and velocity; density and viscosity are functions `μ(c, θ)`, `ρ(c, θ)` of concentration
    and temperature; reaction rate is a function `R(s, c, θ)` of saturation, concentration
    and temperature.
    """
    if eg is None:
        if Omega.geometry.dim == 2:
            eg = (0, -1)
        if Omega.geometry.dim == 3:
            eg = (0, 0, -1)

    STREAMF = isinstance(flow_petsc, tuple)
    CG = lambda stab: stab is None
    SOLUTAL = dispersion_solutal is not None
    THERMAL = dispersion_thermal is not None
    THERMOSOLUTAL = SOLUTAL and THERMAL
    EVOL = epsilon is not None

    solvers = []
    namespace = []

    # time
    order = finite_difference_order(D_adv_solutal, D_diff_solutal, D_reac_solutal, D_adv_thermal, D_diff_thermal, D_reac_evol)
    t = ConstantSeries(Omega, "t", order, ics=0.0)  
    dt = ConstantSeries(Omega, 'dt')

    # flow
    if STREAMF: 
        psi_deg = 2
        psi = FunctionSeries((Omega, 'P', psi_deg), 'psi')
        u = FunctionSeries((Omega, "P", psi_deg - 1, 2), "u", order)
    else:
        u_fam = 'BDMCF' if Omega.topology.cell_name() == CellType.QUADRILATERAL else 'BDM'
        u_deg = 1
        up = FunctionSeries((Omega, [(u_fam, u_deg), ('DP', u_deg - 1)]), ('up', ['u', 'p']), order)
        u = up.sub(0)

    # transport
    if SOLUTAL:
        c = FunctionSeries((Omega, 'P' if CG(c_stabilization) else 'DP', 1), 'c', order, ics=c_ics)
    if THERMAL:
        theta = FunctionSeries((Omega, 'P' if CG(theta_stabilization) else 'DP', 1), 'theta', order, ics=theta_ics)
    if EVOL:
        epsilon = Constant(Omega, epsilon, 'epsilon')
        s = FunctionSeries((Omega, 'P', 1), 's', order, ics=s_ics)
    else:
        s = Function((Omega, 'P', 1), 's', s_ics)

    # constitutive
    varphi = Function((Omega, 'P', 1), rock_porosity, 'varphi')
    phi = ExprSeries(varphi * (1 - s), 'phi')
    k = ExprSeries(permeability(phi), 'k')
    rho = ExprSeries(
        density(c, theta) if THERMOSOLUTAL
        else density(c) if SOLUTAL
        else density(theta),
        'rho',
    )
    mu = ExprSeries(
        viscosity(c, theta) if THERMOSOLUTAL
        else viscosity(c) if SOLUTAL
        else viscosity(theta),
        'mu',
    )
    namespace.extend((epsilon, varphi, phi, k, rho, mu))
    if SOLUTAL:
        d = ExprSeries(dispersion_solutal(phi, u), 'd')
        r = ExprSeries.from_args(reaction, name='r')(
            *(s, c, theta) if THERMOSOLUTAL
            else (s, c)
        )
        namespace.extend((d, r))
    if THERMAL:
        g = ExprSeries(dispersion_thermal(phi, u), 'g')
        namespace.append(g)

    # flow solvers
    if STREAMF:
        psi_bcs = BoundaryConditions(("dirichlet", dOmega.union, 0.0)) if flow_bcs is Ellipsis else flow_bcs
        psi_petsc, u_petsc = flow_petsc
        psi_petsc = OptionsPETSc("gmres", "none") if psi_petsc is None else psi_petsc
        egx, egy = eg
        psi_solver = bvp(darcy_streamfunction, psi_bcs, petsc=psi_petsc)(
            psi, k[0], mu[0], egx * rho[0], egy * rho[0],
        ) 
        if u_petsc is Ellipsis:
            u_solver = interpolation(u, streamfunction_velocity)(psi[0])
        else:
            u_petsc = OptionsPETSc("gmres", "none") if u_petsc is None else u_petsc
            u_solver = projection(u, streamfunction_velocity, petsc=u_petsc)(psi[0])
        
        solvers.extend((psi_solver, u_solver))
    else:
        u_bcs = BoundaryConditions(('essential', dOmega.union, (0.0, 0.0), 0)) if flow_bcs is Ellipsis else flow_bcs[0]
        p_bcs = None if flow_bcs is Ellipsis else flow_bcs[1]
        flow_petsc = OptionsPETSc("gmres", "lu") if flow_petsc is None else flow_petsc
        flow_petsc['pc_factor_mat_solver_type'] = 'mumps'
        u_solver = bvp(darcy_incompressible, u_bcs, petsc=flow_petsc)(
            up, k[0], mu[0], rho[0] * as_vector(eg), p_bcs)
        solvers.append(u_solver)

    # timestep solver
    if not SOLUTAL:
        k_courant = None
    dt_solver = evaluation(dt, cflr_timestep)(
        u[0], r[0], cfl_h, cfl_courant, k_courant, dt_max, dt_min,
    ) 
    solvers.append(dt_solver)

    # thermal solver
    if THERMAL:
        theta_bcs = BoundaryConditions(("neumann", dOmega.union, 0.0)) if theta_bcs is Ellipsis else theta_bcs
        theta_limits = (0, 1) if theta_limits is Ellipsis else theta_limits
        theta_petsc = OptionsPETSc("gmres", "none") if theta_petsc is None else theta_petsc
        if CG(theta_stabilization):
            theta_solver = ibvp(advection_diffusion, bcs=theta_bcs, petsc=theta_petsc, dofs_corrector=theta_limits)(
                c, dt, u, d, D_adv_thermal, D_diff_thermal, phi=phi, supg=theta_stabilization,
            )
        else:
            theta_alpha, theta_gamma = theta_stabilization
            theta_solver = ibvp(advection_diffusion_dg, petsc=theta_petsc, dofs_corrector=theta_limits)(
                c, dt, u, d, theta_alpha, theta_gamma, D_adv_thermal, D_diff_thermal, phi=phi, bcs=theta_bcs,
            )
        solvers.append(theta_solver)

    # solutal solver
    if SOLUTAL:
        c_bcs = BoundaryConditions(("neumann", dOmega.union, 0.0)) if c_bcs is Ellipsis else c_bcs
        c_limits = (0, 1) if c_limits is Ellipsis else c_limits
        c_petsc = OptionsPETSc("gmres", "none") if c_petsc is None else c_petsc
        if CG(c_stabilization):
            c_solver = ibvp(advection_diffusion_reaction, bcs=c_bcs, petsc=c_petsc, dofs_corrector=c_limits)(
                c, dt, u, d, r, D_adv_solutal, D_diff_solutal, D_reac_solutal, phi=phi, supg=c_stabilization,
            )
        else:
            c_alpha, c_gamma = c_stabilization
            c_solver = ibvp(advection_diffusion_reaction_dg, petsc=c_petsc, dofs_corrector=c_limits)(
                c, dt, u, d, r, c_alpha, c_gamma, D_adv_solutal, D_diff_solutal, D_reac_solutal, phi=phi, bcs=c_bcs,
            )
        solvers.append(c_solver)

    # evolution solver
    if EVOL:
        reaction_evol = lambda *args: -epsilon * reaction(*args)
        r_evol = ExprSeries.from_args(reaction_evol, name='rEvol')(
            *(s, c, theta) if THERMOSOLUTAL
            else (s, c)
        )
        namespace.append(r_evol)
        s_limits = (0, max(s.ics.x)) if s_limits is Ellipsis else s_limits
        if s_petsc is Ellipsis:
            s_solver = interpolation(s, evolution_expression, dofs_corrector=s_limits, future=True)(
                s, dt, r_evol, D_reac_evol, phi=varphi,
            )
        else:
            s_petsc = OptionsPETSc("gmres", "none") if s_petsc is None else s_petsc
            s_solver = ivp(evolution_forms, petsc=s_petsc, dofs_corrector=s_limits)(
                s, dt, r_evol, D_reac_evol, phi=varphi
            )
        solvers.append(s_solver)

    # optional solvers
    if secondary:
        uMinMax = ConstantSeries(Omega, "uMinMax", shape=(2,))
        solvers.append(evaluation(uMinMax, extremum)(u[0]))
        uDiv = ConstantSeries(Omega, 'uDiv')
        solvers.append(integration(uDiv, div_norm, 'dx')(u[0], 2))
        dtCFL = ConstantSeries(Omega, "dtCFL")
        solvers.append(evaluation(dtCFL, cfl_timestep)(u[0], cfl_h))
        if SOLUTAL:
            cMinMax = ConstantSeries(Omega, "cMinMax", shape=(2,))
            solvers.append(evaluation(cMinMax, extremum)(c[0]))
            mD = ConstantSeries(Omega, "mD")
            solvers.append(integration(mD, mass_dissolved, 'dx')(c[0], s[0]))
            fBoundary = ConstantSeries(Omega, "fBoundary", shape=(len(dOmega.union), 2))
            solvers.append(integration(fBoundary, flux, 'ds', *dOmega.union)(c[0], u[0], d[0]))
            dtK = ConstantSeries(Omega, "dtK")
            solvers.append(evaluation(dtK, reactive_timestep)(r[0]))
        if EVOL:
            sMinMax = ConstantSeries(Omega, "sMinMax", shape=(2,))
            solvers.append(evaluation(sMinMax, extremum)(s[0]))
            mC = ConstantSeries(Omega, "mC")
            solvers.append(integration(mC, mass_capillary_trapped, 'dx')(s[0], epsilon))
        if THERMAL:
            thetaMinMax = ConstantSeries(Omega, "thetaMinMax", shape=(2,))
            solvers.append(evaluation(thetaMinMax, extremum)(theta[0]))
            jBoundary = ConstantSeries(Omega, "jBoundary", shape=(len(dOmega.union), 2)),
            solvers.append(integration(jBoundary, flux, 'ds', *dOmega.union)(theta[0], u[0], g[0]))

    solvers.extend(secondary_extras)
    namespace.extend(namespace_extras)
    
    return Simulation(solvers, t, dt, namespace)