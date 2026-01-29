from collections.abc import Iterable
from typing import Callable, TypeAlias
from types import EllipsisType

import numpy as np
from ufl import as_vector, inner, sqrt
from ufl.core.expr import Expr
from dolfinx.mesh import Mesh

from lucifex.fem import Function, Constant, SpatialPerturbation
from lucifex.mesh import MeshBoundary
from lucifex.fdm import (
    FunctionSeries, ConstantSeries, FiniteDifference, FE, 
    ExprSeries, FiniteDifferenceArgwise, finite_difference_order,
    cflr_timestep, cfl_timestep, reactive_timestep,
)
from lucifex.solver import (
    BoundaryConditions, OptionsPETSc, Solver,
    bvp, ibvp, ivp, 
    interpolation, projection, evaluation, integration,
    extrema, div_norm, L_norm,
)
from lucifex.sim import Simulation
from lucifex.pde.streamfunction import streamfunction_velocity
from lucifex.pde.darcy import darcy_streamfunction, darcy_incompressible
from lucifex.pde.advection_diffusion import (
    advection_diffusion_reaction, advection_diffusion, 
    advection_diffusion_reaction_dg, advection_diffusion_dg, flux,
)
from lucifex.utils.mesh_utils import CellType
from lucifex.utils.dofs_utils import limits_corrector
from lucifex.utils.py_utils import arity
from lucifex.pde.evolution import evolution, evolution_update

from .utils import mass_dissolved, mass_capillary_trapped, vertical_flux


Phi: TypeAlias = Expr | ExprSeries
C: TypeAlias = FunctionSeries
Theta: TypeAlias = FunctionSeries
S: TypeAlias = FunctionSeries
U: TypeAlias = FunctionSeries
def dns_generic(
    # domain
    Omega: Mesh,
    dOmega: MeshBoundary,
    # gravity
    eg: tuple[Expr | Function | Constant | float, ...] | None = None,
    # physical 
    epsilon: float | None = None,
    # initial conditions
    c_ics: Function | SpatialPerturbation | None = None,
    theta_ics: Function | SpatialPerturbation | None = None,
    s_ics: Function | SpatialPerturbation | None = None,
    # boundary conditions
    flow_bcs: BoundaryConditions | EllipsisType | None | tuple[BoundaryConditions | None, BoundaryConditions | None] = ...,
    c_bcs: BoundaryConditions | EllipsisType | None = None,
    theta_bcs: BoundaryConditions | EllipsisType | None = None,
    # constitutive relations
    rock_porosity: Callable[[np.ndarray], np.ndarray] | float = 1,
    permeability: Callable[[Phi], Expr | ExprSeries] = lambda phi: phi**2,
    dispersion_solutal: Callable[[Phi, U], Expr | ExprSeries]
    | Callable[[Phi], Expr | ExprSeries] = lambda phi: phi,
    dispersion_thermal: Callable[[Phi, U], Expr | ExprSeries]
    | Callable[[Phi], Expr | ExprSeries] = lambda phi: phi,
    density: Callable[[C | Theta], ExprSeries] 
    | Callable[[C, Theta], ExprSeries] 
    | None = lambda c, theta: c - theta,
    viscosity: Callable[[C | Theta], ExprSeries]
    | Callable[[C, Theta], ExprSeries] 
    | None = None, 
    reaction: Callable[[S], Expr | ExprSeries] 
    | Callable[[S, Theta], ExprSeries]
    | None = None,
    source: Callable[[S], Expr | ExprSeries] 
    | Callable[[S, Theta], ExprSeries] 
    | None = None,
    # time step
    dt_min: float = 0.0,
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float | None = 0.5,
    r_courant: float | None = None,
    # time discretization
    D_adv_solutal: FiniteDifference | FiniteDifferenceArgwise = FE,
    D_diff_solutal: FiniteDifference = FE,
    D_reac_solutal: FiniteDifference 
    | FiniteDifferenceArgwise = FE,
    D_src_solutal: FiniteDifference = FE,
    D_adv_thermal: FiniteDifference | FiniteDifferenceArgwise = FE,
    D_diff_thermal: FiniteDifference = FE,
    D_reac_evol: FiniteDifference 
    | FiniteDifferenceArgwise = FE,
    # stabilization
    c_stabilization: str | tuple[float, float] | None = None,
    c_limits: tuple[float, float] | bool = False,
    theta_stabilization: str | tuple[float, float] = None,
    theta_limits: tuple[float, float] | EllipsisType | None = None,
    s_limits: tuple[float, float] | bool = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    theta_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    s_petsc: OptionsPETSc | None = None,
    # optional postprocessing
    diagnostic: bool | Iterable[Solver] = False,   
    fluxes_solutal: Iterable[tuple[str, float | int, float]] = (), 
    fluxes_thermal: Iterable[tuple[str, float | int, float]] = (), 
    namespace: Iterable[Function | Constant | ExprSeries | tuple[str, Expr]] = (),
) -> Simulation:    
    """
    `𝜑∂s/∂t = -ε(R(s,θ)c + J(s, θ))`
    `ϕ = 𝜑(1 - s)` \\
    `ϕ∂c/∂t + 𝐮·∇c =  ∇·(D(ϕ,𝐮)·∇c) + R(s, θ)c + J(s,θ)` \\
    `ϕ∂θ/∂t + 𝐮·∇θ = ∇·(G(ϕ,𝐮)·∇θ)`\\
    `∇⋅𝐮 = 0` \\
    `𝐮 = -(K(ϕ)/μ(c, θ))·(∇p - ρ(c,θ)e₉)` \\
    """
    if eg is None:
        if Omega.geometry.dim == 2:
            eg = (0, -1)
        if Omega.geometry.dim == 3:
            eg = (0, 0, -1)

    STREAMF = isinstance(flow_petsc, tuple)
    SOLUTAL = c_ics is not None
    THERMAL = theta_ics is not None
    THERMOSOLUTAL = SOLUTAL and THERMAL
    EVOL = epsilon is not None
    SOLUTAL_CNTS = c_stabilization is None or isinstance(c_stabilization, str)
    THERMAL_CNTS = theta_stabilization is None or isinstance(theta_stabilization, str)
    SOLUTAL_MECH = arity(dispersion_solutal) == 2
    THERMAL_MECH = arity(dispersion_thermal) == 2

    solvers = []
    namespace = list(namespace)

    # time
    order = finite_difference_order(
        D_adv_solutal, D_diff_solutal, D_reac_solutal, D_src_solutal,
        D_adv_thermal, D_diff_thermal, D_reac_evol,
    )
    t = ConstantSeries(Omega, "t", order, ics=0.0)  
    dt = ConstantSeries(Omega, 'dt')

    # flow
    if STREAMF: 
        psi_deg = 2
        psi = FunctionSeries((Omega, 'P', psi_deg), 'psi')
        u = FunctionSeries((Omega, "P", psi_deg - 1, 2), "u", order)
    else:
        u_fam = 'BDM' if Omega.topology.cell_name() == CellType.TRIANGLE else 'BDMCF'
        u_deg = 1
        up = FunctionSeries((Omega, [(u_fam, u_deg), ('DP', u_deg - 1)]), ('up', ['u', 'p']), order)
        u = up.sub(0)

    # transport
    if SOLUTAL:
        c = FunctionSeries((Omega, 'P' if SOLUTAL_CNTS else 'DP', 1), 'c', order, ics=c_ics)
    if THERMAL:
        theta = FunctionSeries((Omega, 'P' if THERMAL_CNTS else 'DP', 1), 'theta', order, ics=theta_ics)
    if EVOL:
        s = FunctionSeries((Omega, 'P', 1), 's', order, ics=s_ics)
        epsilon = Constant(Omega, epsilon, 'epsilon')
        namespace.append(epsilon)
    else:
        s = Function((Omega, 'P', 1), s_ics, 's')

    # constitutive
    varphi = Function((Omega, 'P', 1), rock_porosity, 'varphi')
    phi = varphi * (1 - s)
    k = permeability(phi)   
    rho = density(
        *(c, theta) if THERMOSOLUTAL
        else (c, ) if SOLUTAL
        else (theta, )
    ) if density else 1
    mu = viscosity(
        *(c, theta) if THERMOSOLUTAL
        else (c, ) if SOLUTAL
        else (theta, )
    ) if viscosity else 1
    namespace.extend((varphi, ('phi', phi), ('k', k), ('rho', rho), ('mu', mu)))
    if SOLUTAL:
        d = dispersion_solutal(
            *(phi, u) if SOLUTAL_MECH
            else (phi, )
        )
        reaction_eff = lambda c, *args: (
            (reaction(*args) * c if reaction else 0) 
            + (source(*args) if source else 0)
        )
        r = reaction(
            *(s, theta) if THERMOSOLUTAL
            else (s, )
        ) if reaction else 0
        j = source(
            *(s, theta) if THERMOSOLUTAL
            else (s, )
        ) if source else 0
        rEff = reaction_eff(
            *(c, s, theta) if THERMOSOLUTAL
            else (c, s)
        )
        namespace.extend([('d', d), ('r', r), ('j', j), ('rEff', rEff)])
    if THERMAL:
        g = dispersion_thermal(
            *(phi, u) if THERMAL_MECH
            else (phi, )
        )
        namespace.append(('g', g))

    # flow solvers
    if STREAMF:
        psi_bcs = BoundaryConditions(("dirichlet", dOmega.union, 0.0)) if flow_bcs is Ellipsis else flow_bcs
        psi_petsc, u_petsc = flow_petsc
        psi_cache_matrix = not(EVOL or viscosity)
        egx, egy = eg
        psi_solver = bvp(darcy_streamfunction, psi_bcs, petsc=psi_petsc, cache_matrix=psi_cache_matrix)(
            psi, FE(k), FE(mu), egx * FE(rho), egy * FE(rho),
        ) 
        if u_petsc is None:
            u_solver = interpolation(u, streamfunction_velocity)(psi[0])
        else:
            u_solver = projection(u, streamfunction_velocity, petsc=u_petsc)(psi[0])
        
        solvers.extend((psi_solver, u_solver))
    else:
        u_bcs = BoundaryConditions(('essential', dOmega.union, (0.0, 0.0), 0)) if flow_bcs is Ellipsis else flow_bcs[0]
        p_bcs = None if flow_bcs is Ellipsis else flow_bcs[1]
        flow_petsc = flow_petsc.replace(pc_factor_mat_solver_type='mumps')
        u_solver = bvp(darcy_incompressible, u_bcs, petsc=flow_petsc)(
            up, FE(k), FE(mu), FE(rho) * as_vector(eg), p_bcs)
        solvers.append(u_solver)

    # timestep solver
    if SOLUTAL:
        dt_solver = evaluation(dt, cflr_timestep)(
            u[0], FE(rEff), cfl_h, cfl_courant, r_courant, dt_max, dt_min,
        ) 
    else:
        dt_solver = evaluation(dt, cfl_timestep)(
            u[0], cfl_h, cfl_courant, dt_max, dt_min, 
        )
    solvers.append(dt_solver)

    # thermal solver
    if THERMAL:
        theta_bcs = BoundaryConditions(("neumann", dOmega.union, 0.0)) if theta_bcs is Ellipsis else theta_bcs
        theta_limits = (0, 1) if c_limits is True else c_limits
        theta_corrector = ('thetaCorr', limits_corrector(*theta_limits)) if theta_limits else None
        if THERMAL_CNTS:
            theta_solver = ibvp(advection_diffusion, bcs=theta_bcs, petsc=theta_petsc, corrector=theta_corrector)(
                theta, dt, u, g, D_adv_thermal, D_diff_thermal, phi=phi, supg=theta_stabilization,
            )
        else:
            theta_alpha, theta_gamma = theta_stabilization
            theta_solver = ibvp(advection_diffusion_dg, petsc=theta_petsc, corrector=theta_corrector)(
                theta, dt, u, g, theta_alpha, theta_gamma, D_adv_thermal, D_diff_thermal, phi=phi, bcs=theta_bcs,
            )
        solvers.append(theta_solver)

    # solutal solver
    if SOLUTAL:
        c_bcs = BoundaryConditions(("neumann", dOmega.union, 0.0)) if c_bcs is Ellipsis else c_bcs
        c_limits = (0, 1) if c_limits is True else c_limits
        c_corrector = ('cCorr', limits_corrector(*c_limits)) if c_limits else None
        if SOLUTAL_CNTS:
            c_solver = ibvp(advection_diffusion_reaction, bcs=c_bcs, petsc=c_petsc, corrector=c_corrector)(
                c, dt, u, d, r, j, D_adv_solutal, D_diff_solutal, D_reac_solutal, D_src_solutal, phi=phi, supg=c_stabilization,
            )
        else:
            c_alpha, c_gamma = c_stabilization
            c_solver = ibvp(advection_diffusion_reaction_dg, petsc=c_petsc, corrector=c_corrector)(
                c, dt, u, d, r, j, c_alpha, c_gamma, D_adv_solutal, D_diff_solutal, D_reac_solutal, phi=phi, bcs=c_bcs,
            )
        solvers.append(c_solver)

    # evolution solver
    if EVOL:
        reaction_evol = lambda c, *args: -epsilon * reaction_eff(c, *args)
        rEvol = ExprSeries.from_expr_func(reaction_evol, name='rEvol')(
            *(c, s, theta) if THERMOSOLUTAL
            else (c, s)
        )
        namespace.append(rEvol)
        s_limits = (np.min(s.ics.x.array), np.max(s.ics.x.array)) if s_limits is True else s_limits
        s_corrector = ('sCorr', limits_corrector(*s_limits)) if s_limits else None
        if s_petsc is None:
            s_solver = interpolation(s, evolution_update, corrector=s_corrector, future=True)(
                s, dt, rEvol, D_reac_evol, phi=varphi, explicit=1,
            )
        else:
            s_solver = ivp(evolution, petsc=s_petsc, corrector=s_corrector)(
                s, dt, rEvol, D_reac_evol, phi=varphi
            )
        solvers.append(s_solver)

    # optional post-processing
    if diagnostic:
        uMinMax = ConstantSeries(Omega, "uMinMax", shape=(2,))
        solvers.append(evaluation(uMinMax, extrema)(u[0]))
        norm = 2
        uRMS = ConstantSeries(Omega, 'uRMS')
        solvers.append(integration(uRMS, L_norm, 'dx', norm=norm)(sqrt(inner(u[0], u[0])), norm))
        uDiv = ConstantSeries(Omega, 'uDiv')
        solvers.append(integration(uDiv, div_norm, 'dx', norm=norm)(u[0], norm))
        dtCFL = ConstantSeries(Omega, "dtCFL")
        solvers.append(evaluation(dtCFL, cfl_timestep)(u[0], cfl_h))
        if THERMAL:
            thetaMinMax = ConstantSeries(Omega, "thetaMinMax", shape=(2,))
            solvers.append(evaluation(thetaMinMax, extrema)(theta[0]))
            qBoundary = ConstantSeries(Omega, "qBoundary", shape=(len(dOmega.markers), 2))
            solvers.append(integration(qBoundary, flux, 'ds', *dOmega.markers)(theta[0], u[0], FE(g)))
        if SOLUTAL:
            cMinMax = ConstantSeries(Omega, "cMinMax", shape=(2,))
            solvers.append(evaluation(cMinMax, extrema)(c[0]))
            mD = ConstantSeries(Omega, "mD")
            solvers.append(integration(mD, mass_dissolved, 'dx')(c[0], FE(phi)))
            fBoundary = ConstantSeries(Omega, "fBoundary", shape=(len(dOmega.markers), 2))
            solvers.append(integration(fBoundary, flux, 'ds', *dOmega.markers)(c[0], u[0], FE(d)))
        if EVOL:
            sMinMax = ConstantSeries(Omega, "sMinMax", shape=(2,))
            solvers.append(evaluation(sMinMax, extrema)(s[0]))
            mC = ConstantSeries(Omega, "mC")
            solvers.append(integration(mC, mass_capillary_trapped, 'dx')(s[0], epsilon))
            dtK = ConstantSeries(Omega, "dtK")
            solvers.append(evaluation(dtK, reactive_timestep)(rEff[0]))
        if isinstance(diagnostic, Iterable):
            solvers.extend(diagnostic)

    for i in (*fluxes_solutal, *fluxes_thermal):
        name, y_target, l = i
        f = ConstantSeries(
            Omega, 
            (name, (name, f'{name}Plus', f'{name}Minus')), 
            shape=(3, 2),
        )
        if i in fluxes_solutal:
           _c = c
           _d = d
        else:
            _c = theta
            _d = g 
        solvers.append(integration(f, vertical_flux)(_c[0], u[0], FE(_d), y_target, l))
    
    return Simulation(solvers, t, dt, namespace)