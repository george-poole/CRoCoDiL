from collections.abc import Iterable
from typing import Callable, TypeAlias
from types import EllipsisType

import numpy as np
from ufl import as_vector, inner, sqrt, min_value
from ufl.core.expr import Expr
from dolfinx.mesh import Mesh

from lucifex.fem import Function, Constant, SpatialPerturbation
from lucifex.mesh import MeshBoundary
from lucifex.fdm import (
    FunctionSeries, ConstantSeries, FiniteDifference, FE, 
    ExprSeries, FiniteDifferenceArgwise, finite_difference_order,
    adr_timestep, advective_timestep, advective_diffusive_timestep, 
    diffusive_timestep, reactive_timestep,
)
from lucifex.fdm.ufl_operators import max_value
from lucifex.solver import (
    BoundaryConditions, OptionsPETSc, Solver,
    bvp, ibvp, ivp, 
    interpolation, projection, evaluation, integration,
    extrema, div_norm, L_norm,
)
from lucifex.sim import Simulation
from lucifex.pde.streamfunction_vorticity import velocity_from_streamfunction
from lucifex.pde.darcy import darcy_streamfunction, darcy_incompressible
from lucifex.pde.advection_diffusion import (
    dg_advection_diffusion_reaction, advection_diffusion, 
    dg_advection_diffusion, advection_diffusion_reaction, flux,
)
from lucifex.utils.fenicsx_utils import is_simplicial, limits_corrector
from lucifex.utils.py_utils import arity
from lucifex.pde.evolution import evolution, evolution_rhs

from .diagnostic import mass_dissolved, mass_capillary, vertical_flux


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
    # timestep
    dt_min: float = 0.0,
    dt_max: float = 0.5,
    dt_h: str | float = "hmin",
    courant_adv: float | None = 1.0,
    courant_diff: float | None = 1.0,
    courant_reac: float | None = 1.0,
    # time discretization
    D_adv_solutal: FiniteDifference | FiniteDifferenceArgwise = FE,
    D_diff_solutal: FiniteDifference = FE,
    D_reac_solutal: FiniteDifference 
    | FiniteDifferenceArgwise = FE,
    D_src_solutal: FiniteDifference = FE,
    D_adv_thermal: FiniteDifference | FiniteDifferenceArgwise = FE,
    D_diff_thermal: FiniteDifference = FE,
    D_evol: FiniteDifference 
    | FiniteDifferenceArgwise = FE,
    # stabilization
    c_stabilization: str | float | tuple[float, float] | None = None,
    c_limits: tuple[float, float] | bool = False,
    theta_stabilization: str | float | tuple[float, float] = None,
    theta_limits: tuple[float, float] | EllipsisType | None = None,
    s_elem: tuple[str, int] = ('P', 1),
    s_limits: tuple[float, float] | bool = False,
    phi_elem: tuple[str, int] = ('P', 1),
    use_max_value: bool = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'hypre'), None),
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
    `𝜑∂s/∂t = -ε(R(θ,s)c + J(θ,s))`
    `ϕ = 𝜑(1 - s)` \\
    `ϕ∂c/∂t + 𝐮·∇c =  ∇·(D(ϕ,𝐮)·∇c) + R(θ,s)c + J(θ,s)` \\
    `ϕ∂θ/∂t + 𝐮·∇θ = ∇·(G(ϕ,𝐮)·∇θ)`\\
    `∇⋅𝐮 = 0` \\
    `𝐮 = -(K(ϕ)/μ(c,θ))·(∇p - ρ(c,θ)e₉)` \\
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
    is_cg = lambda arg: isinstance(arg, (str, type(None)))
    SOLUTAL_CG = is_cg(c_stabilization)
    THERMAL_CG = is_cg(theta_stabilization)
    SOLUTAL_DISP = arity(dispersion_solutal) == 2
    THERMAL_DISP = arity(dispersion_thermal) == 2
    HETEROGENEOUS = callable(rock_porosity)
    if THERMOSOLUTAL:
        assert arity(density) == 2  if density else True, "Density should be `ρ(c,θ)`"
        assert arity(viscosity) == 2  if viscosity else True, "Viscosity should be `μ(c,θ)`"
        assert arity(reaction) == 2  if reaction else True, "Reaction should be `R(θ,s)`"
        assert arity(source) == 2  if source else True, "Source should be `J(θ,s)`"

    solvers = []
    namespace = list(namespace)

    # time
    order = finite_difference_order(
        D_adv_solutal, D_diff_solutal, D_reac_solutal, D_src_solutal,
        D_adv_thermal, D_diff_thermal, D_evol,
    )
    t = ConstantSeries(Omega, "t", order, ics=0.0)  
    dt = ConstantSeries(Omega, 'dt')

    # flow
    if STREAMF: 
        psi_deg = 2
        psi = FunctionSeries((Omega, 'P', psi_deg), 'psi')
        u = FunctionSeries((Omega, "P", psi_deg - 1, 2), "u", order)
    else:
        u_fam = 'BDM' if is_simplicial(Omega) else 'BDMCF'
        u_deg = 1
        up = FunctionSeries((Omega, [(u_fam, u_deg), ('DP', u_deg - 1)]), ('up', ['u', 'p']), order)
        u = up.sub(0)

    # transport
    if SOLUTAL:
        c = FunctionSeries((Omega, 'P' if SOLUTAL_CG else 'DP', 1), 'c', order, ics=c_ics)
    if THERMAL:
        theta = FunctionSeries((Omega, 'P' if THERMAL_CG else 'DP', 1), 'theta', order, ics=theta_ics)
    if EVOL:
        s = FunctionSeries((Omega, *s_elem), 's', order, ics=s_ics)
        epsilon = Constant(Omega, epsilon, 'epsilon')
        namespace.append(epsilon)
    else:
        s = Function((Omega, *s_elem), s_ics, 's')

    # constitutive
    if HETEROGENEOUS:
        varphi = Function((Omega, *phi_elem), rock_porosity, 'varphi')
    else:
        varphi = Constant(Omega, rock_porosity, 'varphi')
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
            *(phi, u) if SOLUTAL_DISP
            else (phi, )
        )
        r = reaction(
            *(theta, s) if THERMOSOLUTAL and arity(reaction) == 2
            else (s, )
        ) if reaction else 0
        j = source(
            *(theta, s) if THERMOSOLUTAL and arity(source) == 2
            else (s, )
        ) if source else 0
        if THERMOSOLUTAL:
            reaction_plus_source = lambda c, theta, s: (
                (reaction(theta, s) * c if reaction else 0) 
                + (source(theta, s) if source else 0)
            )
        else:
            reaction_plus_source = lambda c, s: (
                (reaction(s) * c if reaction else 0) 
                + (source(s) if source else 0)
            )
        if use_max_value:
            _reaction_plus_source = lambda *args: (
                max_value(reaction_plus_source(*args), 0)
            )
            reaction_plus_source = _reaction_plus_source
        Sigma = reaction_plus_source(
            *(c, theta, s) if THERMOSOLUTAL
            else (c, s)
        )
        namespace.extend([('d', d), ('r', r), ('j', j), ('Sigma', Sigma)])
    if THERMAL:
        g = dispersion_thermal(
            *(phi, u) if THERMAL_DISP
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
            u_solver = interpolation(u, velocity_from_streamfunction)(psi[0])
        else:
            u_solver = projection(u, velocity_from_streamfunction, petsc=u_petsc)(psi[0])
        
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
        if THERMAL:
            dt_disp = min_value(FE(d), FE(g))
        else:
            dt_disp = FE(d)
        dt_solver = evaluation(dt, adr_timestep)(
            u[0], dt_disp, FE(Sigma), dt_h, courant_adv, courant_diff, courant_reac, dt_max, dt_min,
        ) 
    else:
        dt_solver = evaluation(dt, advective_diffusive_timestep)(
            u[0], FE(g), dt_h, courant_adv, courant_diff, dt_max, dt_min, 
        )
    solvers.append(dt_solver)

    # thermal solver
    if THERMAL:
        theta_bcs = BoundaryConditions(("neumann", dOmega.union, 0.0)) if theta_bcs is Ellipsis else theta_bcs
        theta_limits = (0, 1) if theta_limits is True else theta_limits
        theta_corrector = ('thetaCorr', limits_corrector(*theta_limits)) if theta_limits else None
        if THERMAL_CG:
            theta_solver = ibvp(advection_diffusion, bcs=theta_bcs, petsc=theta_petsc, corrector=theta_corrector)(
                theta, dt, u, g, D_adv_thermal, D_diff_thermal, phi=phi, supg=theta_stabilization,
            )
        else:
            theta_solver = ibvp(dg_advection_diffusion, petsc=theta_petsc, corrector=theta_corrector)(
                theta, dt, theta_stabilization, u, g, D_adv_thermal, D_diff_thermal, phi=phi, bcs=theta_bcs,
            )
        solvers.append(theta_solver)

    # solutal solver
    if SOLUTAL:
        c_bcs = BoundaryConditions(("neumann", dOmega.union, 0.0)) if c_bcs is Ellipsis else c_bcs
        c_limits = (0, 1) if c_limits is True else c_limits
        c_corrector = ('cCorr', limits_corrector(*c_limits)) if c_limits else None
        if SOLUTAL_CG:
            c_solver = ibvp(advection_diffusion_reaction, bcs=c_bcs, petsc=c_petsc, corrector=c_corrector)(
                c, dt, u, d, r, j, D_adv_solutal, D_diff_solutal, D_reac_solutal, D_src_solutal, phi=phi, supg=c_stabilization,
            )
        else:
            c_solver = ibvp(dg_advection_diffusion_reaction, petsc=c_petsc, corrector=c_corrector)(
                c, dt, c_stabilization, u, d, r, j, D_adv_solutal, D_diff_solutal, D_reac_solutal, phi=phi, bcs=c_bcs,
            )
        solvers.append(c_solver)

    # evolution solver
    if EVOL:
        if THERMOSOLUTAL:
            evol = lambda c, theta, s: -epsilon * reaction_plus_source(c, theta, s)
        else:
            evol = lambda c, s: -epsilon * reaction_plus_source(c, s)
        SigmaEvol = ExprSeries.from_expr_factory(evol, name='SigmaEvol')(
            *(c, theta, s) if THERMOSOLUTAL
            else (c, s)
        )
        namespace.append(SigmaEvol)
        s_limits = (np.min(s.ics.x.array), np.max(s.ics.x.array)) if s_limits is True else s_limits
        s_corrector = ('sCorr', limits_corrector(*s_limits)) if s_limits else None
        if s_petsc is None:
            s_solver = interpolation(s, evolution_rhs, corrector=s_corrector, future=True)(
                s, dt, SigmaEvol, D_evol, phi=varphi, argwise_index=-1,
            )
        else:
            s_solver = ivp(evolution, petsc=s_petsc, corrector=s_corrector)(
                s, dt, SigmaEvol, D_evol, phi=varphi
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
        dtU = ConstantSeries(Omega, "dtU")
        solvers.append(evaluation(dtU, advective_timestep)(u[0], dt_h))
        if THERMAL:
            thetaMinMax = ConstantSeries(Omega, "thetaMinMax", shape=(2,))
            solvers.append(evaluation(thetaMinMax, extrema)(theta[0]))
            qBoundary = ConstantSeries(Omega, "qBoundary", shape=(len(dOmega.markers), 2))
            solvers.append(integration(qBoundary, flux, 'ds', *dOmega.markers)(theta[0], u[0], FE(g)))
            dtG = ConstantSeries(Omega, "dtG")
            solvers.append(evaluation(dtG, diffusive_timestep)(FE(g), dt_h))
        if SOLUTAL:
            mD = ConstantSeries(Omega, "mD")
            solvers.append(integration(mD, mass_dissolved, 'dx')(c[0], FE(phi)))
            cMinMax = ConstantSeries(Omega, "cMinMax", shape=(2,))
            solvers.append(evaluation(cMinMax, extrema)(c[0]))
            fBoundary = ConstantSeries(Omega, "fBoundary", shape=(len(dOmega.markers), 2))
            solvers.append(integration(fBoundary, flux, 'ds', *dOmega.markers)(c[0], u[0], FE(d)))
            dtD = ConstantSeries(Omega, "dtD")
            solvers.append(evaluation(dtD, diffusive_timestep)(FE(d), dt_h))
        if EVOL:
            mC = ConstantSeries(Omega, "mC")
            solvers.append(integration(mC, mass_capillary, 'dx')(s[0], varphi, epsilon))
            sMinMax = ConstantSeries(Omega, "sMinMax", shape=(2,))
            solvers.append(evaluation(sMinMax, extrema)(s[0]))
            dtSigma = ConstantSeries(Omega, "dtSigma")
            solvers.append(evaluation(dtSigma, reactive_timestep)(Sigma[0]))
        if isinstance(diagnostic, Iterable):
            solvers.extend(diagnostic)

    for i in (*fluxes_solutal, *fluxes_thermal):
        name, y_target, lx = i
        f = ConstantSeries(
            Omega, 
            (name, (name, f'{name}Plus', f'{name}Minus')), 
            shape=(3, 2),
        )
        vertical_flux_solver = lambda c, d: (
            integration(f, vertical_flux)(c[0], u[0], FE(d), y_target, lx)
        )
        if i in fluxes_solutal:
           solvers.append(vertical_flux_solver(c, d))
        else:
            solvers.append(vertical_flux_solver(theta, g))
    
    return Simulation(solvers, t, dt, namespace)