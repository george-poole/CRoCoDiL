from typing import Callable, Iterable, TypeAlias
from types import EllipsisType

import numpy as np
from ufl.core.expr import Expr
from dolfinx.mesh import Mesh

from lucifex.fem import LUCiFExFunction as Function, LUCiFExConstant as Constant
from lucifex.mesh import MeshBoundary
from lucifex.fdm import (
    FunctionSeries, ConstantSeries, FiniteDifference, AB1, Series, 
    ExprSeries, finite_difference_order,
    cflr_timestep, cfl_timestep, reactive_timestep,
)
from lucifex.solver import (
    BoundaryConditions, OptionsPETSc, eval_solver, 
    dx_solver, ds_solver,
)
from lucifex.utils import CellType, extremum
from lucifex.sim import Simulation

from ..math.secondary import flux, mass_dissolved, mass_capillary_trapped
from ..math.solvers import (
    streamfunction_solvers, advection_diffusion_solver,
    advection_diffusion_reaction_solver, evolution_solver,
    darcy_solver, streamfunction_method, continuous)


C: TypeAlias = FunctionSeries
Theta: TypeAlias = FunctionSeries
S: TypeAlias = FunctionSeries
Phi: TypeAlias = FunctionSeries
U: TypeAlias = FunctionSeries
def abstract_simulation(
    # domain
    Omega: Mesh,
    dOmega: MeshBoundary,
    egx: Expr | Function | Constant | float | None = None,
    egy: Expr | Function | Constant | float | None = None,
    egz: Expr | Function | Constant | float | None = None,
    # physical 
    Ra: float | None = None,
    Rb: float | None = None,
    Da: float | None = None,
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
    varphi: Callable[[np.ndarray], np.ndarray] | float = 1,
    permeability: Callable[[Phi], Series] = lambda phi: phi**2,
    dispersion_solutal: Callable[[Phi, U], Series] = lambda phi, _: phi,
    dispersion_thermal: Callable[[Phi, U], Series] = lambda phi, _: phi,
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
    D_adv_solutal: FiniteDifference | tuple[FiniteDifference, FiniteDifference] = AB1,
    D_diff_solutal: FiniteDifference = AB1,
    D_reac_solutal: FiniteDifference 
    | tuple[FiniteDifference, FiniteDifference]
    | tuple[FiniteDifference, FiniteDifference, FiniteDifference] = AB1,
    D_adv_thermal: FiniteDifference | tuple[FiniteDifference, FiniteDifference] = AB1,
    D_diff_thermal: FiniteDifference = AB1,
    D_reac_evol: FiniteDifference 
    | tuple[FiniteDifference, FiniteDifference]
    | tuple[FiniteDifference, FiniteDifference, FiniteDifference] = AB1,
    # stabilization
    c_stabilization: str | tuple[float, float] = None,
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
    Default constitutive relations are uniform rock porosity `𝜑 = 1`, 
    isotropic quadratic permeability `K(ϕ)=ϕ²`, isotropic linear solutal
    dispersion `D(ϕ)=ϕ`, isotropic linear thermal dispersion `G(ϕ)=ϕ`, 
    uniform viscosity `μ = 1`.

    Default gravity unit vector is `e₉ = -eʸ` in 2D or `e₉ = -eᶻ` in 3D.
    
    Default boundary conditions are no flux of fluid, solute and heat everywhere on `∂Ω`. 

    In general: rock porosity is a spatial function `𝜑(𝐱)`; permeability `K(ϕ)` is a 
    function of porosity; solutal and thermal dispersions are functions `D(ϕ, 𝐮)`, `G(ϕ, 𝐮)` of porosity
    and velocity; density and viscosity are functions `μ(c, θ)`, `ρ(c, θ)` of concentration
    and temperature; reaction rate is a function `R(s, c, θ)` of saturation, concentration
    and temperature.

    `ϕ = 𝜑(1 - s)` is the effective porosity.
    """
    if all(i is None for i in (egx, egy, egz)):
        if Omega.geometry.dim == 2:
            egy = -1
        if Omega.geometry.dim == 3:
            egx = 0
            egy = 0
            egz = -1

    RA, DA, EPS, RB = (bool(_) for _ in (Ra, Da, epsilon, Rb))
    PSI = streamfunction_method(flow_petsc)
    order = finite_difference_order(D_adv_solutal, D_diff_solutal, D_reac_solutal, D_adv_thermal, D_diff_thermal, D_reac_evol)
    solvers = []
    namespace = []

    # flow fields
    if PSI: 
        psi_deg = 2
        psi = FunctionSeries((Omega, 'P', psi_deg), 'psi')
        u = FunctionSeries((Omega, "P", psi_deg - 1, 2), "u", order)
    else:
        u_fam = 'BDMCF' if Omega.topology.cell_name() == CellType.QUADRILATERAL else 'BDM'
        u_deg = 1
        up = FunctionSeries((Omega, [(u_fam, u_deg), ('DP', u_deg - 1)]), "up", order)
        u = up.sub(0, "u")

    # transport fields
    t = ConstantSeries(Omega, "t", order, ics=0.0)  
    dt = ConstantSeries(Omega, 'dt')
    if RA:
        Ra = Constant(Omega, Ra, 'Ra')
        c = FunctionSeries((Omega, 'P' if continuous(c_stabilization) else 'DP', 1), 'c', order, ics=c_ics)
    if RB:
        Rb = Constant(Omega, Rb, 'Rb')
        theta = FunctionSeries((Omega, 'P' if continuous(theta_stabilization) else 'DP', 1), 'theta', order, ics=theta_ics)
    if DA:
        Da = Constant(Omega, Da, 'Da')
    if EPS:
        epsilon = Constant(Omega, epsilon, 'epsilon')
        namespace.append(epsilon)
    if DA and EPS:
        s = FunctionSeries((Omega, 'P', 1), 's', order, ics=s_ics)
    else:
        s = Function((Omega, 'P', 1), 's', s_ics)

    # constitutive relations
    varphi = Function((Omega, 'P', 1), varphi, 'varphi')
    phi = ExprSeries(varphi * (1 - s), 'phi')
    k = ExprSeries(permeability(phi), 'k')
    if RA:
        d = ExprSeries(dispersion_solutal(phi, u), 'd')
        namespace.extend((Ra, d))
    if RB:
        g = ExprSeries(dispersion_thermal(phi, u), 'g')
        namespace.extend((Rb, g))
    if DA:
        r = ExprSeries.from_relation(reaction, 'r')(
            *(s, c, theta) if RB
            else (s, c)
        )
        namespace.extend((Da, r))
    rho = ExprSeries(
        density(c, theta) if RA and RB
        else density(c) if RA
        else density(theta),
        'rho',
    )
    mu = ExprSeries(
        viscosity(c, theta) if RA and RB
        else viscosity(c) if RA
        else viscosity(theta),
        'mu',
    )
    namespace.extend((varphi, phi, k, rho, mu))

    # flow solvers
    if PSI:
        psi_bcs = BoundaryConditions(("dirichlet", dOmega.union, 0.0)) if flow_bcs is Ellipsis else flow_bcs
        solvers.extend(
            streamfunction_solvers(psi, u, rho[0], k[0], mu[0], egx, egy, psi_bcs, flow_petsc)
        )
    else:
        u_bcs = BoundaryConditions(('essential', dOmega.union, (0.0, 0.0), 0)) if flow_bcs is Ellipsis else flow_bcs[0]
        p_bcs = None if flow_bcs is Ellipsis else flow_bcs[1]
        solvers.append(
            darcy_solver(up, rho[0], k[0], mu[0], egx, egy, egz, u_bcs, p_bcs)
        )

    # timestep solver
    if DA:
        cflr_r = Da * r[0]
    else:
        cflr_r = 0.0
    solvers.append(
        eval_solver(dt, cflr_timestep)(
            u[0], cflr_r, cfl_h, cfl_courant, k_courant, dt_max, dt_min,
        ) 
    )

    # transport solvers
    if RB:
        theta_bcs = BoundaryConditions(("neumann", dOmega.union, 0.0)) if theta_bcs is Ellipsis else theta_bcs
        theta_limits = (0, 1) if theta_limits is Ellipsis else theta_limits
        solvers.append(
            advection_diffusion_solver(
                theta, dt, phi, u, Ra, g, D_adv_thermal, D_diff_thermal, theta_bcs, theta_petsc, theta_limits, theta_stabilization,
            )
        )
    if RA:
        c_bcs = BoundaryConditions(("neumann", dOmega.union, 0.0)) if c_bcs is Ellipsis else c_bcs
        c_limits = (0, 1) if c_limits is Ellipsis else c_limits
        solvers.append(
            advection_diffusion_reaction_solver(
                c, dt, phi, u, Ra, d, Da, r, D_adv_solutal, D_diff_solutal, D_reac_solutal, c_bcs, c_petsc, c_limits, c_stabilization,
            )
        )
    # evolution solver
    if EPS and DA:
        s_limits = (0, max(s.ics.x)) if s_limits is Ellipsis else s_limits
        solvers.append(
            evolution_solver(s, dt, varphi, epsilon, Da, r, D_reac_evol, s_limits, s_petsc)
        )

    # optional solvers
    if secondary:
        solvers.append(
            eval_solver(ConstantSeries(Omega, "uMinMax", shape=(2,)), extremum)(u[0])
        )
        solvers.append(
            eval_solver(ConstantSeries(Omega, "dtCFL"), cfl_timestep)(u[0], cfl_h)
        )
        if DA:
            solvers.append(
                eval_solver(ConstantSeries(Omega, "dtK"), reactive_timestep)(cflr_r)
            )
        if RA:
            solvers.append(
                eval_solver(ConstantSeries(Omega, "cMinMax", shape=(2,)), extremum)(c[0])
            )
            solvers.append(
                dx_solver(ConstantSeries(Omega, "mD"), mass_dissolved)(c[0], s[0])
            )
            solvers.append(
                ds_solver(
                    ConstantSeries(Omega, "fOmega", shape=(len(dOmega.union), 2)), 
                    flux, 
                    dOmega.union,
                )(c[0], u[0], d[0], Ra)
            )
            if EPS:
                solvers.append(
                    dx_solver(ConstantSeries(Omega, "mC"), mass_capillary_trapped)(s[0], epsilon)
                )
        if RB:
            solvers.append(
                eval_solver(ConstantSeries(Omega, "thetaMinMax", shape=(2,)), extremum)(c[0])
            )
            solvers.append(
                ds_solver(
                    ConstantSeries(Omega, "jOmega", shape=(len(dOmega.union), 2)), 
                    flux,
                    dOmega.union,
                )(c[0], u[0], g[0], Rb)
            )

    solvers.extend(secondary_extras)
    namespace.extend(namespace_extras)
    
    return Simulation(solvers, t, dt, namespace)

