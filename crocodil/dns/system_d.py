from typing import Callable, Iterable, Literal

import numpy as np
import scipy.special as sp
from scipy.integrate import quad
from mpi4py import MPI

from lucifex.mesh import rectangle_mesh, mesh_boundary
from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.sim import configure_simulation
from lucifex.solver import OptionsPETSc, OptionsJIT, BoundaryConditions
from lucifex.utils.fenicsx_utils import CellType, limits_corrector
from lucifex.utils.py_utils import FrozenDict

from .generic import dns_generic
from .utils import CROCODIL_JIT_DIR, SCALINGS


SYSTEM_D_REFERENCE = FrozenDict(
    Ra=1000.0,
    Da=100.0,
    epsilon=1e-2,
    Le=1.0,
    Pe=1000.0,
    sr=0.2,
    cr=0.0,
    gamma=1.0,
    delta=0,
    l=0.1,
    inflow=None,
    p_dOmega=None,
    aspect=2.0,
    mirror=0,
)


@configure_simulation(
    jit=OptionsJIT(CROCODIL_JIT_DIR),
    dir_root="./data",
)
def dns_system_d(
    # mesh
    comm: MPI.Comm | str = 'COMM_WORLD',
    aspect: float = 2.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    mirror: Literal[0, 1, -1] = 0,
    # physical
    scaling: str = 'advective',
    Ra: float = 1e3,
    Da: float = 1e2,
    epsilon: float = 1e-2,
    Le: float = 1.0,
    Pe: float = 1e3,
    # initial saturation
    sr: float = 0.2,
    s_ampl: float = 0,
    s_freq: tuple[int, int] = (16, 16),
    s_seed: tuple[int, int] = (1234, 5678),
    # initial concentration
    cr: float = 1.0,
    c_ampl: float = 0,
    c_freq: tuple[int, int] = (16, 16),
    c_seed: tuple[int, int] = (1234, 5678),
    # initial temperature
    theta_buoy: bool = True,
    theta_ampl: float = 0,
    theta_freq: tuple[int, int] = (16, 16),
    theta_seed: tuple[int, int] = (1234, 5678),
    # constitutive relations
    delta: float = 0,
    gamma: float = 1.0,
    permeability: Callable | None = lambda phi: phi**2,
    # inflow profile
    l: float = 0.1,
    inflow: tuple[str, float]  | tuple[str, float, float] 
    | tuple[Callable, float] | None = None,
    p_dOmega: str | tuple[str, ...] | None = None,
    # timestep
    dt_min: float = 0.0,
    dt_max: float = np.inf,
    dt_h: str | float = "hmin",
    dt_Cu: float | None = 0.5,
    dt_Cd: float | None = 0.5,
    dt_Cr: float | None = 0.1,
    # time discretization
    D_adv_solutal: FiniteDifference 
    | FiniteDifferenceArgwise = (AB(1) @ CN),
    D_diff_solutal: FiniteDifference
    | FiniteDifferenceArgwise = (AB(1) @ CN),
    D_adv_thermal: FiniteDifference 
    | FiniteDifferenceArgwise = (AB(1) @ CN),
    D_diff_thermal: FiniteDifference
    | FiniteDifferenceArgwise = (AB(1) @ CN),
    D_reac: FiniteDifference 
    | FiniteDifferenceArgwise = (AB(1) @ CN),
    D_src: FiniteDifference = AB(1),
    D_evol: FiniteDifference 
    | FiniteDifferenceArgwise = (AM(1) @ AM(1) @ AB(1)),
    # stabilization
    c_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    c_limits: bool = False,
    theta_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    theta_limits: bool = False,
    s_limits: bool = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('gmres', 'ilu'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    theta_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    s_petsc: OptionsPETSc | None = None,
    # optional post-processing
    diagnostic: bool = False, 
    fluxes_solutal: Iterable[tuple[str, float | int, float]] = (), 
    fluxes_thermal: Iterable[tuple[str, float | int, float]] = (), 
):
    # scaling
    scaling_map = SCALINGS[scaling](Ra, Da)
    X = scaling_map['X']
    Lx = 0.5 * aspect * X
    if mirror:
        Lx = (0, Lx) if mirror == 1 else (-Lx, 0)
    else:
        Lx = (-Lx, Lx)
    Ly = 1.0 * X
    Ll = l * X
    # inflow conditions
    In = (Pe / Ra) * scaling_map['Bu']
    match inflow:
        case ('gaussian', n):
            xInflow = 0.5 * n * Ll
            uInflow = lambda x: In * np.exp(-(x[0] / (0.5 * Ll))**2)
            uOutflow = 0.25 * np.sqrt(np.pi) * In * Ll * sp.erf(Lx[1] / Ll) # TODO check
        case ('sigmoid', eps, n):
            xInflow = 0.5 * Ll + n / eps
            uInflow = lambda x: In * (
                (1 + np.exp(-(0.5 * Ll + x[0]) / eps)) 
                * (1 + np.exp(-(0.5 * Ll - x[0]) / eps))
            )
            uOutflow = quad(uInflow, Lx[0], Lx[1])[0]
        case None:
            xInflow = 0.5 * Ll
            uInflow = lambda x: np.full_like(x[0], In)
            uOutflow = 0.5 * In * l
        case other:
            uInflow, xInflow = other
            uOutflow = quad(uInflow, Lx[0], Lx[1])[0]
    if mirror:
        xInflow = (0, xInflow) if mirror == 1 else (-xInflow, 0)
    else:
        xInflow = (-xInflow, xInflow)
    is_inflow = lambda x: ((x[0] >= xInflow[0]) & (x[0] <= xInflow[1]))
    uLower = lambda x: 0.0 + (uInflow(x)) * is_inflow(x)    
    # space
    Omega = rectangle_mesh(Lx, Ly, Nx, Ny, cell, 'Omega', comm)
    dOmega = mesh_boundary(
        Omega,
        {
            'left': lambda x: x[0] - Lx[0],
            'right': lambda x: x[0] - Lx[-1],
            'upper': lambda x: x[1] - Ly,
            'lower': lambda x: x[1],
            'lower_imperm': lambda x: (
                np.logical_not(is_inflow(x)) & np.isclose(x[1], 0)
            ),
            'lower_inflow': lambda x: (
                is_inflow(x) & np.isclose(x[1], 0)
            ),
        },
    )
    # constants
    Di, Bu, Ki = scaling_map[Omega, 'Di', 'Bu', 'Ki']
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')
    Le = Constant(Omega, Le, 'Le')
    # initial conditions
    s_ics = sr
    if s_ampl:
        s_ics = SpatialPerturbation(
            s_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], s_freq, s_seed),
            [Lx, Ly],
            s_ampl,
            limits_corrector(0, sr),
        )
    c_ics = cr
    if c_ampl:
        c_ics = SpatialPerturbation(
            c_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed),
            [Lx, Ly],
            c_ampl,
            limits_corrector(0, 1),
        )  
    if theta_buoy:
        theta_ics = 0.0
        theta_dbc = 1.0
    else:
        theta_ics = 1.0 - theta_ampl
        theta_dbc = 0.0
    if theta_ampl:
        theta_ics = SpatialPerturbation(
            theta_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], theta_freq, theta_seed),
            [Lx, Ly],
            theta_ampl,
            limits_corrector(0, 1),
        ) 
    # boundary conditions
    c_bcs = BoundaryConditions(
        # ('neumann', dOmega['left', 'right', 'upper', 'lower'], 0.0),
        ('dirichlet', dOmega['lower_inflow'], 0.0),
    )
    theta_bcs = BoundaryConditions(
        # ('neumann', dOmega['left', 'right', 'upper', 'lower'], 0.0),
        ('dirichlet', dOmega['lower_inflow'], theta_dbc),

    )
    if isinstance(flow_petsc, tuple):
        if p_dOmega:
            raise RuntimeError(
                'Cannot apply pressure boundary condition in streamfunction formulation.'
            )
        psi_bcs = BoundaryConditions(
            ('dirichlet', dOmega['upper', 'lower'], ...), #TODO
            ('dirichlet', dOmega['left', 'right'], ...),
            ('dirichlet', dOmega['inflow'], ...),
        )
        flow_bcs = psi_bcs
    else:
        if not p_dOmega:
            u_bcs = BoundaryConditions(
                ('essential', dOmega['upper'], (0.0, 0.0), 0),
                ('essential', dOmega['left'], (0.0 if mirror == +1 else -uOutflow, 0.0), 0),
                ('essential', dOmega['right'], (0.0 if mirror == -1 else uOutflow, 0.0), 0),
                ('essential', dOmega['lower'], (0.0, uLower), 0),
            )
            p_bcs = None
        else:
            if isinstance(p_dOmega, str):
                p_dOmega = (p_dOmega, )
            p_bcs = BoundaryConditions(
                ('natural', dOmega[p_dOmega], 0.0, 1),
            )
            u_dOmega = tuple(
                i for i in ('left', 'right', 'upper') if not i in p_dOmega
            )
            u_bcs = BoundaryConditions(
                ('essential', dOmega['lower'], (0.0, uLower), 0),
                *([('essential', dOmega[u_dOmega], (0.0, 0.0), 0)] if u_dOmega else []),
            )
        flow_bcs = (u_bcs, p_bcs)
    # constitutive relations
    gamma = Constant(Omega, gamma, 'gamma')
    delta = Constant(Omega, delta, 'delta')
    dispersion_solutal = lambda phi: Di * phi
    dispersion_thermal = lambda phi: (Di/Le) * phi
    density = lambda c, theta: Bu * (c - gamma * theta)
    reaction = lambda theta, s: -Ki * s * (1 + delta * theta)
    source = lambda theta, s: Ki * s * (1 + delta * theta)

    if diagnostic:
        fluxes_thermal = [('q', 0.1 * Ly, Lx), *fluxes_thermal]
        fluxes_solutal = [('f', 0.1 * Ly, Lx), *fluxes_solutal] 

    auxiliary = [
        Ra, Da, Le, Di, Bu, Ki, 
        gamma, delta,
        ('X', X), ('Lx', Lx), ('Ly', Ly), 
        ('l', l), ('Ll', Ll),
        ('In', In), ('xInflow', xInflow),
        ('sr', sr), ('cr', cr), 
    ]

    return dns_generic(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # physical
        epsilon=epsilon,
        # initial conditions
        c_ics=c_ics,
        theta_ics=theta_ics,
        s_ics=s_ics, 
        # boundary conditions
        flow_bcs=flow_bcs,
        c_bcs=c_bcs,
        theta_bcs=theta_bcs,
        # constitutive relations
        permeability=permeability,
        density=density,
        reaction=reaction,
        source=source,
        dispersion_solutal=dispersion_solutal,
        dispersion_thermal=dispersion_thermal,
        # timestep
        dt_min=dt_min,
        dt_max=dt_max,
        dt_h=dt_h,
        dt_Cu=dt_Cu,
        dt_Cd=dt_Cd,
        dt_Cr=dt_Cr,
        # time discretization
        D_adv_solutal=D_adv_solutal,
        D_diff_solutal=D_diff_solutal,
        D_reac_solutal=D_reac,
        D_src_solutal=D_src,
        D_adv_thermal=D_adv_thermal,
        D_diff_thermal=D_diff_thermal,
        D_evol=D_evol,
        # stabilization
        c_stabilization=c_stabilization,
        c_limits=c_limits,
        theta_stabilization=theta_stabilization,
        theta_limits=theta_limits,
        s_limits=s_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        theta_petsc=theta_petsc,
        c_petsc=c_petsc,
        s_petsc=s_petsc,
        # optional solvers
        diagnostic=diagnostic,
        fluxes_solutal=fluxes_solutal,
        fluxes_thermal=fluxes_thermal,
        auxiliary=auxiliary,
    )

