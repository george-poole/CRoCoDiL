import numpy as np
from mpi4py import MPI

from lucifex.mesh import rectangle_mesh, mesh_boundary
from lucifex.fem import Function, Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.sim import configure_simulation
from lucifex.solver import OptionsPETSc, OptionsJIT, BoundaryConditions
from lucifex.utils.fenicsx_utils import CellType, limits_corrector
from lucifex.utils.py_utils import FrozenDict

from .generic import dns_generic
from .utils import DEFAULT_JIT_DIR, SCALINGS, heaviside, rectangle_mesh_closure


SYSTEM_D_REFERENCE = FrozenDict(
    Ra=1000.0,
    Da=100.0,
    epsilon=1e-2,
    Le=1.0,
    zeta0=0.9,
    sr=0.2,
    cr=0.0,
    gamma=1.0,
    delta=0,
    lIn=0.1,
    uI=1.0,
    thetaIn=1.0,
    aspect=2.0,
)


@configure_simulation(
    jit=OptionsJIT(DEFAULT_JIT_DIR),
    dir_root="./data",
)
def dns_system_d(
    # mesh
    comm: MPI.Comm | str = 'COMM_WORLD',
    aspect: float = 2.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    scaling: str = 'advective',
    Ra: float = 1e3,
    Da: float = 1e2,
    epsilon: float = 1e-2,
    Le: float = 1.0,
    # initial saturation
    sr: float = 0.2,
    s_ampl: float = 0,
    s_freq: tuple[int, int] = (16, 16),
    s_seed: tuple[int, int] = (1234, 5678),
    # initial concentration
    cr: float = 1.0,
    c_ampl: float = 1e-6,
    c_freq: tuple[int, int] = (16, 16),
    c_seed: tuple[int, int] = (1234, 5678),
    # initial temperature
    theta_ampl: float = 0,
    theta_freq: tuple[int, int] = (16, 16),
    theta_seed: tuple[int, int] = (1234, 5678),
    # constitutive relations
    delta: float = 0,
    gamma: float = 1.0,
    # inflow and outflow boundary conditions
    lIn: float = 0.1,
    uIn: float = 1.0,
    thetaIn: float = 1.0,
    pOut: float | None = None,
    # timestep
    dt_min: float = 0.0,
    dt_max: float = np.inf,
    dt_h: str | float = "hmin",
    courant_adv: float | None = 0.5,
    courant_diff: float | None = 0.5,
    courant_reac: float | None = 0.1,
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
    # fluxes_solutal: Iterable[tuple[str, float | int, float]] = (), 
    # fluxes_thermal: Iterable[tuple[str, float | int, float]] = (), 
):
    # space
    scaling_map = SCALINGS[scaling](Ra, Da)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    LlIn = lIn * X
    Omega = rectangle_mesh((-0.5 * Lx, 0.5 * Lx), Ly, Nx, Ny, cell, 'Omega', comm)
    dOmega = mesh_boundary(
        Omega,
        {
            'upper': lambda x: x[1] - Ly,
            'left': lambda x: x[0] + 0.5 * Lx,
            'right': lambda x: x[0] - 0.5 * Lx,
            'inflow': lambda x: (
                (x[0] < 0.5 * LlIn) &
                (x[0] > -0.5 * LlIn) &
                (np.isclose(x[1], 0))
            ),
            'lower': lambda x: (
                ((x[0] >= 0.5 * LlIn) | (x[0] <= -0.5 * LlIn)) &
                (np.isclose(x[1], 0))
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
    theta_ics = 0
    _theta_limits = (min(0, thetaIn), max(0, thetaIn))
    if theta_ampl:
        theta_ics = SpatialPerturbation(
            theta_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], theta_freq, theta_seed),
            [Lx, Ly],
            theta_ampl,
            limits_corrector(*_theta_limits),
        ) 
    if theta_limits:
        theta_limits = _theta_limits
    # boundary conditions
    c_bcs = BoundaryConditions(
        ('neumann', dOmega['left', 'right', 'upper', 'lower'], 0.0),
        ('dirichlet', dOmega['inflow'], 0.0),
    )
    theta_bcs = BoundaryConditions(
        ('neumann', dOmega['left', 'right', 'upper', 'lower'], 0.0),
        ('dirichlet', dOmega['inflow'], thetaIn),

    )
    if isinstance(flow_petsc, tuple):
        if pOut is not None:
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
        if pOut is None:
            u_bcs = BoundaryConditions(
                ('essential', dOmega['upper', 'lower'], (0.0, 0.0), 0),
                ('essential', dOmega['left', 'right'], (0.0, 0.5 * uIn * lIn), 0),
                ('essential', dOmega['inflow'], (0.0, -uIn), 0),
            )
            p_bcs = None
        else:
            u_bcs = BoundaryConditions(
                ('essential', dOmega['upper', 'lower'], (0.0, 0.0), 0),
                ('essential', dOmega['inflow'], (0.0, -uIn), 0),
            )
            p_bcs = BoundaryConditions(
                ('natural', dOmega['left', 'right'], pOut, 1),
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
        density=density,
        reaction=reaction,
        source=source,
        dispersion_solutal=dispersion_solutal,
        dispersion_thermal=dispersion_thermal,
        # timestep
        dt_min=dt_min,
        dt_max=dt_max,
        dt_h=dt_h,
        courant_adv=courant_adv,
        courant_diff=courant_diff,
        courant_reac=courant_reac,
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
        # fluxes_solutal=fluxes_solutal,
        # fluxes_thermal=fluxes_thermal,
        # auxiliary=auxiliary,
    )

