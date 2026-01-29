import scipy.special as sp

from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB
from lucifex.utils import CellType
from lucifex.solver import BoundaryConditions, OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation

from .generic import dns_generic
from .utils import CONVECTION_REACTION_SCALINGS, rectangle_mesh_closure, heaviside


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_darcy_rayleigh_benard(
    # domain
    aspect: float = 2.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    scaling: str = 'advective',
    Ra: float = 1e3,
    # initial conditions
    theta_ampl: float = 1e-6,
    theta_freq: tuple[int, int] = (8, 8),
    theta_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    # time discretization
    D_adv: FiniteDifference | FiniteDifferenceArgwise = (AB(2) @ CN),
    D_diff: FiniteDifference | FiniteDifferenceArgwise = CN,
    # stabilization
    theta_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    theta_limits: bool | tuple[float, float] = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    theta_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    # optional post-processing
    diagnostic: bool = False,
    fluxes = (),
):
    # space
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)
    # constants
    Di, Bu = scaling_map[Omega, 'Di', 'Bu']
    Ra = Constant(Omega, Ra, 'Ra')
    # initial and boundary conditions
    theta_ics = SpatialPerturbation(
        lambda x: 1 - x[1] / Ly,
        cubic_noise(['neumann', 'dirichlet'], [Lx, Ly], theta_freq, theta_seed),
        [Lx, Ly],
        theta_ampl,
    )
    theta_bcs = BoundaryConditions(
        ("dirichlet", dOmega['lower'], 1.0),
        ("dirichlet", dOmega['upper'], 0.0),
        ('neumann', dOmega['left', 'right'], 0.0)
    )
    # constitutive
    dispersion = lambda phi: Di * phi
    density = lambda theta: -Bu * theta

    if diagnostic:
        fluxes = [('q', int(0.5 * Ny), Lx), *fluxes]

    return dns_generic(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # initial conditions
        theta_ics=theta_ics,
        # boundary conditions
        theta_bcs=theta_bcs,
        # constitutive relations
        density=density,
        dispersion_thermal=dispersion,
        # time step
        dt_max=dt_max,
        cfl_h=cfl_h,
        cfl_courant=cfl_courant,
        # time discretization
        D_adv_thermal=D_adv,
        D_diff_thermal=D_diff,
        # stabilization
        theta_stabilization=theta_stabilization,
        theta_limits=theta_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        theta_petsc=theta_petsc,
        # optional solvers
        diagnostic=diagnostic,
        fluxes_thermal=fluxes,
        namespace=[Ra, Di, Bu, ('X', X)],
    )


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_darcy_evolving(
    # domain
    aspect: float = 2.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    scaling: str ='advective',
    Ra: float = 1e3,
    # initial conditions
    theta_ampl: float = 1e-6,
    theta_freq: tuple[int, int] = (8, 8),
    theta_seed: tuple[int, int] = (1234, 5678),
    theta_eps: float = 1e-2,
    # time step
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    # time discretization
    D_adv: FiniteDifference | FiniteDifferenceArgwise = (AB(2) @ CN),
    D_diff: FiniteDifference = CN,
    # stabilization
    theta_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    theta_limits: bool | tuple[float, float] = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    theta_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    # optional post-processing
    diagnostic: bool = False,
    fluxes = (),
):
    # space
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)
    # constants
    Di, Bu = scaling_map[Omega, 'Di', 'Bu']
    Ra = Constant(Omega, Ra, 'Ra')
    # initial and boundary conditions
    theta_ics = SpatialPerturbation(
        lambda x: 1 + sp.erf((x[1] - Ly) / (Ly * theta_eps)),
        cubic_noise(['neumann', ('neumann', 'dirichlet')], [Lx, Ly], theta_freq, theta_seed),
        [Lx, Ly],
        theta_ampl,  
    )
    theta_bcs = BoundaryConditions(
        ("dirichlet", dOmega['upper'], 1.0),
        ('neumann', dOmega['lower', 'upper', 'right'], 0.0)
    )
    # constitutive relations
    density = lambda theta: Bu * theta
    dispersion = lambda phi: Di * phi

    if diagnostic:
        fluxes = [('q', -2, Lx), *fluxes]

    return dns_generic(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # initial conditions
        theta_ics=theta_ics,
        # boundary conditions
        theta_bcs=theta_bcs,
        # constitutive relations
        density=density,
        dispersion_thermal=dispersion,
        # time step
        dt_max=dt_max,
        cfl_h=cfl_h,
        cfl_courant=cfl_courant,
        # time discretization
        D_adv_thermal=D_adv,
        D_diff_thermal=D_diff,
        # stabilization
        theta_stabilization=theta_stabilization,
        theta_limits=theta_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        theta_petsc=theta_petsc,
        # optional solvers
        diagnostic=diagnostic,
        fluxes_thermal=fluxes,
        namespace=[Ra, Di, Bu, ('X', X)],
    )


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_darcy_rayleigh_taylor(
    # domain
    aspect: float = 2.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    scaling: str ='advective',
    Ra: float = 1e3,
    # initial conditions
    cr: float = 1.0,
    h0: float = 0.5,
    h0_eps: float | None = None,
    c_ampl: float = 1e-6,
    c_freq: tuple[int, int] = (8, 8),
    c_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    # time discretization
    D_adv: FiniteDifference 
    | FiniteDifferenceArgwise = (AB(2) @ CN),
    D_diff: FiniteDifference = CN,
    # stabilization
    c_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    c_limits: bool | tuple[float, float] = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    # optional post-processing
    diagnostic: bool = False,
    fluxes = (),
):
    # space
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)
    # constants
    Di, Bu = scaling_map[Omega, 'Di', 'Bu']
    Ra = Constant(Omega, Ra, 'Ra')
    # initial conditions
    c_ics = SpatialPerturbation(
        heaviside(lambda x: x[1] - h0 * X, cr, eps=h0_eps * X if h0_eps is not None else None),
        cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed),
        [Lx, Ly],
        c_ampl,
    )
    c_bcs = BoundaryConditions(
        ("neumann", dOmega.union, 0.0),
    )   
    # constitutive
    density = lambda c: Bu * c
    dispersion = lambda phi: Di * phi

    if diagnostic:
        fluxes = [('f', int(0.5 * Ny), Lx), *fluxes]

    return dns_generic(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # initial conditions
        c_ics=c_ics,
        # boundary conditions
        c_bcs=c_bcs,
        # constitutive relations
        density=density,
        dispersion_solutal=dispersion,
        # time step
        dt_max=dt_max,
        cfl_h=cfl_h,
        cfl_courant=cfl_courant,
        # time discretization
        D_adv_solutal=D_adv,
        D_diff_solutal=D_diff,
        # stabilization
        c_stabilization=c_stabilization,
        c_limits=c_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        c_petsc=c_petsc,
        # optional solvers
        diagnostic=diagnostic,
        fluxes_solutal=fluxes,
        namespace=[Ra, Di, Bu, ('X', X)],
    )


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_darcy_thermosolutal(
    # domain
    aspect: float = 2.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    scaling: str ='advective',
    Ra: float = 1e3,
    Le: float = 1.0,
    gamma: float = 1.0,
    # initial conditions
    c_ampl: float = 1e-3,
    c_freq: tuple[int, int] = (8, 8),
    c_seed: tuple[int, int] = (1234, 5678),
    theta_ampl: float = 1e-3,
    theta_freq: tuple[int, int] = (8, 8),
    theta_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    # time discretization
    D_adv_solutal: FiniteDifference 
    | FiniteDifferenceArgwise = (AB(2) @ CN),
    D_diff_solutal: FiniteDifference
    | FiniteDifferenceArgwise = CN,
    D_adv_thermal: FiniteDifference 
    | FiniteDifferenceArgwise = (AB(2) @ CN),
    D_diff_thermal: FiniteDifference
    | FiniteDifferenceArgwise = CN,
    # stabilization
    c_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    c_limits: bool | tuple[float, float] = False,
    theta_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    theta_limits: bool | tuple[float, float] = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    theta_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    # optional post-processing
    diagnostic: bool = False,
):
    # space
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)
    # constants
    Di, Bu = scaling_map[Omega, 'Di', 'Bu']
    Ra = Constant(Omega, Ra, 'Ra')
    Le = Constant(Omega, Le, 'Le')
    gamma = Constant(Omega, gamma, 'gamma')
    # initial conditions
    c_ics = SpatialPerturbation(
        lambda x: x[1] / Ly,
        cubic_noise(['neumann', 'dirichlet'], [Lx, Ly], c_freq, c_seed),
        [Lx, Ly],
        c_ampl,
    )
    theta_ics = SpatialPerturbation(
        lambda x: 1 - x[1] / Ly,
        cubic_noise(['neumann', 'dirichlet'], [Lx, Ly], theta_freq, theta_seed),
        [Lx, Ly],
        theta_ampl,
    )     
    # boundary conditions
    c_bcs = BoundaryConditions(
        ("dirichlet", dOmega['lower'], 0.0),
        ("dirichlet", dOmega['upper'], 1.0),
        ('neumann', dOmega['left', 'right'], 0.0)
    ) 
    theta_bcs = BoundaryConditions(
        ("dirichlet", dOmega['lower'], 1.0),
        ("dirichlet", dOmega['upper'], 0.0),
        ('neumann', dOmega['left', 'right'], 0.0)
    ) 
    # constitutive relations
    density = lambda c, theta: Bu * (c - gamma * theta)
    dispersion_solutal = lambda phi: Di * phi
    dispersion_thermal = lambda phi: (Di / Le) * phi

    return dns_generic(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # initial conditions
        c_ics=c_ics,
        theta_ics=theta_ics,
        # boundary conditions
        c_bcs=c_bcs,
        theta_bcs=theta_bcs,
        # constitutive relations
        density=density,
        dispersion_solutal=dispersion_solutal,
        dispersion_thermal=dispersion_thermal,
        # time step
        dt_max=dt_max,
        cfl_h=cfl_h,
        cfl_courant=cfl_courant,
        # time discretization
        D_adv_solutal=D_adv_solutal,
        D_diff_solutal=D_diff_solutal,
        D_adv_thermal=D_adv_thermal,
        D_diff_thermal=D_diff_thermal,
        # stabilization
        c_stabilization=c_stabilization,
        c_limits=c_limits,
        theta_stabilization=theta_stabilization,
        theta_limits=theta_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        c_petsc=c_petsc,
        theta_petsc=theta_petsc,
        # optional solvers
        diagnostic=diagnostic,
        namespace=[Ra, Le, Di, Bu, gamma, ('X', X)],
    )