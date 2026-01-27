import scipy.special as sp

from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB2
from lucifex.utils import CellType
from lucifex.solver import BoundaryConditions, OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation

from .generic import dns_generic
from .utils import CONVECTION_REACTION_SCALINGS, rectangle_mesh_closure, heaviside


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_rayleigh_benard_2s(
    # mesh
    aspect: float = 2.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    scaling: str = 'advective',
    Ra: float = 1e3,
    # initial conditions
    theta_eps: float = 1e-6,
    theta_freq: tuple[int, int] = (8, 8),
    theta_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_max: float = 0.1,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.5,
    # time discretization
    D_adv: FiniteDifference | FiniteDifferenceArgwise = (AB2 @ CN),
    D_diff: FiniteDifference = CN,
    # stabilization
    theta_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    theta_limits: bool | tuple[float, float] = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    # secondary
    secondary: bool = False,
):
    """
    `c(x,y,t=0) = ... y ... ` \\
    `θ(x,y,t=0) = 1 - y/Ly + N(x,y`)
    """
    scaling_mapping = CONVECTION_REACTION_SCALINGS[scaling].mapping(Ra)

    Xl = scaling_mapping('Xl')
    Lx = aspect * Xl
    Ly = 1.0 * Xl
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)
    
    Ra = Constant(Omega, Ra, 'Ra')
    theta_bcs = BoundaryConditions(
        ("dirichlet", dOmega['lower'], 1.0),
        ("dirichlet", dOmega['upper'], 0.0),
        ('neumann', dOmega['left', 'right'], 0.0)
    )
    theta_ics = SpatialPerturbation(
        lambda x: 1 - x[1] / Ly,
        cubic_noise(['neumann', 'dirichlet'], [Lx, Ly], theta_freq, theta_seed),
        [Lx, Ly],
        theta_eps,
        )   
    
    density = lambda theta: -theta
    dispersion = lambda phi: (1/Ra) * phi

    simulation = dns_generic(
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
        c_petsc=c_petsc,
        # optional solvers
        diagnostic=secondary,
    )

    return simulation


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_rayleigh_benard_1s(
    # mesh
    aspect: float = 2.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    scaling: str ='advective',
    Ra: float = 1e3,
    # initial conditions
    erf_eps: float = 1e-2,
    theta_eps: float = 1e-6,
    theta_freq: tuple[int, int] = (8, 8),
    theta_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_max: float = 0.1,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.5,
    # time discretization
    D_adv: FiniteDifference | FiniteDifferenceArgwise = (AB2 @ CN),
    D_diff: FiniteDifference = CN,
    # stabilization
    theta_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    theta_limits: bool | tuple[float, float] = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    theta_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    # secondary
    secondary: bool = False,
):
    scaling_mapping = CONVECTION_REACTION_SCALINGS[scaling].mapping(Ra)
    Xl = scaling_mapping('Xl')
    Lx = aspect * Xl
    Ly = 1.0 * Xl
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)

    Ra = Constant(Omega, Ra, 'Ra')
    Pe, Bu = scaling_mapping('Pe', 'Bu', mesh=Omega)

    theta_bcs = BoundaryConditions(
        ("dirichlet", dOmega['upper'], 1.0),
        ('neumann', dOmega['lower', 'upper', 'right'], 0.0)
    )
    theta_ics = SpatialPerturbation(
        lambda x: 1 + sp.erf((x[1] - Ly) / (Ly * erf_eps)),
        cubic_noise(['neumann', ('neumann', 'dirichlet')], [Lx, Ly], theta_freq, theta_seed),
        [Lx, Ly],
        theta_eps,  
        )

    density = lambda theta: Bu * theta
    dispersion = lambda phi: (1/Pe) * phi

    simulation = dns_generic(
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
        diagnostic=secondary,
    )

    return simulation


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_rayleigh_taylor(
    # mesh
    aspect: float = 2.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    scaling: str ='advective',
    Ra: float = 1e3,
    # initial conditions
    cr: float = 0.0,
    h0: float = 0.0,
    h0_eps: float | None = None,
    c_eps: float = 1e-6,
    c_freq: tuple[int, int] = (8, 8),
    c_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    # time discretization
    D_adv: FiniteDifference 
    | FiniteDifferenceArgwise = (AB2 @ CN),
    D_diff: FiniteDifference = CN,
    # stabilization
    c_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    c_limits: bool | tuple[float, float] = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    # secondary
    secondary: bool = False,
):
    """
    `c₀(x,y) = cr · H(y - h₀) + N(x, y)`
    """
    scaling_mapping = CONVECTION_REACTION_SCALINGS[scaling].mapping(Ra)
    Xl = scaling_mapping('Xl')
    Lx = aspect * Xl
    Ly = 1.0 * Xl
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)

    c_ics = SpatialPerturbation(
        heaviside(lambda x: x[1] - h0, cr, eps=h0_eps),
        cubic_noise(['dirichlet', 'dirichlet'], [Lx, Ly], c_freq, c_seed),
        c_eps,
        [Lx, Ly],
        )   
    
    density = lambda c: c

    simulation = dns_generic(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # physical
        Ra=Ra,
        # initial conditions
        c_ics=c_ics,
        # constitutive relations
        density=density,
        # time step
        dt_max=dt_max,
        cfl_h=cfl_h,
        cfl_courant=cfl_courant,
        # time discretization
        D_adv_thermal=D_adv,
        D_diff_thermal=D_diff,
        # stabilization
        c_stabilization=c_stabilization,
        c_limits=c_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        c_petsc=c_petsc,
        # optional solvers
        diagnostic=secondary,
    )

    return simulation
