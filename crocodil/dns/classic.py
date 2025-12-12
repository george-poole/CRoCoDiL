from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB2
from lucifex.utils import CellType
from lucifex.solver import BoundaryConditions, OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation

from .generic import dns_generic
from .utils import rectangle_domain, heaviside


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_rayleigh_benard_rectangle(
    # mesh
    Lx: float = 2.0,
    Ly: float = 1.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # temperature SpatialPerturbation
    theta_eps: float = 1e-6,
    theta_freq: tuple[int, int] = (8, 8),
    theta_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    # time discretization
    D_adv: FiniteDifference | FiniteDifferenceArgwise = (AB2, CN),
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
    `θ(x,y,t=0) = 1 - y + N(x,y`)
    """
    Omega, dOmega = rectangle_domain(Lx, Ly, Nx, Ny, cell)
    Ra = Constant(Omega, Ra, 'Ra')
    theta_bcs = BoundaryConditions(
        ("dirichlet", dOmega['lower'], 1.0),
        ("dirichlet", dOmega['upper'], 0.0),
        ('neumann', dOmega['left', 'right'], 0.0)
    )
    # θ₀ = 1 - y + N(x, y)
    theta_ics = SpatialPerturbation(
        lambda x: 1 - x[1] / Ly,
        cubic_noise(['neumann', 'dirichlet'], [Lx, Ly], theta_freq, theta_seed),
        theta_eps,
        [Lx, Ly],
        )   
    
    density = lambda theta: -theta
    dispersion = lambda phi, _: (1/Ra) * phi

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
        dispersion_solutal=None,
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
        secondary=secondary,
    )

    return simulation



@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
    dir_base="./data",
)
def dns_rayleigh_taylor_rectangle(
    # mesh
    Lx: float = 4.0,
    Ly: float = 1.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    Ra: float = 1e3,
    # initial concentration
    cr: float = 0.0,
    h0: float = 0.0,
    H_eps: float | None = None,
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
    Omega, dOmega = rectangle_domain(Lx, Ly, Nx, Ny, cell)
    c_ics = SpatialPerturbation(
        heaviside(lambda x: x[1] - h0, cr, eps=H_eps),
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
        secondary=secondary,
    )

    return simulation
