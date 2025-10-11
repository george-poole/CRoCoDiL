from types import EllipsisType

from lucifex.fem import LUCiFExConstant as Constant
from lucifex.fdm import FiniteDifference, CN, AB
from lucifex.utils import CellType, SpatialPerturbation, cubic_noise
from lucifex.solver import OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation

from .create import create_simulation, create_rectangle_domain

@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
    dir_base="./data",
)
def thermosolutal_dissolution_2d(
    # mesh
    Lx: float = 2.0,
    Ly: float = 1.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    Ra: float = 1e3,
    Rb: float = 1e3,
    Da: float = 1e2,
    epsilon: float = 1e-2,
    # initial saturation
    sr = 0.1,
    # temperature SpatialPerturbation
    theta_eps: float = 1e-6,
    theta_freq: tuple[int, int] = (8, 8),
    theta_seed: tuple[int, int] = (1234, 5678),
    # constitutive relations
    delta: float = 1.0,
    gamma: float = 1.0,
    # time step
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    k_courant: float = 0.1,
    # time discretization
    D_adv: FiniteDifference | tuple[FiniteDifference, FiniteDifference] = (AB(2), CN),
    D_diff: FiniteDifference = CN,
    D_reac: FiniteDifference | tuple[FiniteDifference, FiniteDifference, FiniteDifference] = (AB(2), CN, CN),
    # stabilization
    c_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    c_limits: bool | tuple[float, float] = False,
    theta_stabilization = None,
    theta_limits = False,
    s_limits = None,
    # linear algebra
    flow_petsc: OptionsPETSc | None | tuple[OptionsPETSc | None, OptionsPETSc | EllipsisType | None] = (None, ...),
    c_petsc: OptionsPETSc | None = None,
    s_petsc: OptionsPETSc | EllipsisType | None = ...,
    # secondary
    secondary: bool = False,   
):
    """
    `c(x,y,t=0) = c₀(y) ` \\
    `θ(x,y,t=0) = 1 - y + N(x,y`) \\
    `s(x,y,t=0) = sr`
    """
    Omega, dOmega = create_rectangle_domain(Lx, Ly, Nx, Ny, cell)

    c_ics = lambda x: x[1]
    theta_ics = SpatialPerturbation(
        lambda x: 1 - x[1] / Ly,
        cubic_noise(['neumann', 'dirichlet'], [Lx, Ly], theta_freq, theta_seed),
        theta_eps,
        [Lx, Ly],
    )   
    s_ics = sr

    gamma = Constant(Omega, gamma)
    density = lambda c, theta: c - gamma * theta
    reaction = lambda s, c, theta: s * (1 + delta * theta - c)

    if c_limits is Ellipsis:
        c_limits = (0, 1 + delta)

    return create_simulation(
        # domain
        Omega=Omega, 
        dOmeg=dOmega, 
        # physical
        Ra=Ra,
        Rb=Rb,
        Da=Da,
        epsilon=epsilon,
        # initial conditions
        c_ics=c_ics,
        theta_ics=theta_ics,
        s_ics=s_ics, 
        # constitutive relations
        density=density,
        reaction=reaction,
        # time step
        dt_max=dt_max,
        cfl_h=cfl_h,
        cfl_courant=cfl_courant,
        k_courant=k_courant,
        # time discretization
        D_adv_solutal=D_adv,
        D_diff_solutal=D_diff,
        D_reac_solutal=D_reac,
        D_adv_thermal=D_adv,
        D_diff_thermal=D_diff,
        D_reac_evol=D_reac,
        # stabilization
        c_stabilization=c_stabilization,
        c_limits=c_limits,
        theta_stabilization=theta_stabilization,
        theta_limits=theta_limits,
        s_limits=s_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        c_petsc=c_petsc,
        s_petsc=s_petsc,
        # optional solvers
        secondary=secondary,
        namespace_extras=[gamma],
    )
