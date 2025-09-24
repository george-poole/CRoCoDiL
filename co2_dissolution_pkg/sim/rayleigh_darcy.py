from types import EllipsisType

from lucifex.fdm import ConstantSeries, FiniteDifference, CN, AB2
from lucifex.utils import CellType, Perturbation, cubic_noise
from lucifex.solver import BoundaryConditions, OptionsPETSc, OptionsJIT, dS_solver
from lucifex.sim import configure_simulation

from ..math.secondary import flux

from .domain import create_rectangle_domain
from .abstract import abstract_simulation


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def rayleigh_darcy_2d(
    # mesh
    Lx: float = 4.0,
    Ly: float = 1.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    Rb: float = 1e3,
    # temperature perturbation
    theta_eps: float = 1e-6,
    theta_freq: tuple[int, int] = (8, 8),
    theta_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    # time discretization
    D_adv: FiniteDifference | tuple[FiniteDifference, FiniteDifference] = (AB2, CN),
    D_diff: FiniteDifference = CN,
    # stabilization
    theta_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    theta_limits: bool | tuple[float, float] = False,
    # linear algebra
    flow_petsc: OptionsPETSc | None | tuple[OptionsPETSc | None, OptionsPETSc | EllipsisType | None] = (None, ...),
    c_petsc: OptionsPETSc | None = None,
    # secondary
    secondary: bool = False,
):
    """
    `c(x,y,t=0) = ... y ... ` \\
    `θ(x,y,t=0) = 1 - y + N(x,y`)
    """
    Omega, dOmega = create_rectangle_domain(Lx, Ly, Nx, Ny, cell)
    theta_bcs = BoundaryConditions(
        ("dirichlet", dOmega['lower'], 1.0),
        ("dirichlet", dOmega['upper'], 0.0),
        ('neumann', dOmega['left', 'right'], 0.0)
    )
    # θ₀ = 1 - y + N(x, y)
    theta_ics = Perturbation(
        lambda x: 1 - x[1],
        cubic_noise(['neumann', 'dirichlet'], [Lx, Ly], theta_freq, theta_seed),
        theta_eps,
        [Lx, Ly],
        )   
    
    density = lambda theta: -theta

    simulation = abstract_simulation(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # physical
        Rb=Rb,
        # initial conditions
        theta_ics=theta_ics,
        # boundary conditions
        theta_bcs=theta_bcs,
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
        theta_stabilization=theta_stabilization,
        theta_limits=theta_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        c_petsc=c_petsc,
        # optional solvers
        secondary=secondary,
    )

    theta, u, g = simulation['theta', 'u', 'g']
    j = ConstantSeries(Omega, "j", shape=(2, ))
    simulation.solvers.append(
        dS_solver(j, flux, lambda x: x[1] - Ly / 2, facet_side="+")(
            theta[0], u[0], g[0], Rb,
        )
    )

    return simulation

