from types import EllipsisType

from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB
from lucifex.utils import CellType, SpatialPerturbation, cubic_noise
from lucifex.solver import OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation

from .generic import thermosolutal_convection_generic
from .utils import heaviside, rectangle_domain


# FIXME FIXME FIXME
@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
    dir_base="./data",
)
def rayleigh_taylor_rectangle(
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
    | FiniteDifferenceArgwise = (AB(2) @ CN),
    D_diff: FiniteDifference = CN,
    # stabilization
    c_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    c_limits: bool | tuple[float, float] = False,
    # linear algebra
    flow_petsc: OptionsPETSc | None | tuple[OptionsPETSc | None, OptionsPETSc | EllipsisType | None] = (None, ...),
    c_petsc: OptionsPETSc | None = None,
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

    simulation = thermosolutal_convection_generic(
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