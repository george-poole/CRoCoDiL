import numpy as np
from matplotlib.axes import Axes
from matplotlib.patches import FancyArrowPatch
from mpi4py import MPI
from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.pde.constitutive import permeability_cross_bedded
from lucifex.utils.fenicsx_utils import CellType, limits_corrector
from lucifex.solver import OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation
from lucifex.utils.py_utils import FrozenDict

from .generic import dns_generic
from .utils import DEFAULT_JIT_DIR, heaviside, rectangle_mesh_closure
from .theory import CONVECTION_REACTION_SCALINGS


SYSTEM_C_REFERENCE = FrozenDict(
    Ra=1000.0,
    Da=100.0,
    epsilon=1e-2,
    beta=10.0,
    kappa=1.0,
    vartheta=0.0,
    zeta_min=0.8,
    sr=0.2,
    cr=0.0,
    aspect=2.0,
)


def plot_gravity_arrow(
    ax: Axes,
    start: tuple[float, float],
    length: float,    
    beta: float,
    **kws,
) -> None:
    arrow_start =np.asarray(start)
    arrow_end = arrow_start - length * np.array((np.sin(np.radians(beta)), np.cos(np.radians(beta))))
    _kws = dict(
        arrowstyle='-|>',
        mutation_scale=15,
        color='cyan',
        linewidth=2.0,
    )
    _kws.update(kws)
    arrow = FancyArrowPatch(
        arrow_start,
        arrow_end,
        **_kws,
    )
    ax.add_patch(arrow)


def upper_subdomain_dimensions_initial(
    beta: float,
    zeta0_min: float,
    Lx: float,
    Ly: float,
) -> tuple[float, float]:
    vert = (1 - zeta0_min) * Ly
    horz = vert / np.tan(np.radians(beta))
    if horz > Lx:
        raise ValueError('Not a triangle')
    return horz, vert


def upper_subdomain_hypotenuse_initial(
    beta: float,
    zeta0_min: float,
    Lx: float,
    Ly: float,
) -> tuple[float, float]:
    horz, vert = upper_subdomain_dimensions_initial(beta, zeta0_min, Lx, Ly)
    return np.hypot(horz, vert)


def upper_subdomain_area_initial(
    beta: float,
    zeta0_min: float,
    Lx: float,
    Ly: float,
) -> tuple[float, float]:
    horz, vert = upper_subdomain_dimensions_initial(beta, zeta0_min, Lx, Ly)
    return 0.5 * horz * vert


def mass_dissolved_initial(
    beta: float,
    zeta0_min: float,
    sr: float,
    cr: float,
    Lx: float,
    Ly: float,
)-> float:
    area = upper_subdomain_area_initial(beta, zeta0_min, Lx, Ly)
    return area * (1 - sr) * cr


def mass_capillary_initial(
    epsilon: float,
    beta: float,
    zeta0_min: float,
    sr: float,
    Lx: float,
    Ly: float,
)-> float:
    area = upper_subdomain_area_initial(beta, zeta0_min, Lx, Ly)
    return area * sr / epsilon


@configure_simulation(
    jit=DEFAULT_JIT_DIR,
)
def dns_system_c(
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
    # constitutive relations
    beta: float = 10.0,
    kappa: float = 1.0,
    vartheta: float = 0.0,
    # initial front
    zeta0_min: float = 0.5,
    zeta0_eps: float | tuple[float, float] | None = None,
    # initial conditions
    sr: float = 0.2,
    s_ampl: float = 0,
    s_freq: tuple[int, int] = (16, 16),
    s_seed: tuple[int, int] = (1234, 5678),
    # initial concentration
    cr: float = 1.0,
    c_ampl: float = 1e-6,
    c_freq: tuple[int, int] = (16, 16),
    c_seed: tuple[int, int] = (1234, 5678),
    # timestep
    dt_min: float = 0.0,
    dt_max: float = 0.1,
    dt_h: str | float = "hmin",
    courant_adv: float | None = 0.5,
    courant_diff: float | None = 0.5,
    courant_reac: float | None = 0.1,
    # time discretization
    D_adv: FiniteDifference
    | FiniteDifferenceArgwise = (AB(2) @ CN),
    D_diff: FiniteDifference
    | FiniteDifferenceArgwise = (AB(1) @ CN),
    D_reac: FiniteDifference 
    | FiniteDifferenceArgwise = (AB(1) @ AM(1)),
    D_src: FiniteDifference = AB(1),
    D_evol: FiniteDifference 
    | FiniteDifferenceArgwise = (AM(1) @ AB(1)),
    # stabilization
    c_stabilization: str | tuple[float, float] | None = None,
    c_limits: bool = False,
    s_elem: tuple[str, int] = ('DP', 0),
    s_limits: bool = False,
    phi_elem: tuple[str, int] = ('P', 1),
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'hypre'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    s_petsc: OptionsPETSc | None = None,
    # optional postprocessing
    diagnostic: bool = True,
):
    # space
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra, Da)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    Lzeta0_min = zeta0_min * X
    Lzeta0_eps = zeta0_eps * X if zeta0_eps is not None else None
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell, comm=comm)
    # constants
    Di, Ki, Bu = scaling_map[Omega, 'Di', 'Ki', 'Bu']
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')
    # initial conditions
    zeta0 = lambda x: np.tan(np.radians(beta)) * (Lx - x[0]) + Lzeta0_min * Ly
    s_ics = heaviside(lambda x: x[1] - zeta0(x), max(0, sr - s_ampl), eps=Lzeta0_eps) 
    if s_ampl:
        s_ics = SpatialPerturbation(
            s_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], s_freq, s_seed),
            [Lx, Ly],
            s_ampl,
            limits_corrector(0, sr),
            )   
    c_ics = heaviside(lambda x: x[1] - zeta0(x), max(0, cr - c_ampl), eps=Lzeta0_eps)
    if c_ampl:
        c_ics = SpatialPerturbation(
            c_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed),
            [Lx, Ly],
            c_ampl,
            limits_corrector(0, 1),
        ) 
    # constitutive relations
    vartheta = Constant(Omega, np.radians(beta), name='vartheta')
    kappa = Constant(Omega, kappa, name='kappa')
    if np.isclose(float(kappa), 1):
        permeability = lambda phi: phi**2
    else:
        permeability = lambda phi: permeability_cross_bedded(phi**2, kappa, vartheta)
    dispersion = lambda phi: Di * phi
    reaction = lambda s: -Ki * s
    source = lambda s: Ki * s
    density = lambda c: Bu * c
    egx = Constant(Omega, -np.sin(np.radians(beta)))
    egy = Constant(Omega, -np.cos(np.radians(beta)))

    auxiliary = [
        Ra, Da, Di, Bu, Ki, vartheta, kappa, ('beta', beta),
        ('X', X), ('Lx', Lx), ('Ly', Ly), ('zeta0_min', zeta0_min), 
        ('sr', sr), ('cr', cr),
    ]

    return dns_generic(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        eg=(egx, egy),
        # physical
        epsilon=epsilon,
        # initial conditions
        s_ics=s_ics, 
        c_ics=c_ics, 
        # constitutive relations
        permeability=permeability,
        density=density,
        reaction=reaction,
        source=source,
        dispersion_solutal=dispersion,
        # timestep
        dt_min=dt_min,
        dt_max=dt_max,
        dt_h=dt_h,
        courant_adv=courant_adv,
        courant_diff=courant_diff,
        courant_reac=courant_reac,
        # time discretization
        D_adv_solutal=D_adv,
        D_diff_solutal=D_diff,
        D_reac_solutal=D_reac,
        D_src_solutal=D_src,
        D_evol=D_evol,
        # stabilization
        c_stabilization=c_stabilization,
        c_limits=c_limits,
        s_elem=s_elem,
        s_limits=s_limits,
        phi_elem=phi_elem,
        # linear algebra
        flow_petsc=flow_petsc,
        c_petsc=c_petsc,
        s_petsc=s_petsc,
        # optional solvers
        diagnostic=diagnostic,
        auxiliary=auxiliary,
    )


