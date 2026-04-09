import numpy as np
from mpi4py import MPI
from lucifex.mesh import annulus_sector_mesh, mesh_boundary
from lucifex.fem import Constant, SpatialPerturbation
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.utils.fenicsx_utils import CellType
from lucifex.solver import OptionsPETSc, OptionsJIT, BoundaryConditions
from lucifex.sim import configure_simulation
from lucifex.utils.fenicsx_utils import limits_corrector
from lucifex.utils.py_utils import FrozenDict

from .generic import dns_generic
from .utils import CROCODIL_JIT_DIR, heaviside, SCALINGS


SYSTEM_B_REFERENCE = FrozenDict(
    Ra=1000.0,
    Da=100.0,
    epsilon=1e-2,
    zeta0=0.8,
    sr=0.2,
    cr=0.0,
    aspect=2.0,
)
"""
Reference parameters `Ra, Da, epsilon, zeta0, sr, cr` and `aspect`
governing the physical (as opposed to numerical) behaviour of system A.
"""


@configure_simulation(
    jit=OptionsJIT(CROCODIL_JIT_DIR),
)
def dns_system_b(
    # mesh
    comm: MPI.Comm | str = 'COMM_WORLD',
    aspect: float = 0.5,
    Nr: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    scaling: str = 'advective',
    Ra: float = 1e3,
    Da: float = 1e2,
    epsilon: float = 1e-2,
    Pe: float = 1e3,
    # initial front
    zeta0: float = 0.8,
    zeta0_eps: float | tuple[float, float] | None = None,
    # initial saturation
    sr: float = 0.2,
    s_ampl: float = 0,
    s_freq: int = 8,
    # initial concentration
    cr: float = 1.0,
    c_ampl: float = 1e-6,
    c_freq: int = 8,
    # timestep
    dt_min: float = 0.0,
    dt_max: float = np.inf,
    dt_h: str | float = "hmin",
    dt_Cu: float | None = 0.5, # TODO rename dt_Cu
    dt_Cd: float | None = 0.5,
    dt_Cr: float | None = 0.1,
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
    scaling_map = SCALINGS[scaling](Ra, Da)
    X = scaling_map['X']
    Router = X
    Rinner = aspect * Router
    Lzeta0 = zeta0 * X
    Lzeta0_eps = zeta0_eps * X if zeta0_eps is not None else None
    dr = (Router - Rinner) / Nr
    Omega = annulus_sector_mesh(dr, cell, 'Omega', comm)(Rinner, Router, 180)
    r2 = lambda x: x[0]**2 + x[1]**2
    dOmega = mesh_boundary(
        Omega, 
        {
            "inner": lambda x: r2(x) - Rinner**2,
            "outer": lambda x: r2(x) - Router**2,
            "left": lambda x: np.logical_and(np.isclose(x[1], 0.0), x[0] < 0),
            "right": lambda x: np.logical_and(np.isclose(x[1], 0.0), x[0] > 0),
        },
    )
    # constants
    uIn = (Pe / Ra) * scaling_map['Bu'] # TODO check this
    Di, Ki, Bu = scaling_map(Omega)['Di', 'Ki', 'Bu']
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')
    # initial conditions
    s_ics = heaviside(lambda x: x[1] - Lzeta0, max(0, sr - s_ampl), eps=Lzeta0_eps) 
    if s_ampl:
        radial_noise = lambda x: (s_ampl 
            * np.cos(s_freq * np.pi * (np.sqrt(r2(x) - Rinner) / Router))
        )
        s_ics = SpatialPerturbation(
            s_ics,
            radial_noise,
            Omega.geometry.x,
            s_ampl,
            limits_corrector(0, sr),
        )
    c_ics = heaviside(lambda x: x[1] - Lzeta0, max(0, cr - c_ampl), eps=Lzeta0_eps)
    if c_ampl:
        radial_noise = lambda x: (c_ampl 
            * np.cos(c_freq * np.pi * (np.sqrt(r2(x) - Rinner) / Router))
        )
        c_ics = SpatialPerturbation(
            c_ics,
            radial_noise,
            Omega.geometry.x,
            c_ampl,
            limits_corrector(0, 1),
        )  
    # boundary conditions
    if isinstance(flow_petsc, tuple):
        psi_bcs = BoundaryConditions(
            ('dirichlet', dOmega['outer'], 0.0),
            ('dirichlet', dOmega['left'], lambda x: -uIn * (x[0] + Router)),
            ('dirichlet', dOmega['inner'], -uIn * (Router - Rinner)),
            ('dirichlet', dOmega['right'], lambda x: uIn * (x[0] - Router)),
        )
        flow_bcs = psi_bcs
    else:
        u_bcs = BoundaryConditions(
            ('dirichlet', dOmega['outer', 'inner'], (0.0, 0.0), 0),
            ('dirichlet', dOmega['left'], (0.0, uIn), 0),
            ('dirichlet', dOmega['left'], (0.0, -uIn), 0),
        )
        p_bcs = None
        flow_bcs = (u_bcs, p_bcs)
    c_bcs = BoundaryConditions(
        ('dirichlet', dOmega.union, 0.0),
    )
    # constitutive relations
    dispersion = lambda phi: Di * phi
    reaction = lambda s: -Ki * s
    source = lambda s: Ki * s
    density = lambda c: Bu * c

    auxiliary = [
        Ra, Da, Di, Bu, Ki, 
        ('X', X), ('Rinner', Rinner), ('Router', Router), 
        ('zeta0', zeta0), ('Lzeta0', Lzeta0),
        ('sr', sr), ('cr', cr),
    ]

    return dns_generic(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # physical
        epsilon=epsilon,
        # initial conditions
        s_ics=s_ics, 
        c_ics=c_ics,
        # boundary conditions
        c_bcs=c_bcs,
        flow_bcs=flow_bcs,
        # constitutive relations
        density=density,
        reaction=reaction,
        source=source,
        dispersion_solutal=dispersion,
        # timestep
        dt_min=dt_min,
        dt_max=dt_max,
        dt_h=dt_h,
        dt_Cu=dt_Cu,
        dt_Cd=dt_Cd,
        dt_Cr=dt_Cr,
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