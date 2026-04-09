from typing import Iterable

import numpy as np
from mpi4py import MPI
from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.fdm.ufl_operators import exp
from lucifex.utils.fenicsx_utils import CellType
from lucifex.solver import OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation
from lucifex.utils.fenicsx_utils import limits_corrector
from lucifex.utils.py_utils import FrozenDict

from .generic import dns_generic
from .utils import CROCODIL_JIT_DIR, heaviside, rectangle_mesh_closure, SCALINGS
from .system_a import SYSTEM_A_REFERENCE


SYSTEM_A_MU_REFERENCE = FrozenDict(
    **SYSTEM_A_REFERENCE,
    gamma=0.5,
)


@configure_simulation(
    jit=OptionsJIT(CROCODIL_JIT_DIR),
)
def dns_system_aMu(
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
    # initial front
    zeta0: float = 0.9,
    zeta0_eps: float | tuple[float, float] | None = None,
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
    # constitutive relations
    viscosity_type: str = 'exponential',
    gamma: float = 0.5,
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
    fluxes: Iterable[tuple[str, float | int, float]] = (),
):
    # space
    scaling_map = SCALINGS[scaling](Ra, Da)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    Lzeta0 = zeta0 * X
    Lzeta0_eps = zeta0_eps * X if zeta0_eps is not None else None
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell, comm=comm)
    # constants
    Di, Ki, Bu = scaling_map(Omega)['Di', 'Ki', 'Bu']
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')
    gamma = Constant(Omega, gamma, 'gamma')
    # initial conditions
    s_ics = heaviside(lambda x: x[1] - Lzeta0, max(0, sr - s_ampl), eps=Lzeta0_eps) 
    if s_ampl:
        s_ics = SpatialPerturbation(
            s_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], s_freq, s_seed),
            [Lx, Ly],
            s_ampl,
            limits_corrector(0, sr),
        )
    c_ics = heaviside(lambda x: x[1] - Lzeta0, max(0, cr - c_ampl), eps=Lzeta0_eps)
    if c_ampl:
        c_ics = SpatialPerturbation(
            c_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed),
            [Lx, Ly],
            c_ampl,
            limits_corrector(0, 1),
        )  
    # constitutive relations
    dispersion = lambda phi: Di * phi
    reaction = lambda s: -Ki * s
    source = lambda s: Ki * s
    density = lambda c: Bu * c
    match viscosity_type:
        case 'linear':
            viscosity = lambda c: 1 + gamma * c
        case 'exponential':
            viscosity = lambda c: exp(gamma * c)
        case _:
            raise ValueError(viscosity_type)
        
    if diagnostic:
        fluxes = [('f', Lzeta0, Lx), *fluxes]

    auxiliary = [
        Ra, Da, Di, Bu, Ki, gamma,
        ('X', X), ('Lx', Lx), ('Ly', Ly), 
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
        # constitutive relations
        density=density,
        viscosity=viscosity,
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
        fluxes_solutal=fluxes,
        diagnostic=diagnostic,
        auxiliary=auxiliary,
    )
