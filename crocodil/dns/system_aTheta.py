from typing import Iterable

import numpy as np
from mpi4py import MPI
from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.solver import OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation
from lucifex.utils.fenicsx_utils import CellType, limits_corrector
from lucifex.utils.py_utils import FrozenDict

from .generic import dns_generic
from .utils import CROCODIL_JIT_DIR, SCALINGS, heaviside, rectangle_mesh_closure
from .system_a import SYSTEM_A_REFERENCE


SYSTEM_A_THETA_REFERENCE = FrozenDict(
    **SYSTEM_A_REFERENCE,
    Le=1.0,
    gamma=1.0,
    delta=0,
)


def thermal_rayleigh(
    Ra: float,
    Le: float,
    gamma: float,
) -> float:
    """
    `(Ra)thermal = (γ / Le)·(Ra)solutal` given density `ρ(c,θ) = c - γθ`
    """
    return gamma * Ra / Le 


@configure_simulation(
    jit=OptionsJIT(CROCODIL_JIT_DIR),
    dir_root="./data",
)
def dns_system_aTheta(
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
    # initial temperature
    theta_buoy: bool = True, 
    theta_ampl: float = 1e-6,
    theta_freq: tuple[int, int] = (16, 16),
    theta_seed: tuple[int, int] = (1234, 5678),
    # constitutive relations
    delta: float = 0,
    gamma: float = 1.0,
    # timestep
    dt_min: float = 0.0,
    dt_max: float = np.inf,
    dt_h: str | float = "hmin",
    dt_Cu: float | None = 0.5,
    dt_Cd: float | None = 0.5,
    dt_Cr: float | None = 0.1,
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
    fluxes_solutal: Iterable[tuple[str, float | int, float]] = (), 
    fluxes_thermal: Iterable[tuple[str, float | int, float]] = (), 
):
    """
    `Ω = [0, A·X] × [0, X]` \\
    `𝜑∂s/∂t = -εKi s(1 + δθ)(1 - c)` \\
    `ϕ∂c/∂t + 𝐮·∇c =  Di ∇·(ϕ∇c) + Ki s(1 + δθ)(1 - c)` \\
    `ϕ∂θ/∂t + 𝐮·∇θ =  LeDi ∇·(ϕ∇θ) \\
    `∇⋅𝐮 = 0` \\
    `𝐮 = -(∇p + Bu(c - γθ)𝐞ʸ)` \\

    `𝐧⋅∇c = 0` on `∂Ω` \\
    `𝐧⋅𝐮 = 0` on `∂Ω`
    """
    # space
    scaling_map = SCALINGS[scaling](Ra, Da)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    Lzeta0 = zeta0 * X
    Lzeta0_eps = zeta0_eps * X if zeta0_eps is not None else None
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell, comm=comm)
    # constants
    Di, Bu, Ki = scaling_map[Omega, 'Di', 'Bu', 'Ki']
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')
    Le = Constant(Omega, Le, 'Le')
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
    if theta_buoy:
        thetaPlus = 0.0
        thetaMinus = 1.0
    else:
        thetaPlus = 1.0 - theta_ampl
        thetaMinus = 0.0
    theta_ics = heaviside(lambda x: x[1] - Lzeta0, thetaPlus, thetaMinus, eps=Lzeta0_eps)
    if theta_ampl:
        theta_ics = SpatialPerturbation(
            theta_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], theta_freq, theta_seed),
            [Lx, Ly],
            theta_ampl,
            limits_corrector(0, 1),
        )
    # constitutive relations
    gamma = Constant(Omega, gamma, 'gamma')
    delta = Constant(Omega, delta, 'delta')
    dispersion_solutal = lambda phi: Di * phi
    dispersion_thermal = lambda phi: (Di/Le) * phi
    density = lambda c, theta: Bu * (c - gamma * theta)
    reaction = lambda theta, s: -Ki * s * (1 + delta * theta)
    source = lambda theta, s: Ki * s * (1 + delta * theta)

    if diagnostic:
        fluxes_thermal = [('q', Lzeta0, Lx), *fluxes_thermal]
        fluxes_solutal = [('f', Lzeta0, Lx), *fluxes_solutal] 

    auxiliary = [
        Ra, Da, Le, Di, Bu, Ki,
        gamma, delta,
        ('X', X), ('Lx', Lx), ('Ly', Ly), 
        ('sr', sr), ('cr', cr), ('zeta0', zeta0),
    ]

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
        dt_Cu=dt_Cu,
        dt_Cd=dt_Cd,
        dt_Cr=dt_Cr,
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
        fluxes_solutal=fluxes_solutal,
        fluxes_thermal=fluxes_thermal,
        auxiliary=auxiliary,
    )
