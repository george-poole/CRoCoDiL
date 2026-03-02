from typing import Iterable, Callable

import numpy as np
from mpi4py import MPI
from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.utils.fenicsx_utils import CellType, limits_corrector
from lucifex.solver import BoundaryConditions, OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation
from lucifex.utils.py_utils import FrozenDict

from .generic import dns_generic
from .utils import rectangle_mesh_closure
from .theory import CONVECTION_REACTION_SCALINGS


SYSTEM_B_REFERENCE = FrozenDict(
    Ra=1000.0,
    Da=100.0,
    epsilon=1e-2,
    Le=1.0,
    sr=0.2,
    gamma=1.0,
    delta=1e-1,
    aspect=2.0,
)


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
    dir_root="./data",
)
def dns_system_b(
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
    # initial saturation
    sr: float = 0.2,
    s_ampl: float = 0,
    s_freq: tuple[int, int] = (16, 16),
    s_seed: tuple[int, int] = (1234, 5678),
    # initial concentration
    c_ics: float | Callable[[np.ndarray], np.ndarray] = 0.0,
    c_ampl: float = 0,
    c_freq: tuple[int, int] = (16, 16),
    c_seed: tuple[int, int] = (1234, 5678),
    # initial temperature
    theta_ampl: float = 1e-6,
    theta_freq: tuple[int, int] = (16, 16),
    theta_seed: tuple[int, int] = (1234, 5678),
    # constitutive relations
    delta: float = 0.1,
    gamma: float = 1.0,
    # timestep
    dt_min: float = 0.0,
    dt_max: float = 0.1,
    dt_h: str | float = "hmin",
    courant_adv: float | None = 0.5,
    courant_diff: float | None = 0.5,
    courant_reac: float | None = 0.1,
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
    c_limits: bool | tuple[float, float] = False,
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
    `𝜑∂s/∂t = -εKi s(1 + δθ - c)` \\
    `ϕ∂c/∂t + 𝐮·∇c =  Di ∇·(ϕ∇c) + Ki s(1 + δθ - c)` \\
    `ϕ∂θ/∂t + 𝐮·∇θ =  Di/Le ∇·(ϕ∇θ) \\
    `∇⋅𝐮 = 0` \\
    `𝐮 = -(∇p + Bu(c -γθ)𝐞ʸ)` \\

    `s₀(𝐱) = sᵣ + N(𝐱)` \\
    `c₀(𝐱) = 1 + δ(1 - θ) + (δ / √Λ) · sinh(√Λ · (y − X/2)) / cosh(√Λ X/2) + N(𝐱)` \\
    `θ₀(𝐱) = 1 - y + N(𝐱)` \\
    `𝐧⋅∇c = 0` on `∂Ω` \\
    `𝐧⋅𝐮 = 0` on `∂Ω`
    """
    # space
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra, Da)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell, comm=comm)
    # constants
    Di, Bu, Ki = scaling_map[Omega, 'Di', 'Bu', 'Ki']
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')
    Le = Constant(Omega, Le, 'Le')
    # initial conditions
    s_ics = sr
    if s_ampl:
        s_ics = SpatialPerturbation(
            s_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], s_freq, s_seed),
            [Lx, Ly],
            s_ampl,
            limits_corrector(0, sr),
        )
    # sqrtLmbda = np.sqrt(float(Ki) * sr / (float(Di) * (1 - sr)))
    # c_ics = lambda x: (
    #     1 + delta * (1 - x[1]) 
    #     + (delta / sqrtLmbda) * np.sinh(sqrtLmbda * (x[1] - 0.5 * Ly)) / np.cosh(0.5 * sqrtLmbda * Ly)
    # )
    if c_ics is None:
        c_ics = 1.0
    if c_ampl:
        c_ics = SpatialPerturbation(
            c_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed),
            [Lx, Ly],
            c_ampl,
            limits_corrector(0, 1 + delta),
        ) 
    if c_limits is True:
        c_limits = (0, 1 + delta)
    theta_ics = lambda x: 1 - x[1] / Ly,
    if theta_ampl:
        SpatialPerturbation(
            theta_ics,
            cubic_noise(['neumann', 'dirichlet'], [Lx, Ly], theta_freq, theta_seed),
            [Lx, Ly],
            theta_ampl,
            limits_corrector(0, 1),
        )   
    # boundary conditions
    theta_bcs = BoundaryConditions(
        ("dirichlet", dOmega['lower'], 1.0),
        ("dirichlet", dOmega['upper'], 0.0),
        ('neumann', dOmega['left', 'right'], 0.0)
    )
    # constitutive relations
    gamma = Constant(Omega, gamma, 'gamma')
    delta = Constant(Omega, delta, 'delta')
    dispersion_solutal = lambda phi: Di * phi
    dispersion_thermal = lambda phi: (Di/Le) * phi
    density = lambda c, theta: Bu * (c - gamma * theta)
    # reaction = lambda _, s: -Ki * s
    # source = lambda theta, s: Ki * s * (1 + delta * theta)
    reaction = lambda theta, s: -Ki * s * (1 + delta * theta)
    source = lambda theta, s: Ki * s * (1 + delta * theta)

    if diagnostic:
        fluxes_thermal = [('q', 0.5 * Ly, Lx), *fluxes_thermal]
        fluxes_solutal = [('f', 0.5 * Ly, Lx), *fluxes_solutal] 

    namespace=[
        Ra, Da, Le, Di, Bu, Ki, 
        ('X', X), ('Lx', Lx), ('Ly', Ly), 
        ('sr', sr),
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
        # boundary conditions
        theta_bcs=theta_bcs,
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
        courant_adv=courant_adv,
        courant_diff=courant_diff,
        courant_reac=courant_reac,
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
        namespace=namespace,
    )
