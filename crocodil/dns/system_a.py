from typing import Iterable

from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.utils import CellType
from lucifex.solver import OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation
from lucifex.utils import limits_corrector, frozen_dict

from .generic import dns_generic
from .utils import heaviside, rectangle_mesh_closure, CONVECTION_REACTION_SCALINGS


SYSTEM_A_REFERENCE = frozen_dict(
    aspect=2.0,
    Ra=1000.0,
    Da=100.0,
    epsilon=1e-2,
    zeta0=0.9,
    sr=0.2,
    cr=0.0,
)
"""
Dictionary containing reference parameters `aspect, Ra, Da, epsilon, h0, sr, cr` 
governing the physical (as opposed to numerical) behaviour of system A.
"""

def critical_sr(
    zeta0: float,
    cr: float,
    epsilon: float,
) -> float:
    """
    `sᵣ = ε ( 1 / (1 - ζ₀) - cᵣ) / (1 - εcᵣ)`
    """
    return epsilon * (-cr + 1 / (1 - zeta0)) / (1 - epsilon * cr)


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_system_a(
    # mesh
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
    zeta_eps: float | tuple[float, float] | None = None,
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
    # time step
    dt_min: float = 0.0,
    dt_max: float = 0.5,
    dt_h: str | float = "hmin",
    adv_courant: float | None = 0.5,
    diff_courant: float = 0.5,
    reac_courant: float = 0.1,
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
    c_stabilization: str | tuple[float, float] = None,
    c_limits: bool = False,
    s_limits: bool = False,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('gmres', 'ilu'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    s_petsc: OptionsPETSc | None = None,
    # optional postprocessing
    diagnostic: bool = True,
    fluxes: Iterable[tuple[str, float | int, float]] = (),
):
    """
    `Ω = [0, A·X] × [0, X]` \\
    `𝜑∂s/∂t = -εKi s(1 - c)` \\
    `ϕ∂c/∂t + 𝐮·∇c =  Di ∇·(ϕ∇c) + Ki s(1 - c)` \\
    `∇⋅𝐮 = 0` \\
    `𝐮 = -(∇p + Bu c 𝐞ʸ)` \\

    `s₀ = sᵣH(y - ζ₀) + N(𝐱)` \\
    `c₀ = cᵣH(y - ζ₀) + N(𝐱)`\\
    `𝐧⋅∇c = 0` on `∂Ω` \\
    `𝐧⋅𝐮 = 0` on `∂Ω`
    """
    # space
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra, Da)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    Lzeta = zeta0 * X
    Lzeta_eps = zeta_eps * X if zeta_eps is not None else None
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)
    # constants
    Di, Bu, Ki = scaling_map[Omega, 'Di', 'Bu', 'Ki']
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')
    # initial conditions
    s_ics = heaviside(lambda x: x[1] - Lzeta, sr - s_ampl, eps=Lzeta_eps) 
    if s_ampl:
        s_ics = SpatialPerturbation(
            s_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], s_freq, s_seed),
            [Lx, Ly],
            s_ampl,
            limits_corrector(0, sr),
        )
    c_ics = heaviside(lambda x: x[1] - Lzeta, cr - c_ampl, eps=Lzeta_eps),
    if c_ampl:
        c_ics = SpatialPerturbation(
            c_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed),
            [Lx, Ly],
            c_ampl,
            limits_corrector(0, 1),
            )  
    # constitutive
    density = lambda c: Bu * c
    dispersion = lambda phi: Di * phi
    reaction = lambda s: -Ki * s
    source = lambda s: Ki * s

    if diagnostic:
        fluxes = [('f', Lzeta, Lx), *fluxes]

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
        reaction=reaction,
        source=source,
        dispersion_solutal=dispersion,
        # time step
        dt_min=dt_min,
        dt_max=dt_max,
        dt_h=dt_h,
        adv_courant=adv_courant,
        diff_courant=diff_courant,
        reac_courant=reac_courant,
        # time discretization
        D_adv_solutal=D_adv,
        D_diff_solutal=D_diff,
        D_reac_solutal=D_reac,
        D_src_solutal=D_src,
        D_reac_evol=D_evol,
        # stabilization
        c_stabilization=c_stabilization,
        c_limits=c_limits,
        s_limits=s_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        c_petsc=c_petsc,
        s_petsc=s_petsc,
        # optional solvers
        diagnostic=diagnostic,
        fluxes_solutal=fluxes,
        namespace=[Ra, Da, Di, Bu, Ki, ('X', X)],
    )


