from typing import Iterable

from lucifex.fem import Constant, SpatialPerturbation, cubic_noise
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.utils import CellType
from lucifex.solver import OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation
from lucifex.utils.dofs_utils import limits_corrector

from .generic import dns_generic
from .utils import heaviside, rectangle_mesh_closure, CONVECTION_REACTION_SCALINGS


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
    h0: float = 0.9,
    h0_eps: float | tuple[float, float] | None = None,
    # initial saturation
    sr: float = 0.2,
    s_ampl: float | None = None,
    s_freq: tuple[int, int] | None = None,
    s_seed: tuple[int, int] | None = None,
    # initial concentration
    cr: float = 1.0,
    c_ampl: float | None = 1e-6,
    c_freq: tuple[int, int] | None = (16, 16),
    c_seed: tuple[int, int] | None = (1234, 5678),
    # time step
    dt_min: float = 0.0,
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.5,
    r_courant: float = 0.1,
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
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
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

    `s₀ = sᵣH(y - h₀) + N(𝐱)` \\
    `c₀ = cᵣH(y - h₀) + N(𝐱)`\\
    `𝐧⋅∇c = 0` on `∂Ω` \\
    `𝐧⋅𝐮 = 0` on `∂Ω`
    """
    # space
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra, Da)
    X = scaling_map['X']
    Lx = aspect * X
    Ly = 1.0 * X
    h0_X = h0 * X
    h0_eps_X = h0_eps * X if h0_eps is not None else None
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)
    # constants
    Di, Bu, Ki = scaling_map[Omega, 'Di', 'Bu', 'Ki']
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')
    # initial conditions
    s_ics = heaviside(lambda x: x[1] - h0_X, sr, eps=h0_eps_X) 
    if s_ampl:
        s_ics = SpatialPerturbation(
            s_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], s_freq, s_seed),
            [Lx, Ly],
            s_ampl,
            limits_corrector(0, sr),
        )

    c_ics = heaviside(lambda x: x[1] - h0_X, cr, eps=h0_eps_X),
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
        fluxes = [('f', h0, Lx), *fluxes]

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
        cfl_h=cfl_h,
        cfl_courant=cfl_courant,
        r_courant=r_courant,
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


