from types import EllipsisType

import numpy as np

from lucifex.fem import Constant, Function
from lucifex.fdm import ConstantSeries, FiniteDifference, FiniteDifferenceArgwise, CN, AB
from lucifex.utils import CellType, SpatialPerturbation, cubic_noise, as_index, mesh_axes
from lucifex.solver import OptionsPETSc, OptionsJIT, integration
from lucifex.sim import configure_simulation
from lucifex.pde.transport import flux

from co2_pkg.sim import dns_generic, heaviside, rectangle_domain, ScalingType


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
    scaling: ScalingType = ScalingType.ADVECTIVE,
    Ra: float = 1e3,
    Da: float = 1e2,
    epsilon: float = 1e-2,
    # initial saturation
    sr: float = 0.2,
    h0: float = 0.9,
    heaviside_eps: float | tuple[float, float] | None = None,
    # initial concentration
    cr: float = 0.0,
    c_eps: float = 1e-6,
    c_freq: tuple[int, int] = (16, 16),
    c_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_min: float = 0.0,
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    k_courant: float = 0.1,
    # time discretization
    D_adv: FiniteDifference
    | FiniteDifferenceArgwise = (AB(2) @ CN),
    D_diff: FiniteDifference = CN,
    D_reac: FiniteDifference 
    | FiniteDifferenceArgwise = (AB(2) @ CN),
    D_reac_evol: FiniteDifference 
    | FiniteDifferenceArgwise = AB(1),
    # stabilization
    c_stabilization: str | tuple[float, float] = None,
    c_limits: EllipsisType | None = None,
    s_limits: EllipsisType | None = None,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    s_petsc: OptionsPETSc | None = None,
    # secondary
    secondary: bool = True,
):
    """
    `scaling` determines `Pe, Ki, Bu, Xl` from `Ra, Da`.

    `Ω = [0, aspect·Xl] × [0, Xl]`

    `ϕ∂c/∂t + 𝐮·∇c =  1/Pe ∇·(D(ϕ,𝐮)·∇c) + KiR(s,c,θ)` \\ 
    `ϕ∂θ/∂t + 𝐮·∇θ = 1/LePe ∇·(G(ϕ,𝐮)·∇θ)`\\
    `∇⋅𝐮 = 0` \\
    `𝐮 = -(∇p + Bu ρ(c,θ)e₉)` \\
    `𝜑∂s/∂t = -εKiR(s,c,θ)`

    `s₀ = sᵣ · H(y - h₀)` \\
    `c₀ = cᵣ · H(y - h₀) + N(x, y)`
    """
    Pe, Ki, Bu, Xl = ScalingType(scaling).mapping(Ra, Da, ['Pe', 'Ki', 'Bu', 'Xl'])

    Lx = aspect * Xl
    Ly = 1.0 * Xl
    Omega, dOmega = rectangle_domain(Lx, Ly, Nx, Ny, cell)
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')
    s_ics = heaviside(lambda x: x[1] - h0, sr, eps=heaviside_eps) 
    c_ics = SpatialPerturbation(
        heaviside(lambda x: x[1] - h0, cr, eps=heaviside_eps),
        cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed),
        [Lx, Ly],
        c_eps,
        )   
    density = lambda c: Bu * c
    dispersion = lambda phi, _: (1/Pe) * phi
    reaction = lambda s, c: Ki * s * (1 - c)

    simulation = dns_generic(
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
        dispersion_solutal=dispersion,
        dispersion_thermal=None,
        # time step
        dt_min=dt_min,
        dt_max=dt_max,
        cfl_h=cfl_h,
        cfl_courant=cfl_courant,
        k_courant=k_courant,
        # time discretization
        D_adv_solutal=D_adv,
        D_diff_solutal=D_diff,
        D_reac_solutal=D_reac,
        D_reac_evol=D_reac_evol,
        # stabilization
        c_stabilization=c_stabilization,
        c_limits=c_limits,
        s_limits=s_limits,
        # linear algebra
        flow_petsc=flow_petsc,
        c_petsc=c_petsc,
        s_petsc=s_petsc,
        # optional solvers
        secondary=secondary,
        namespace_extras=[Ra, Da, ],
    )

    if secondary:
        c, u, d = simulation['c', 'u', 'd']
        f = ConstantSeries(
            Omega, 
            ('f', ['fInterface', 'fPlus', 'fMinus', 'fMid']), 
            shape=(4, 2),
        )
        flux_solver = integration(f, interfacial_flux)(c[0], u[0], d[0], h0)
        simulation.solvers.append(flux_solver)

    return simulation


def interfacial_flux(
    c: Function,
    u: Function,
    d: Function,
    h0: float,
    Lx: float = 1.0,
) -> np.ndarray:
    """
    Evaluates the advective and diffusive fluxes per unit length
     
    `Fᵁ = 1 / Lx ∫ (𝐧·𝐚)u ds` \\
    `Fᴰ = 1 / Lx ∫ 𝐧·(D·∇u) ds`

    at heights
    
    `y ≃ h₀, h₀⁺, h₀⁻, h₀/2`
    """
    mesh = c.function_space.mesh
    y_axis = mesh_axes(mesh)[1]
    h0_index = as_index(y_axis, h0, less_than=True)
    h0_approx = y_axis[h0_index]
    h0_plus = y_axis[h0_index + 2]
    h0_minus = y_axis[h0_index - 2]
    h0_half = y_axis[int(0.5 * h0_index)]
    return (1 / Lx) * flux(
        'dS', 
        lambda x: x[1] - h0_approx, 
        lambda x: x[1] - h0_plus, 
        lambda x: x[1] - h0_minus, 
        lambda x: x[1] - h0_half,
        facet_side="+",
    )(c, u, d)

