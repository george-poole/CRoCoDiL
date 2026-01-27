import numpy as np

from lucifex.fem import Constant, Function, SpatialPerturbation, cubic_noise
from lucifex.fdm import ConstantSeries, FiniteDifference, FiniteDifferenceArgwise, CN, AB, AM
from lucifex.utils import CellType, as_index, mesh_axes
from lucifex.solver import OptionsPETSc, OptionsJIT, integration
from lucifex.sim import configure_simulation
from lucifex.pde.advection_diffusion import flux
from lucifex.utils.dofs_utils import limits_corrector

from crocodil.dns import dns_generic, heaviside, rectangle_mesh_closure, CONVECTION_REACTION_SCALINGS


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
    D_diff: FiniteDifference = CN,
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
    # secondary
    secondary: bool = True,
):
    """
    `ϕ∂c/∂t + 𝐮·∇c =  1/Pe ∇·(D(ϕ,𝐮)·∇c) + KiR(s,c)` \\ 
    `∇⋅𝐮 = 0` \\
    `𝐮 = -(∇p + Bu ρ(c)e₉)` \\
    `𝜑∂s/∂t = -εKiR(s,c)`

    `scaling` determines `Pe, Ki, Bu, Xl` from `Ra, Da`. \\
    `Ω = [0, aspect·Xl] × [0, Xl]`

    `s₀ = sᵣ · H(y - h₀)` plus optional noise \\
    `c₀ = cᵣ · H(y - h₀)` plus optional noise
    """
    scaling_map = CONVECTION_REACTION_SCALINGS[scaling](Ra, Da)

    Xl = scaling_map['Xl']
    Lx = aspect * Xl
    Ly = 1.0 * Xl
    Omega, dOmega = rectangle_mesh_closure(Lx, Ly, Nx, Ny, cell)
    Pe, Bu, Ki = scaling_map[Omega, 'Pe', 'Bu', 'Ki']
    Ra = Constant(Omega, Ra, 'Ra')
    Da = Constant(Omega, Da, 'Da')

    s_ics = heaviside(lambda x: x[1] - h0, sr, eps=h0_eps) 
    if s_ampl:
        s_ics = SpatialPerturbation(
            s_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], s_freq, s_seed),
            [Lx, Ly],
            s_ampl,
            limits_corrector(0, sr),
        )

    c_ics = heaviside(lambda x: x[1] - h0, cr, eps=h0_eps),
    if c_ampl:
        c_ics = SpatialPerturbation(
            c_ics,
            cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed),
            [Lx, Ly],
            c_ampl,
            limits_corrector(0, 1),
            )  
         
    density = lambda c: Bu * c
    dispersion = lambda phi: (1/Pe) * phi
    reaction = lambda s: -Ki * s
    source = lambda s: Ki * s

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
        diagnostic=secondary,
        namespace=(Ra, Da),
    )

    if secondary:
        c, u, d = simulation['c', 'u', 'd']
        f = ConstantSeries(
            Omega, 
            ('f', ['fInterface', 'fPlus', 'fMinus', 'fMid']), 
            shape=(4, 2),
        )
        flux_solver = integration(f, interfacial_flux)(c[0], u[0], d[0], h0, Lx)
        simulation.solvers.append(flux_solver)

    return simulation


def interfacial_flux(
    c: Function,
    u: Function,
    d: Function,
    h0: float,
    Lx: float,
    tol: float | None = 1e-6,
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

    ineq = lambda aprx, trgt: aprx <= trgt and np.abs(aprx - trgt) < tol
    ineq_msg = 'Mesh resolution must be chosen such that `h0` is aligned with cell facets.'
    h0_index = as_index(y_axis, h0, ineq, ineq_msg)
    h0_approx = y_axis[h0_index]
    h0_plus = y_axis[h0_index + 1]
    h0_minus = y_axis[h0_index - 1]
    h0_mid = y_axis[int(0.5 * h0_index)]
    return (1 / Lx) * flux(
        'dS', 
        lambda x: x[1] - h0_approx, 
        lambda x: x[1] - h0_plus, 
        lambda x: x[1] - h0_minus, 
        lambda x: x[1] - h0_mid,
        facet_side="+",
    )(c, u, d)
