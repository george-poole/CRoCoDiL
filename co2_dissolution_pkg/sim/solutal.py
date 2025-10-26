from types import EllipsisType

import numpy as np

from lucifex.fem import Constant, Function
from lucifex.fdm import ConstantSeries, FiniteDifference, CN, AB
from lucifex.utils import CellType, SpatialPerturbation, cubic_noise, as_index, mesh_axes
from lucifex.solver import OptionsPETSc, OptionsJIT, integration
from lucifex.sim import configure_simulation
from lucifex.pde.constitutive import permeability_cross_bedded
from lucifex.pde.transport import flux

from .generic import thermosolutal_convection_generic
from .utils import heaviside, rectangle_domain


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def solutal_rectangle(
    # mesh
    Lx: float = 2.0,
    Ly: float = 1.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    Ra: float = 1e3,
    Da: float = 1e2,
    epsilon: float = 1e-2,
    # heaviside
    h0: float = 0.9,
    heaviside_eps: float | tuple[float, float] | None = None,
    # initial saturation
    sr: float = 0.2,
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
    D_adv: FiniteDifference | tuple[FiniteDifference, FiniteDifference] = (AB(2), CN),
    D_diff: FiniteDifference = CN,
    D_reac: FiniteDifference | tuple[FiniteDifference, FiniteDifference] = (AB(2), CN),
    D_reac_evol: FiniteDifference | tuple[FiniteDifference, FiniteDifference] = AB(1),
    # stabilization
    c_stabilization: str | tuple[float, float] = None,
    c_limits: EllipsisType | None = None,
    s_limits: EllipsisType | None = None,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc | None, OptionsPETSc | EllipsisType | None] | OptionsPETSc | None = (None, ...),
    c_petsc: OptionsPETSc | None = None,
    s_petsc: OptionsPETSc | EllipsisType | None = ...,
    # secondary
    secondary: bool = True,
):
    """
    `s(x,y,t=0) = sr · H(y - h₀)` \\
    `c(x,y,t=0) = cr · H(y - h₀) + N(x, y)`
    """
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
    density = lambda c: c
    dispersion = lambda phi, _: (1/Ra) * phi
    reaction = lambda s, c: Da * s * (1 - c)

    simulation = thermosolutal_convection_generic(
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
            ('f', ['fInterface', 'fPlus', 'fMinus', 'fHalf']), 
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
) -> np.ndarray:
    """
    Evaluates the fluxes
     
    `Fᵁ = ∫ (𝐧·𝐚)u ds`, `Fᴰ = ∫ 𝐧·(D·∇u) ds`

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
    contours = (
        lambda x: x[1] - h0_approx, 
        lambda x: x[1] - h0_plus, 
        lambda x: x[1] - h0_minus, 
        lambda x: x[1] - h0_half,
    )
    f = ConstantSeries(
        mesh, 
        ('f', ['fInterface', 'fPlus', 'fMinus', 'fHalf']), 
        shape=(len(contours), 2),
    )
    flux('dS', *contours, facet_side="+")(c[0], u[0], d[0])
    return integration(f, flux, 'dS', *contours, facet_side="+")(c[0], u[0], d[0])


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def solutal_rectangle_inclined(
    # mesh
    Lx: float = 5.0,
    Ly: float = 1.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    beta: float = 45.0,
    # physical
    Ra: float = 1e3,
    Da: float = 1e2,
    epsilon: float = 1e-2,
    # permeability constitutive relation
    kappa: float = 0.1,
    vartheta: float = 45.0,
    # initial conditions
    cr: float = 0.0,
    sr: float = 0.2,
    h0: float | tuple[float, float, float] = 4.0,
    H_eps: float | None = None,
    # concentration SpatialPerturbation
    c_eps: float = 1e-6,
    c_freq: tuple[int, int] = (8, 8),
    c_seed: tuple[int, int] = (1234, 5678),
    # time step
    dt_max: float = 0.5,
    cfl_h: str | float = "hmin",
    cfl_courant: float = 0.75,
    k_courant: float = 0.1,
    # time discretization
    D_adv: FiniteDifference | tuple[FiniteDifference, FiniteDifference] = (AB(2), CN),
    D_diff: FiniteDifference = CN,
    D_reac: FiniteDifference | tuple[FiniteDifference, FiniteDifference] = (AB(2), CN),
    # stabilization
    c_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    c_limits: bool | tuple[float, float] = False,
    s_limits = None,
    # linear algebra
    u_petsc: OptionsPETSc | None | tuple[OptionsPETSc | None, OptionsPETSc | EllipsisType | None] = (None, ...),
    c_petsc: OptionsPETSc | None = None,
    s_petsc: OptionsPETSc | EllipsisType | None = ...,
    # secondary
    secondary: bool = False,
):
    """
    `s(x,y,t=0) = sr · H(x - h₀)` \\
    `c(x,y,t=0) = cr · H(x - h₀) + N(x, y)`
    """
    Omega, dOmega = rectangle_domain(Lx, Ly, Nx, Ny, cell)
    beta_rad = beta * np.pi / 180
    beta = Constant(Omega, beta_rad, name='beta')
    egx = Constant(Omega, -np.sin(beta_rad))
    egy = Constant(Omega, -np.cos(beta_rad))

    s_ics = heaviside(lambda x: x[0] - h0, sr, eps=H_eps) 
    c_ics = SpatialPerturbation(
        heaviside(lambda x: x[0] - h0, cr, eps=H_eps),
        cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed),
        [Lx, Ly],
        c_eps,
        )   

    vartheta_rad = vartheta * np.pi / 180
    vartheta = Constant(Omega, vartheta_rad, name='vartheta')
    permeability = lambda phi: permeability_cross_bedded(phi**2, kappa, vartheta)
    density = lambda c: c
    reaction = lambda s, c: s * (1 - c)

    return thermosolutal_convection_generic(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        egx=egx,
        egy=egy,
        # physical
        Ra=Ra,
        Da=Da,
        epsilon=epsilon,
        # initial conditions
        s_ics=s_ics, 
        c_ics=c_ics, 
        # constitutive relations
        permeability=permeability,
        density=density,
        reaction=reaction,
        # time step
        dt_max=dt_max,
        cfl_h=cfl_h,
        cfl_courant=cfl_courant,
        k_courant=k_courant,
        # time discretization
        D_adv_solutal=D_adv,
        D_diff_solutal=D_diff,
        D_reac_solutal=D_reac,
        # stabilization
        c_stabilization=c_stabilization,
        c_limits=c_limits,
        s_limits=s_limits,
        # linear algebra
        flow_petsc=u_petsc,
        c_petsc=c_petsc,
        s_petsc=s_petsc,
        # optional solvers
        secondary=secondary,
        namespace_extras=[beta, vartheta],
    )


