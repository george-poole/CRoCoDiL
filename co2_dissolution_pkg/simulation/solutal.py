from types import EllipsisType

import numpy as np

from lucifex.fdm import ConstantSeries, FiniteDifference, CN, AB
from lucifex.utils import CellType, Perturbation, cubic_noise
from lucifex.solver import OptionsPETSc, OptionsJIT, dS_solver
from lucifex.sim import create_simulation

from ..formulae.constitutive import permeability_cross_bedded
from ..formulae.secondary import flux
from ..formulae.utils import heaviside

from .domain import create_rectangle_domain
from .astract_thermosolutal import thermosolutal_convection_dissolution


@create_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def solutal_convective_dissolution_2d(
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
    Omega, dOmega = create_rectangle_domain(Lx, Ly, Nx, Ny, cell)
    s_ics = heaviside(lambda x: x[1] - h0, sr, eps=heaviside_eps) 
    c_ics = Perturbation(
        heaviside(lambda x: x[1] - h0, cr, eps=heaviside_eps),
        cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed, (0, 1)),
        [Lx, Ly],
        c_eps,
        )   
    density = lambda c: c
    reaction = lambda s, c: s * (1 - c)

    simulation = thermosolutal_convection_dissolution(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # physical
        Ra=Ra,
        Da=Da,
        epsilon=epsilon,
        # initial conditions
        s_ics=s_ics, 
        c_ics=c_ics,
        # constitutive relations
        density=density,
        reaction=reaction,
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
    )

    if secondary:
        c, u, d = simulation['c', 'u', 'd']
        Gamma = (lambda x: x[1] - h0, lambda x: x[1] - h0 / 2)
        fGamma = ConstantSeries(Omega, "f", shape=(len(Gamma), 2))
        simulation.solvers.append(
            dS_solver(fGamma, flux, Gamma, facet_side="+")(
                c[0], u[0], d[0], Ra,
            )
        )

    return simulation


@create_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def solutal_convective_dissolution_inclined(
    # mesh
    Lx: float = 5.0,
    Ly: float = 1.0,
    Nx: int = 100,
    Ny: int = 100,
    cell: str = CellType.QUADRILATERAL,
    # physical
    Ra: float = 1e3,
    Da: float = 1e2,
    epsilon: float = 1e-2,
    beta: float = 45.0,
    # permeability constitutive relation
    kappa: float = 0.1,
    vartheta: float = 45.0,
    # initial conditions
    cr: float = 0.0,
    sr: float = 0.2,
    h0: float | tuple[float, float, float] = 4.0,
    H_eps: float | None = None,
    # concentration perturbation
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
    Omega, dOmega = create_rectangle_domain(Lx, Ly, Nx, Ny, cell)
    s_ics = heaviside(lambda x: x[0] - h0, sr, eps=H_eps) 
    c_ics = Perturbation(
        heaviside(lambda x: x[0] - h0, cr, eps=H_eps),
        cubic_noise(['neumann', 'neumann'], [Lx, Ly], c_freq, c_seed, (0, 1)),
        [Lx, Ly],
        c_eps,
        )   

    permeability = lambda phi: permeability_cross_bedded(phi, kappa, vartheta * np.pi / 180, 2)
    density = lambda c: c
    reaction = lambda s, c: s * (1 - c)

    return thermosolutal_convection_dissolution(
        # domain
        Omega=Omega, 
        dOmega=dOmega, 
        # physical
        Ra=Ra,
        Da=Da,
        epsilon=epsilon,
        beta=beta,
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
    )


# f = ConstantSeries(Omega, "f", shape=(len(dOmega.union), 2))
# f_solver = ds_solver(flux, f, dOmega.union)(
#     c[0], u[0], d[0], Ra,
# )