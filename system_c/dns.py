from types import EllipsisType

import numpy as np

from lucifex.fem import Constant
from lucifex.fdm import FiniteDifference, FiniteDifferenceArgwise, CN, AB2
from lucifex.utils import CellType, SpatialPerturbation, cubic_noise
from lucifex.solver import OptionsPETSc, OptionsJIT
from lucifex.sim import configure_simulation
from lucifex.pde.constitutive import permeability_cross_bedded

from crocodil.dns.generic import dns_generic
from crocodil.dns.utils import heaviside, rectangle_domain


@configure_simulation(
    jit=OptionsJIT("./__jit__/"),
)
def dns_model_c(
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
    r_courant: float = 0.1,
    # time discretization
    D_adv: FiniteDifference
    | FiniteDifferenceArgwise = (AB2 @ CN),
    D_diff: FiniteDifference = CN,
    D_reac: FiniteDifference
    | FiniteDifferenceArgwise = (AB2 @ CN),
    # stabilization
    c_stabilization: str | tuple[str, float] | tuple[float, float] = None,
    c_limits: bool | tuple[float, float] = False,
    s_limits = None,
    # linear algebra
    flow_petsc: tuple[OptionsPETSc, OptionsPETSc | None] 
    | OptionsPETSc = (OptionsPETSc('cg', 'gamg'), None),
    c_petsc: OptionsPETSc = OptionsPETSc('gmres', 'ilu'),
    s_petsc: OptionsPETSc | None = None,
    # secondary
    secondary: bool = False,
):
    """
    `s(x,y,t=0) = sᵣ · H(x - h₀)` \\
    `c(x,y,t=0) = cᵣ · H(x - h₀) + N(x, y)`
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

    return dns_generic(
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
        r_courant=r_courant,
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
        namespace_extras=[beta, vartheta],
    )


