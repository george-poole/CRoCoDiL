import numpy as np

from lucifex.utils.py_utils import FloatEnum
from lucifex.pde.scaling import ScalingChoice


class CONVECTION_CONSTANTS(FloatEnum):
    """
    Constants from the theory of convection in porous media
    """

    RA_CRITICAL = 4 * np.pi **2
    """
    Critical Rayleigh for the onset of Rayleigh-Benard convection
    """
    FLUX_FACTOR = 0.008
    """
    TODO
    """


CONVECTION_REACTION_SCALINGS = ScalingChoice(
    ('Ad', 'Di', 'Ki', 'Bu', 'X'),
    lambda Ra, Da=0: {
        'advective': (1, 1/Ra, Da, 1, 1),
        'diffusive': (1, 1, Ra * Da, Ra, 1),
        'advective_diffusive': (1, 1, Da/Ra, 1, Ra),
        'reactive': (1, 1, 1, np.sqrt(Ra / Da) if Da else np.inf, np.sqrt(Ra * Da)),
    }
)
"""
Choice of length scale `ℒ`, velocity scale `𝒰`
and time scale `𝒯` in the non-dimensionalization.

`'advective'` \\
`ℒ` = domain size \\
`𝒰` = advective speed

`'diffusive'` \\
`ℒ` = domain size \\
`𝒰` = diffusive speed

`'advective_diffusive'` \\
`ℒ` = diffusive length \\
`𝒰` = advective speed

`'reactive'` \\
`ℒ` = diffusive length \\
`𝒯` = reactive time
"""


def threshold_wavelength(
    Ra: float,
    Ly: float,
    factor: float = 90,
) -> float:
    """
    `λ = 90 Ly / Ra`
    """
    return factor * Ly / Ra


def threshold_dx(
    Ra: float,
    Ly: float,
    n_per_cell: int,
    factor: float = 90,
):
    """
    `Δx ≤ λ / n` to resolve instabilities
    """
    return threshold_wavelength(Ra, Ly, factor) / n_per_cell


def threshold_Nx(
    Ra: float,
    Lx: float,
    Ly: float,
    n_per_cell: int,
    factor: float = 90,
):
    """
    `Nx ≥ n Lx / λ` to resolve instabilities
    """
    return np.ceil(Lx / threshold_dx(Ra, Ly, n_per_cell, factor))


def threshold_rayleigh(
    Lx: float,
    Ly: float,
    Nx: int,
    n_per_cell: int,
    factor: float = 90,
):
    """
    `Ra ≤ 90 Ly Nx / n Lx` to resolve instabilities
    """
    return factor * Ly * Nx / (n_per_cell * Lx)


def diffusive_Nx(
    Ra: float,
    Lx: float,
    courant: float = 1.0, 
) -> float:
    """
    `Nx ≥ Lx (Ra / 2 cD)` for stability where `cD` is the diffusive Courant number.
    """
    return np.ceil(Lx * np.sqrt(0.5 * Ra / courant))

