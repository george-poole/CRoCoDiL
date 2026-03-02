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


# TODO effects of scaling? other versions?
def threshold_wavelength(
    Ra: float,
    Ly: float,
) -> float:
    """
    `λ = 90 Ly / Ra`
    """
    return 90.0 * Ly / Ra


def threshold_dx(
    Ra: float,
    Ly: float,
    n_per_cell: int,
):
    """
    `Δx ≤ λ / N` to resolve instabilities
    """
    return threshold_wavelength(Ra, Ly) / n_per_cell


def threshold_Nx(
    Ra: float,
    Lx: float,
    Ly: float,
    n_per_cell: int,
):
    """
    `Nₓ ≥ n Lₓ / λ` to resolve instabilities
    """
    return np.ceil(Lx / threshold_dx(Ra, Ly, n_per_cell))


def threshold_rayleigh(
    Lx: float,
    Ly: float,
    Nx: int,
    n_per_cell: int,
):
    """
    `Ra ≤ 90 Ly Nₓ / n Lₓ` to resolve instabilities
    """
    return 90.0 * Ly * Nx / (n_per_cell * Lx)

