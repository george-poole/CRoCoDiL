import numpy as np


def mass_dissolved_initial(
    zeta0: float,
    sr: float,
    cr: float,
    Router: float,
)-> float:
    return upper_area(zeta0, Router) * (1 - sr) * cr


def mass_capillary_initial(
    epsilon: float,
    zeta0: float,
    sr: float,
    Router: float,
) -> float:
    return upper_area(zeta0, Router) * sr / epsilon


def upper_area(
    zeta0: float,
    Router: float,
):
    return (
        (Router ** 2) * (
            np.arccos(zeta0) 
            - zeta0 * np.sqrt(1 - zeta0 ** 2)
        )
    )