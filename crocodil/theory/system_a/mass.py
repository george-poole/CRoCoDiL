def critical_saturation(
    zeta0: float,
    cr: float,
    epsilon: float,
) -> float:
    """
    `sᵣ = ε ( 1 / (1 - ζ₀) - cᵣ) / (1 - εcᵣ)`
    """
    return epsilon * (-cr + 1 / (1 - zeta0)) / (1 - epsilon * cr)


def mass_dissolved_asymptote(
    mass_initial: float,
    epsilon: float,
    Lx: float,
    Ly: float,
) -> float:
    """
    `mᴰ -> ... ` as `t -> ∞`
    
    with `vol(Ω) = LxLy`.
    """
    return min(mass_initial, (mass_initial - Lx * Ly / epsilon) / (1 - 1 / epsilon))


def mass_capillary_asymptote(
    mass_initial: float,
    epsilon: float,
    Lx: float,
    Ly: float,
) -> float:
    """
    `mᶜ -> ...` as `t -> ∞`

    with `vol(Ω) = LxLy`.
    """
    return max(0, (mass_initial - Lx * Ly ) / (1 - epsilon))


def mass_dissolved_initial(
    zeta0: float,
    sr: float,
    cr: float,
    Lx: float,
    Ly: float,
)-> float:
    return Lx * Ly * (1 - zeta0) * (1 - sr) * cr


def mass_capillary_initial(
    epsilon: float,
    zeta0: float,
    sr: float,
    Lx: float,
    Ly: float,
) -> float:
    return Lx * Ly * (1 - zeta0) * sr / epsilon

