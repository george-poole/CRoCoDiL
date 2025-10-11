from typing import Callable
import numpy as np


def heaviside(
    fx: Callable[[np.ndarray], np.ndarray],
    f_plus: float = 1.0,
    f_minus: float = 0.0,
    eps: float | tuple[float, float] | None = None,
):
    """
    `H(f(x))` \\
    `= f₊` if `f(x) > 0` \\
    `= f₋` otherwise

    Optionally to smooth out discontinuities

    `H(f(x))` \\
    `= f₊ tanh(f(x) / ϵ)` if `f(x) > 0` \\
    `= f₋` otherwise

    `H(f(x))` \\
    `= (f₊ - f₋) tanh(f(x) / ϵ₊) / 2 + (f₊ + f₋) / 2` if `f(x) > 0` \\
    `= (f₊ - f₋) tanh(f(x) / ϵ₋) / 2 + (f₊ + f₋) / 2` otherwise
    """
    ind = lambda x: (fx(x) >= 0)

    if isinstance(eps, float):
        return lambda x: f_plus * np.tanh(fx(x) / eps) * ind(x) + f_minus
    elif isinstance(eps, tuple):
        eps_lt, eps_gt = eps
        return lambda x: (
            (0.5 * (f_plus - f_minus) * np.tanh(fx(x) / eps_gt) + 0.5 * (f_plus + f_minus)) * ind(x)
            + 0.5 * (f_plus - f_minus) * np.tanh(fx(x) / eps_lt) + 0.5 * (f_plus + f_minus)
        )
    else:
        return lambda x: f_plus * ind(x) + f_minus