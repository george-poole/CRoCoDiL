from typing import Callable
from typing_extensions import Unpack

import numpy as np
from dolfinx.mesh import Mesh
from ufl.core.expr import Expr

from lucifex.fem import Function, Constant
from lucifex.mesh import MeshBoundary, mesh_boundary, rectangle_mesh
from lucifex.utils import integral
from lucifex.utils.py_utils import StrEnum


class ScalingType(StrEnum):
    """
    Choice of length scale `ℒ`, velocity scale `𝒰`
    and time scale `𝒯` in the non-dimensionalization.
    """
    ADVECTIVE = 'advective'
    """
    `ℒ` = domain size \\
    `𝒰` = advective speed
    """
    DIFFUSIVE = 'diffusive'
    """
    `ℒ` = domain size \\
    `𝒰` = advective speed
    """
    ADVECTIVE_DIFFUSIVE = 'advective_diffusive'
    """
    `ℒ` = diffusive length \\
    `𝒰` = advective speed
    """
    REACTIVE = 'reactive'
    """
    `ℒ` = diffusive length \\
    `T` = reactive time
    """

    def mapping(
        self,
        Ra: float,
        Da: float,
    ) -> Callable[[str], int | float] | Callable[[Unpack[tuple[str, ...]]], list[int | float]]:
        if self.value == ScalingType.ADVECTIVE:
            Ad = 1
            Pe = Ra
            Ki = Da
            Bu = 1
            Xl = 1
        elif self.value == ScalingType.ADVECTIVE_DIFFUSIVE:
            Ad = 1
            Pe = 1
            Ki = Da / Ra
            Bu = 1
            Xl = Ra
        elif self.value == ScalingType.DIFFUSIVE:
            Ad = 1
            Pe = 1
            Ki = Ra * Da
            Bu = Ra
            Xl = 1
        elif self.value == ScalingType.REACTIVE:
            Ad = 1
            Pe = 1
            Ki = 1
            Bu = np.sqrt(Ra / Da)
            Xl = np.sqrt(Ra * Da)
        else:
            raise NotImplementedError(f'{self.value}')
        
        _names = ('Ad', 'Pe', 'Ki', 'Bu', 'Xl')

        _mapping = dict(
            zip(
                _names,
                (Ad, Pe, Ki, Bu, Xl),
            )
        )

        def _inner(*args):
            if len(args) == 1:
                key = args[0]
                try:
                    return _mapping[key]
                except KeyError:
                    raise KeyError(f"'{key}' is not in the names {_names}")
            else:
                return [_mapping[i] for i in args]
            
        return _inner
            

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
    

@integral
def mass_capillary_trapped(
    s: Function, 
    epsilon: Constant | float,
) -> Expr:
    """
    Mass capillary-trapped per unit per unit area (2D) or volume (3D)
    
    `ρᶜ = s / ε`

    for the mass integral

    `mᶜ = ∫ ρᶜ dx` 
    """
    return s / epsilon


@integral
def mass_dissolved(
    c: Function, 
    s: Function,
) -> Expr:
    """
    Mass dissolved per unit per unit area (2D) or volume (3D)

    `ρᴰ = ∫ (1 - s)·c dx` 

    for the mass integral

    `mᴰ = ∫ ρᴰ dx` 
    """
    return (1 - s) * c


def rectangle_domain(
    Lx: float,
    Ly: float,
    Nx: int,
    Ny: int,
    cell: str,
    name: str = 'LxLy',
    clockwise_names: tuple[str, str, str, str] = ('upper', 'right', 'lower', 'left'),
) -> tuple[Mesh, MeshBoundary]:
    
    mesh = rectangle_mesh(Lx, Ly, Nx, Ny, cell, name)
    boundary = mesh_boundary(
        mesh,
        {
            clockwise_names[0]: lambda x: x[1] - Ly,
            clockwise_names[1]: lambda x: x[0] - Lx,
            clockwise_names[2]: lambda x: x[1],
            clockwise_names[3]: lambda x: x[0],
        },
    )

    return mesh, boundary
