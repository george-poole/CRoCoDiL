from typing import Iterable, Callable, Any, TypeVar, Literal
from dataclasses import dataclass
from inspect import signature

import numpy as np
from scipy.integrate import solve_ivp


class FluxModelEquations:

    @staticmethod
    def ode_solve(
        t: np.ndarray, 
        Da: float,
        epsilon: float,
        zeta0: float,
        sr: float,
        cr: float,
        flux: Callable[[float, float], float],
        gamma: float | tuple[float, float],
        constraint: Literal['zeta', 'c_plus'] | None,
        ics: tuple[float, float, float] | tuple[float, float] | None,
        method: str,
    ) -> tuple:
        """
        Solve coupled ODE system for `c⁺(t), c⁻(t), s⁺(t)`
        """
        if isinstance(gamma, (float, int)):
            gamma = (gamma, gamma)
        gamma_c, gamma_s = gamma
        if constraint is None:
            rhs = lambda _, y, *args: FluxModelEquations.ode_rhs(*y, *args)
            ics = ([cr, 0.0, sr]) if ics is None else ics
        elif constraint == 'zeta':
            rhs = lambda _, y, *args: FluxModelEquations.ode_rhs_zeta_constraint(*y, *args)
            ics = ([cr, sr]) if ics is None else ics
        elif constraint == 'c_plus':
            rhs = lambda _, y, *args: FluxModelEquations.ode_rhs_c_plus_constraint(*y, *args)
            ics = ([0.0]) if ics is None else ics
        else:
            raise ValueError(constraint)

        solution = solve_ivp(
            rhs, 
            (t[0], t[-1]), 
            ics, 
            method, 
            t, 
            args=(Da, epsilon, zeta0, sr, cr, gamma_c, gamma_s, flux),
        )
        if solution.success:
            return solution
        else:
            raise RuntimeError(f'{solution.message}')
    
    @staticmethod
    def ode_rhs(
        c_plus,
        c_minus,
        s_plus,
        Da,
        epsilon,
        zeta0,
        sr, 
        cr,
        gamma_c,
        gamma_s,
        flux,
    ):
        """
        `dc⁺(t) / dt = -F(c⁺, c⁻) / (1 - ζ(c⁺, c⁻, s⁺)) + Da Σ(s⁺, c⁺)` / (1 - s⁺)`

        `dc⁻(t) / dt = F(c⁺, c⁻) / (1 - ζ(c⁺, c⁻, s⁺))`
        
        `ds⁺(t) / dt = -εDa Σ(s⁺, c⁺)`
        """
        zeta = FluxModelEquations.zeta(c_plus, c_minus, s_plus, epsilon, zeta0, sr, cr)
        r = FluxModelEquations.Sigma(c_plus, s_plus)
        f = flux(c_plus, c_minus)
        dc_plus = -f / (1 - zeta) + gamma_c * Da * r / (1 - s_plus)
        dc_minus = f / zeta
        ds_plus = -epsilon * gamma_s * Da * s_plus *(1 - c_plus)
        return [dc_plus, dc_minus, ds_plus]
    
    @staticmethod
    def ode_rhs_zeta_constraint(
        c_plus,
        s_plus,
        Da,
        epsilon,
        zeta0,
        sr, 
        cr,
        gamma_c,
        gamma_s,
        flux,
    ):
        """
        `dc⁺(t) / dt = -F(c⁺, c⁻) / (1 - ζ₀) + Da Σ(s⁺, c⁺)` / (1 - s⁺)`
        
        `ds⁺(t) / dt = -εDa Σ(s⁺, c⁺)`
        """
        c_minus = FluxModelEquations.c_minus_height_constaint(c_plus, s_plus, epsilon, zeta0, sr, cr)
        r = FluxModelEquations.Sigma(c_plus, s_plus)
        f = flux(c_plus, c_minus)
        dc_plus = -f / (1 - zeta0) + gamma_c * Da * r / (1 - s_plus)
        ds_plus = -epsilon * gamma_s * Da * s_plus *(1 - c_plus)
        return [dc_plus, ds_plus]
    
    @staticmethod
    def c_minus_height_constaint(c_plus, s_plus, epsilon, zeta0, sr, cr):
        m0_per_vol = FluxModelEquations.m0(epsilon, zeta0, sr, cr, 1)
        c_minus = m0_per_vol
        c_minus += -(1 - zeta0) * (1 - s_plus) * c_plus 
        c_minus += -(1 - zeta0) * s_plus / epsilon
        c_minus *= (1 / zeta0)
        return c_minus
    
    @staticmethod
    def ode_rhs_c_plus_constraint(
        c_minus,
        Da,
        epsilon,
        zeta0,
        sr, 
        cr,
        gamma_c,
        gamma_s,
        flux,
    ):
        ...
    
    @staticmethod
    def Sigma(
        c_plus,
        s_plus,
    ):
        """
        `Σ(s⁺, c⁺) = s⁺(1 - c⁺)`
        """
        return s_plus * (1 - c_plus)

    @staticmethod
    def zeta(
        c_plus,
        c_minus,
        s_plus,
        epsilon,
        zeta0,
        sr,
        cr,
    ):
        """
        `ζ(c⁺, c⁻, s⁺) = ... / ...` horizontally-averaged interface height
        """ 
        m0_per_vol = FluxModelEquations.m0(epsilon, zeta0, sr, cr, 1)
        numerator = m0_per_vol - s_plus / epsilon - (1 - s_plus) * c_plus
        denom = c_minus - s_plus/epsilon - (1 - s_plus) * c_plus
        return numerator / denom
    
    @staticmethod
    def mD(
        c_plus,
        c_minus,
        s_plus,
        epsilon,
        zeta0,
        sr,
        cr,
        vol: float,
    ):
        """
        `mᴰ(c⁺, c⁻, s⁺) = ...` dissolved mass
        """
        zeta = FluxModelEquations.zeta(c_plus, c_minus, s_plus, epsilon, zeta0, sr, cr)
        m_per_vol = (1 - zeta) * (1 - s_plus) * c_plus + zeta * c_minus
        return m_per_vol * vol

    @staticmethod
    def mC(
        c_plus,
        c_minus,
        s_plus,
        epsilon,
        zeta0,
        sr,
        cr,
        vol: float,
    ):
        """
        `mᶜ(c⁺, c⁻, s⁺) = ...` capillary-trapped mass
        """
        zeta = FluxModelEquations.zeta(c_plus, c_minus, s_plus, epsilon, zeta0, sr, cr)
        m_per_vol = (1 - zeta) * s_plus / epsilon 
        return m_per_vol * vol
    
    @staticmethod
    def m_sum(
        c_plus,
        c_minus,
        s_plus,
        epsilon,
        zeta0,
        sr,
        cr,
        vol: float,
    ):
        """
        `m = mᴰ + mᶜ` total mass
        """
        return sum((
            FluxModelEquations.mD(c_plus, c_minus, s_plus, epsilon, zeta0, sr, cr, vol),
            FluxModelEquations.mC(c_plus, c_minus, s_plus, epsilon, zeta0, sr, cr, vol),
        ))

    @staticmethod
    def m0(
        epsilon,
        zeta0,
        sr, 
        cr,
        vol: float,
    ):
        """
        `m₀ = mᴰ(cᵣ, 0, sᵣ) + mᶜ(cᵣ, 0, sᵣ)` initial mass (conserved)
        """
        m_per_vol = (1 - zeta0) * (1- sr) * cr + (1 - zeta0) * sr / epsilon
        return m_per_vol * vol
    

C_PLUS = TypeVar('C_PLUS')
C_MINUS = TypeVar('C_MINUS')
FLUX = TypeVar('FLUX')
@dataclass
class FluxModel:
    """
    If `ics=None`, then default initial conditions 
    `(c⁺, c⁻, s⁺) = (cᵣ, 0, sᵣ)` are assumed.
    """
    t: Iterable[float]
    vol: float
    Da: float
    epsilon: float
    zeta0: float
    sr: float
    cr: float
    gamma: float | tuple[float, float]
    flux: Callable[[C_PLUS, C_MINUS], FLUX]
    reduced: bool = False
    ics: tuple[float, float, float] | None = None
    method: str = 'RK45'

    def __post_init__(self):
        self._ode_solution = self._pass_kws(
            FluxModelEquations.ode_solve,
        )

    def _get_kws(
        self,
        clbl: Callable,
    ) -> dict[str, Any]:
        return {k: getattr(self, k) for k in signature(clbl).parameters}
    
    def _pass_kws(
        self,
        clbl: Callable,
        **kws_extra,
    ):
        kws = self._get_kws(clbl)
        return clbl(**kws, **kws_extra)

    @property
    def c_plus(self) -> np.ndarray:
        """
        `c⁺(t)` upper concentration
        """
        return self._ode_solution.y[0]
    
    @property
    def c_minus(self) -> np.ndarray:
        """
        `c⁻(t)`lower concentration
        """
        if self.reduced:
            return self._pass_kws(FluxModelEquations.c_minus_height_constaint)
        else:
            return self._ode_solution.y[1]
    
    @property
    def s_plus(self) -> np.ndarray:
        """
        s⁺(t)` upper saturation
        """
        if self.reduced:
            index = 1
        else:
            index = 2
        return self._ode_solution.y[index]
    
    @property
    def zeta(self) -> float | np.ndarray:
        """
        `ζ(t)` horizontally-averaged interface height 
        """
        if self.reduced:
            return self.zeta0
        else:
            return self._pass_kws(FluxModelEquations.zeta)
    
    @property
    def mC(self) -> np.ndarray:
        """
        `mᶜ(t)` capillary-trapped mass
        """
        return self._pass_kws(FluxModelEquations.mC)
    
    @property
    def mD(self) -> np.ndarray:
        """
        `mᴰ(t)` dissolved mass 
        """
        return self._pass_kws(FluxModelEquations.mD)
    
    @property
    def f(self) -> np.ndarray:
        """
        `F(t) = F(c⁺(t), c⁻(t))` flux 
        """
        return self.flux(self.c_plus, self.c_minus)

    
# class DaInftyFluxModel(FluxModel):

#     def __init__(self, t, Sr, Cr, h0, epsilon, lmbda=0.001, n=2, method='RK45'):
#         Da = np.inf
#         super().__init__(t, Da, Sr, Cr, h0, epsilon, lmbda, n, method)

#     @property
#     def c_plus(self):
#         """
#         `c⁺(t) ≈ 1`
#         """
#         return np.full_like(np.array(self.t), 1.0)
    
#     @property
#     def s_plus(self):
#         """
#         S⁺(t)
#         """
#         m_per_L = self.mass_total(
#             self.parameters.Cr, 0, 
#             self.parameters.Sr, 
#             self.parameters.h0, 
#             self.parameters.epsilon, L=1)
#         numerator = m_per_L - self.h * self.c_minus - (1-self.h) * self.c_plus
#         denominator = (1 - self.h) * (1/self.parameters.epsilon - self.c_plus)
#         return numerator / denominator
    
#     @property
#     def h(self):
#         """
#         `ζ(t) ≈ ζ₀`
#         """
#         return np.full_like(np.array(self.t), self.parameters.h0)
    
#     @classmethod
#     def ode_system(
#         cls,
#         _, 
#         y,
#         h0,
#         lmbda,
#         n,
#     ):
#         """
#         `f(t, y)` right-hand side function for the decoupled ODE

#         `dc⁻(t) / dt = F(c⁻) / ζ₀`        
#         """
#         c_minus = y[0]
#         f = cls.flux(1, c_minus, lmbda, n)
#         dc_minus = f/h0
#         return dc_minus


# class EpsilonZeroModel(FluxModel):
#     ...


    # @staticmethod
    # def mass_capillary_trapped(
    #     theta: float | np.ndarray,
    #     Sr: float,
    #     h0: float,
    #     epsilon: float,
    # )  -> float | np.ndarray:
    #     """
    #     `mᶜ/L = (1 - h)Sᵣ/ε`
    #     """
    #     h = SharpFrontmodelling.height(theta, Sr, h0, epsilon)
    #     return (1 / epsilon) * (1 - h) * Sr
    

    # @classmethod
    # def concentration_upper(
    #     t, 
    #     Sr, 
    #     h0,
    #     epsilon
    # ):
    #     t_lim = (t[0], t[-1])
    #     theta_ics = (0.0, )
    #     h = cls.height
    #     dh = modelling.height_derivative
    #     rhs = lambda _, theta: flux(theta, *args, **kwargs) / (h(theta, Sr, h0, epsilon) + theta * dh(theta, Sr, h0, epsilon))
    #     solution = solve_ivp(rhs, t_lim, theta_ics, ode_method, t)
    #     if solution.success:
    #         theta = solution.y[0]
    #     else:
    #         raise RuntimeError(f'{solution.message}')
    
    # @staticmethod
    # def height(
    #     theta: float | np.ndarray,
    #     Sr: float,
    #     h0: float,
    #     epsilon: float,
    # ) -> float | np.ndarray:
    #     """
    #     `h = ...`
    #     """
    #     const = 1 + Sr * (1/epsilon - 1)
    #     return h0 * (1 - (theta / const))**-1
    
    # @staticmethod
    # def height_derivative(
    #     theta: float | np.ndarray,
    #     Sr: float,
    #     h0: float,
    #     epsilon: float,
    # ) -> float | np.ndarray:
    #     """
    #     `dh/dΘ = ...`
    #     """
    #     const = 1 + Sr * (1/epsilon - 1)
    #     return (h0 / const) * (1 - (theta / const))**-2
    

# class _FrontModel(ABC):
#     def __init__(self, theta, Sr, h0, epsilon):
#         self._theta = theta
#         self._h = self.height(theta, Sr, h0, epsilon)
#         self._mD = self.mass_dissolved(theta, Sr, h0, epsilon)
#         self._mC = self.mass_capillary_trapped(theta, Sr, h0, epsilon)

#     @property
#     def theta(self) -> np.ndarray:
#         return self._theta
    
#     @property
#     def h(self) -> np.ndarray:
#         return self._h
    
#     @property
#     def mD(self) -> np.ndarray:
#         return self._mD
    
#     @property
#     def mC(self) -> np.ndarray:
#         return self._mC
    
#     @property
#     def m_total(self) -> np.ndarray:
#         return self.mD + self.mC
    
#     @property
#     def fD(self) -> np.ndarray:
#         return self.mD / self.m_total
    
#     @property
#     def fC(self) -> np.ndarray:
#         return self.mC / self.m_total
    
#     @staticmethod
#     @abstractmethod
#     def mass_dissolved(
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#     ):
#         """
#         mass per unit horizontal length
#         `mᶜ(Θ; Sᵣ, ζ₀, ε) / L`
#         """
#         ...

#     @staticmethod
#     @abstractmethod
#     def mass_capillary_trapped(
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#     ):
#         """
#         mass per unit horizontal length
#         `mᴰ(Θ; Sᵣ, ζ₀, ε) / L`
#         """
#         ...

#     @staticmethod
#     @abstractmethod
#     def height(
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#     ):
#         """
#         `h(Θ; Sᵣ, ζ₀, ε)`
#         """
#         ...

#     @staticmethod
#     @abstractmethod
#     def height_derivative(
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#     ):
#         """
#         `dh/dΘ(Θ; Sᵣ, ζ₀, ε)`
#         """
#         ...


# class _SharpFrontModel(_FrontModel):

#     @staticmethod
#     def mass_dissolved(
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#         vol: float = 1.0,
#     ) -> float | np.ndarray:
#         """
#         `mᴰ/L = hΘ + (1 − h)(1 − Sᵣ)`
#         """
#         h = SharpFrontmodelling.height(theta, Sr, h0, epsilon)
#         return h * theta + (1 - h) * (1 - Sr)
    
#     @staticmethod
#     def mass_capillary_trapped(
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#     )  -> float | np.ndarray:
#         """
#         `mᶜ/L = (1 - h)Sᵣ/ε`
#         """
#         h = SharpFrontmodelling.height(theta, Sr, h0, epsilon)
#         return (1 / epsilon) * (1 - h) * Sr
    
#     @staticmethod
#     def height(
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#     ) -> float | np.ndarray:
#         """
#         `h = ...`
#         """
#         const = 1 + Sr * (1/epsilon - 1)
#         return h0 * (1 - (theta / const))**-1
    
#     @staticmethod
#     def height_derivative(
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#     ) -> float | np.ndarray:
#         """
#         `dh/dΘ = ...`
#         """
#         const = 1 + Sr * (1/epsilon - 1)
#         return (h0 / const) * (1 - (theta / const))**-2

    
# T = TypeVar('T', bound=FrontModel)
# P = ParamSpec('P')
# def mass_model_solver(
#     model: type[T],
#     flux: Callable[Concatenate[np.ndarray, P], np.ndarray],
#     Sr: float,
#     h0: float,
#     epsilon: float,
# ):
#     def _inner(t: Iterable[float], ode_method: str = 'RK45', *args: P.args, **kwargs: P.kwargs,) -> T:
#         t_lim = (t[0], t[-1])
#         theta_ics = (0.0, )
#         h = modelling.height
#         dh = modelling.height_derivative
#         rhs = lambda _, theta: flux(theta, *args, **kwargs) / (h(theta, Sr, h0, epsilon) + theta * dh(theta, Sr, h0, epsilon))
#         solution = solve_ivp(rhs, t_lim, theta_ics, ode_method, t)
#         if solution.success:
#             theta = solution.y[0]
#         else:
#             raise RuntimeError(f'{solution.message}')
        
#         return model(theta, Sr, h0, epsilon)
    
#     return _inner


# class ThetaFunc(Protocol):
#     def __call__(
#         self, 
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#     ) -> float | np.ndarray:
#         ...


# def fit_flux_parameter(
#     t_dns: Iterable[float],
#     quantity_dns: Iterable[float],
#     quantity_func: ThetaFunc,
#     Sr: float,
#     h0: float,
#     epsilon: float,
#     t_shift: float = 0.0 # TODO or extra fitting parameter?
# ) -> float:
#     """
#     Determines the parameter `λ` by fitting model to data from DNS.
#     """
    
#     def _mD_model(t, lmbda):
#         theta = concentration_sharp_front(t, lmbda, Sr, h0, epsilon)
#         return quantity_func(theta, Sr=Sr, zeta0=h0, epsilon=epsilon)
    
#     assert len(t_dns) == len(quantity_dns)
#     t_dns = np.asarray(t_dns)
#     quantity_dns = np.asarray(quantity_dns)
    
#     t_dns_shifted = t_dns - t_shift
#     indices = np.where(t_dns_shifted >= 0)
#     t_dns_shifted = t_dns_shifted[indices]
#     quantity_dns_shifted = quantity_dns[indices]

#     (lmbda, ), _ = curve_fit(_mD_model, t_dns_shifted, quantity_dns_shifted, bounds=(0, np.inf))
#     return lmbda

