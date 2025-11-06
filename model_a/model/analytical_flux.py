from abc import ABC, abstractmethod
from typing import Iterable, Generic, TypeVar
from dataclasses import dataclass
from functools import cached_property
from typing_extensions import Self

import numpy as np
from scipy.integrate import solve_ivp


class FluxModelEquations:
    
    @staticmethod
    def ode_rhs(
        y,
        Da,
        epsilon,
        h0,
        sr, 
        cr,
        lmbda,
        n,
    ):
        """
        `dc⁺(t) / dt = -F(c⁺, c⁻) / (1 - h(t)) + Da R(s⁺, c⁺)` / (1 - s⁺)`

        `dc⁻(t) / dt = F(c⁺, c⁻) / (1 - h(t))`
        
        `ds⁺(t) / dt = -εDa R(s⁺, c⁺)`
        """
        c_plus, c_minus, s_plus = y
        f = FluxModelEquations.f(c_plus, c_minus, lmbda, n)
        h = FluxModelEquations.h(c_plus, c_minus, s_plus, sr, cr, h0, epsilon)
        r = FluxModelEquations.r(s_plus, c_plus)
        dc_plus = -f/(1 - h) + Da * r / (1 - s_plus)
        dc_minus = f/h
        ds_plus = -epsilon * Da * s_plus *(1 - c_plus)
        return [dc_plus, dc_minus, ds_plus]
    
    @staticmethod
    def f(
        c_plus,
        c_minus,
        lmbda, 
        n,
    ):
        """
        `F(c⁺, c⁻) = λ(c⁺ - c⁻)ⁿ`
        """
        return lmbda * (c_plus - c_minus) **n
    
    @staticmethod
    def r(
        s_plus,
        c_plus,
    ):
        """
        `R(s⁺, c⁺) = s⁺(1 - c⁺)`
        """
        return s_plus * (1 - c_plus)

    @staticmethod
    def h(
        c_plus,
        c_minus,
        s_plus,
        epsilon,
        sr,
        cr,
        h0,
    ):
        """
        horizontally-averaged interface height \\
        `h = ... / ...`
        """ 
        m_per_Lx = FluxModelEquations.mT(epsilon, h0, sr, cr)
        numerator = m_per_Lx - s_plus / epsilon - (1 - s_plus) * c_plus
        denom = c_minus - s_plus/epsilon - (1 - s_plus) * c_plus
        return numerator / denom
    
    @staticmethod
    def mD(
        c_plus,
        c_minus,
        s_plus,
        epsilon,
        sr,
        cr,
        h0,
        Lx=1,
    ):
        """
        dissolved mass \\
        `mᴰ(c⁺, c⁻, s⁺) = ...`
        """
        h = FluxModelEquations.h(c_plus, c_minus, s_plus, epsilon, sr, cr, h0)
        m_per_Lx = (1 - h) * (1 - s_plus) * c_plus + h * c_minus
        return m_per_Lx * Lx

    @staticmethod
    def mC(
        c_plus,
        c_minus,
        s_plus,
        epsilon,
        sr,
        cr,
        h0,
        Lx=1,
    ):
        """
        capillary-trapped mass \\
        `mᶜ(c⁺, c⁻, s⁺) = ...`
        """
        h = FluxModelEquations.h(c_plus, c_minus, s_plus, epsilon, sr, cr, h0)
        m_per_Lx = (1 - h) * s_plus / epsilon 
        return m_per_Lx * Lx

    @staticmethod
    def mT(
        epsilon,
        h0,
        sr, 
        cr,
        Lx=1,
    ):
        """
        total mass (conserved) \\
        `mᵀ = mᴰ(cᵣ, 0, sᵣ) + mᶜ(cᵣ, 0, sᵣ)`
        """
        m_per_Lx  =(1 - h0) * (1- sr) * cr + (1 - h0) * sr / epsilon
        return m_per_Lx * Lx
    

    @staticmethod
    def ode_solve(
        t, 
        Da,
        epsilon,
        h0,
        sr,
        cr, 
        lmbda,
        n,
        method: str = 'RK45',
    ) -> tuple:
        """
        Solves coupled ODE system for `c⁺(t), c⁻(t), S⁺(t)`
        """

        ics = ([cr, 0.0, sr])
        solution = solve_ivp(
            FluxModelEquations.ode_rhs, 
            (t[0], t[-1]), 
            ics, 
            method, 
            t, 
            args=(Da, epsilon, h0, sr, cr, lmbda, n),
        )
        if solution.success:
            return solution
        else:
            raise RuntimeError(f'{solution.message}')
    



@dataclass(frozen=True)
class FluxModel:
    t: Iterable[float]
    y: Iterable[float]
    Lx: float
    Da: float
    epsilon: float
    h0: float
    sr: float
    cr: float
    lmbda: float
    n: float

    def __post_init__(self):
        self._ode_solution = FluxModelEquations.ode_solve(
            self.t,
            self.Da,
            self.epsilon,
            self.h0,
            self.sr,
            self.cr,
            self.lmbda,
            self.n

        )

    @property
    def c_plus(self) -> np.ndarray:
        """
        `c⁺(t)`
        """
        return self._ode_solution.y[0]
    
    @property
    def c_minus(self) -> np.ndarray:
        """
        `c⁻(t)`
        """
        return self._ode_solution.y[1]
    
    @property
    def s_plus(self) -> np.ndarray:
        """
        `s⁺(t)`
        """
        return self._ode_solution.y[2]
    
    @property
    def h(self) -> np.ndarray:
        """
        `h(t)`
        """
        return FluxModelEquations.h(
            self.c_plus, 
            self.c_minus, 
            self.s_plus,
            self.epsilon,
            self.sr,
            self.cr,
            self.h0,
        )
    
    @property
    def mC(self) -> np.ndarray:
        """
        `mᶜ(t)`
        """
        return FluxModelEquations.mC(
            self.c_plus,
            self.c_minus,
            self.s_plus,
            self.epsilon,
            self.sr,
            self.cr,
            self.h0,
            self.Lx,
        )
    
    @property
    def mD(self) -> np.ndarray:
        """
        `mᴰ(t)`
        """
        return FluxModelEquations.mD(
            self.c_plus,
            self.c_minus,
            self.s_plus,
            self.epsilon,
            self.sr,
            self.cr,
            self.h0,
            self.Lx,
        )

    
class DaInftyFluxModel(FluxModel):

    def __init__(self, t, Sr, Cr, h0, epsilon, lmbda=0.001, n=2, method='RK45'):
        Da = np.inf
        super().__init__(t, Da, Sr, Cr, h0, epsilon, lmbda, n, method)

    @property
    def c_plus(self):
        """
        `c⁺(t) ≈ 1`
        """
        return np.full_like(np.array(self.t), 1.0)
    
    @property
    def s_plus(self):
        """
        S⁺(t)
        """
        m_per_L = self.mass_total(
            self.parameters.Cr, 0, 
            self.parameters.Sr, 
            self.parameters.h0, 
            self.parameters.epsilon, L=1)
        numerator = m_per_L - self.h * self.c_minus - (1-self.h) * self.c_plus
        denominator = (1 - self.h) * (1/self.parameters.epsilon - self.c_plus)
        return numerator / denominator
    
    @property
    def h(self):
        """
        `h(t) ≈ h₀`
        """
        return np.full_like(np.array(self.t), self.parameters.h0)
    
    @classmethod
    def ode_system(
        cls,
        _, 
        y,
        h0,
        lmbda,
        n,
    ):
        """
        `f(t, y)` right-hand side function for the decoupled ODE

        `dc⁻(t) / dt = F(c⁻) / h₀`        
        """
        c_minus = y[0]
        f = cls.flux(1, c_minus, lmbda, n)
        dc_minus = f/h0
        return dc_minus


class EpsilonZeroModel(FluxModel):
    ...

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
    #     h = SharpFrontModel.height(theta, Sr, h0, epsilon)
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
    #     dh = model.height_derivative
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
#         `mᶜ(Θ; Sᵣ, h₀, ε) / L`
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
#         `mᴰ(Θ; Sᵣ, h₀, ε) / L`
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
#         `h(Θ; Sᵣ, h₀, ε)`
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
#         `dh/dΘ(Θ; Sᵣ, h₀, ε)`
#         """
#         ...


# class _SharpFrontModel(_FrontModel):

#     @staticmethod
#     def mass_dissolved(
#         theta: float | np.ndarray,
#         Sr: float,
#         h0: float,
#         epsilon: float,
#         Lx: float = 1.0,
#     ) -> float | np.ndarray:
#         """
#         `mᴰ/L = hΘ + (1 − h)(1 − Sᵣ)`
#         """
#         h = SharpFrontModel.height(theta, Sr, h0, epsilon)
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
#         h = SharpFrontModel.height(theta, Sr, h0, epsilon)
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
#         h = model.height
#         dh = model.height_derivative
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
#         return quantity_func(theta, Sr=Sr, h0=h0, epsilon=epsilon)
    
#     assert len(t_dns) == len(quantity_dns)
#     t_dns = np.asarray(t_dns)
#     quantity_dns = np.asarray(quantity_dns)
    
#     t_dns_shifted = t_dns - t_shift
#     indices = np.where(t_dns_shifted >= 0)
#     t_dns_shifted = t_dns_shifted[indices]
#     quantity_dns_shifted = quantity_dns[indices]

#     (lmbda, ), _ = curve_fit(_mD_model, t_dns_shifted, quantity_dns_shifted, bounds=(0, np.inf))
#     return lmbda

