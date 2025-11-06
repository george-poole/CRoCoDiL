from abc import ABC, abstractmethod
from typing import Iterable, Generic, TypeVar
from dataclasses import dataclass
from functools import cached_property
from typing_extensions import Self

import numpy as np
from scipy.integrate import solve_ivp


class ExprFluxModel:
    
    @staticmethod
    def ode_rhs(
        y,
        Da,
        Sr, 
        Cr,
        h0,
        epsilon,
        lmbda,
        n,
    ):
        """
        `dc⁺(t) / dt = ...`

        `dc⁻(t) / dt = ...`
        
        `dS⁺(t) / dt = ...`
        
        """
        c_plus, c_minus, s_plus = y
        f = ExprFluxModel.f(c_plus, c_minus, lmbda, n)
        h = ExprFluxModel.h(c_plus, c_minus, s_plus, Sr, Cr, h0, epsilon)
        dc_plus = -f/(1 - h) + Da * s_plus * (1 - c_plus) / (1 - s_plus)
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
        `f = λ(c⁺ - c⁻)ⁿ`
        """
        return lmbda * (c_plus - c_minus) **n
    
    @staticmethod
    def h(
        cls,
        c_plus,
        c_minus,
        s_plus,
        Sr,
        Cr,
        h0,
        epsilon,
    ):
        """
        `h = ... / ...
        """ 
        m_per_L = cls.mass_total(Cr, 0, Sr, h0, epsilon, L=1)
        numerator = m_per_L - s_plus / epsilon - (1 - s_plus) * c_plus
        denom = c_minus - s_plus/epsilon - (1 - s_plus) * c_plus
        return numerator / denom
    
    @staticmethod
    def mD():
        ...

    @staticmethod
    def mC():
        ...


@dataclass
class FluxModel:
    t: Iterable[float]
    y: Iterable[float]
    Da: float
    h0: float
    sr: float
    cr: float

    @property
    def c_plus():
        ...










P = TypeVar('P')
class Model(ABC, Generic[P]):

    def __init__(self, t, parameters: P):
        self._t = t
        self._parameters = parameters

    @property
    def t(self):
        """
        `t` timeseries array
        """
        return self._t

    @property
    def parameters(self) -> P:
        return self._parameters

    @property
    @abstractmethod
    def c_plus(self): # np.ndarray
        ...
    
    @property
    @abstractmethod
    def c_minus(self):
        ...
    
    @property
    @abstractmethod
    def s_plus(self):
        ...
    
    @property
    @abstractmethod
    def mC(self):
        """
        mᶜ(t)
        """
        ...
    
    @property
    @abstractmethod
    def mD(self):
        """
        mᴰ(t)
        """
        ...

    @property
    @abstractmethod
    def f(self):
        """
        F(t)
        """
        ...


@dataclass
class ModelParams:
    """
    `Da, Sᵣ, cᵣ, h₀, ε, L`
    """
    Da: float
    Sr: float
    Cr: float
    h0: float
    epsilon: float
    L: float


@dataclass
class FluxModelParams(ModelParams):
    """
    `Da, Sᵣ, cᵣ, h₀, ε, L, λ, n`
    """
    lmbda: float
    n: float
    ode_method: str


@dataclass
class FluxModel(Model[FluxModelParams]):

    def __init__(
        self, 
        t: Iterable[float], 
        Da: float, Sr: float, Cr, h0, epsilon, L, lmbda, n, 
        method='RK45',
    ):
        super().__init__(t, FluxModelParams(Da, Sr, Cr, h0, epsilon, L, lmbda, n, method))
    
    @property
    def c_plus(self):
        """
        `c⁺(t)`
        """
        return self.ode_solution.y[0]
    
    @property
    def c_minus(self):
        """
        `c⁻(t)`
        """
        return self.solution.y[1]
    
    @property
    def s_plus(self):
        """
        `S⁺(t)`
        """
        return self.solution.y[2]
        
    @property
    def h(self):
        """
        `h(t)`
        """
        return self.height(
            self.c_plus, self.c_minus, self.s_plus, 
            self.parameters.Sr, self.parameters.Cr, self.parameters.h0, self.parameters.epsilon)
    
    @property
    def mD(self):
        return self.mass_dissolved(
            self.c_plus, self.c_minus, self.s_plus, 
            self.parameters.Sr, self.parameters.Cr, self.parameters.h0, self.parameters.epsilon, self.parameters.L,
        )
    
    @classmethod
    def mass_dissolved(
        cls,
        c_plus,
        c_minus,
        s_plus,
        Sr,
        Cr,
        h0,
        epsilon,
        L,
    ) -> float | np.ndarray:
        """
        `mᴰ(c⁺, c⁻, S⁺)` expression
        """
        h = cls.height(c_plus, c_minus, s_plus, Sr, Cr, h0, epsilon)
        return L * (1 - h) * (1 - s_plus) * c_plus + L * h * c_minus
    
    @property
    def mC(self):
        return self.mass_capillary_trapped(
            self.c_plus, self.c_minus, self.s_plus, 
            self.parameters.Sr, self.parameters.Cr, self.parameters.h0, self.parameters.epsilon, self.parameters.L)
    
    @classmethod
    def mass_capillary_trapped(
        cls,
        c_plus,
        c_minus,
        s_plus,
        Sr,
        Cr,
        h0,
        epsilon,
        L,
    ) -> float | np.ndarray:
        """
        `mᴰ(c⁺, c⁻, S⁺)` expression
        """
        h = cls.height(c_plus, c_minus, s_plus, Sr, Cr, h0, epsilon)
        return L * (1 - h) * s_plus / epsilon
    
    @classmethod
    def mass_total(
        cls,
        c_plus,
        c_minus,
        s_plus,
        h,
        epsilon,
        L,
    ):
        rhs = (1 - h) * (1 - s_plus) * c_plus + h * c_minus + (1 - h) * s_plus / epsilon
        return L * rhs
    
    @property
    def f(self):
        return self.flux(self.c_plus, self.c_minus, self.parameters.lmbda, self.parameters.n)
    
    @classmethod
    def ode_solve(
        cls,
        t, 
        Da,
        Sr,
        Cr, 
        h0,
        epsilon,
        lmbda,
        n,
        method: str = 'RK45',
    ) -> tuple:
        """
        Solves coupled ODE system for `c⁺(t), c⁻(t), S⁺(t)`
        """

        ics = ([Cr, 0.0, Sr])
        solution = solve_ivp(
            cls.ode_system, 
            (t[0], t[-1]), 
            ics, method, t, 
            args=(Da, Sr, Cr, h0, epsilon, lmbda, n),
        )
        if solution.success:
            return solution
        else:
            raise RuntimeError(f'{solution.message}')
        
    @cached_property
    def ode_solution(self):
        """
        Solution `c⁺(t), c⁻(t), S⁺(t)` to the coupled ODE
        """
        return self.ode_solve(
            self.t, 
            self.parameters.Da, 
            self.parameters.Sr, 
            self.parameters.Cr, 
            self.parameters.h0, 
            self.parameters.epsilon, 
            self.parameters.lmbda, 
            self.parameters.n,
            self.parameters.ode_method, 
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

