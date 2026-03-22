from typing import Iterable, Callable, TypeVar, Literal
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from lucifex.utils.py_utils import arity

from .model import Model


class LateTimeFormulae:

    @staticmethod
    def ode_solve(
        t: np.ndarray, 
        Da: float,
        epsilon: float,
        zeta0: float,
        sr: float,
        cr: float,
        flux: Callable[[float, float], float] | Callable[[float, float, float], float],
        alpha_cs: float | tuple[float, float],
        constraint: Literal['zeta', 'cPlus'] | None,
        ics: tuple[float, float, float],
        method: str,
    ) -> tuple:
        """
        Solve coupled ODE system for `c⁺(t), c⁻(t), s⁺(t)`
        """
        if isinstance(alpha_cs, (float, int)):
            alpha_cs = (alpha_cs, alpha_cs)
        alpha_c, alpha_s = alpha_cs
        if constraint is None:
            rhs = lambda _, y, *args: LateTimeFormulae.ode_rhs(*y, *args)
            args = (Da, epsilon, zeta0, sr, cr, alpha_c, alpha_s, flux)
        elif constraint == 'zeta':
            rhs = lambda _, y, *args: LateTimeFormulae.ode_rhs_zeta_constraint(*y, *args)
            ics = (ics[0], ics[2])
            args = (Da, epsilon, zeta0, sr, cr, alpha_c, alpha_s, flux)
        elif constraint == 'cPlus':
            rhs = lambda _, y, *args: LateTimeFormulae.ode_rhs_cPlus_constraint(*y, *args)
            ics_cPlus, ics_cMinus, ics_sPlus = ics
            ics = (ics_cMinus, )
            args = (ics_cPlus, ics_sPlus, epsilon, zeta0, sr, cr, flux)
        else:
            raise ValueError(constraint)

        solution = solve_ivp(
            rhs, 
            (t[0], t[-1]), 
            ics, 
            method, 
            t, 
            args=args,
        )
        if solution.success:
            return solution
        else:
            raise RuntimeError(f'{solution.message}')
    
    @staticmethod
    def ode_rhs(
        cPlus,
        cMinus,
        sPlus,
        Da,
        epsilon,
        zeta0,
        sr, 
        cr,
        alpha_c,
        alpha_s,
        flux,
    ):
        """
        `dc⁺(t) / dt = -F(c⁺, c⁻) / (1 - ζ(c⁺, c⁻, s⁺)) + Da Σ(s⁺, c⁺)` / (1 - s⁺)`

        `dc⁻(t) / dt = F(c⁺, c⁻) / (1 - ζ(c⁺, c⁻, s⁺))`
        
        `ds⁺(t) / dt = -εDa Σ(s⁺, c⁺)`
        """
        zeta = LateTimeFormulae.zeta(cPlus, cMinus, sPlus, epsilon, zeta0, sr, cr)
        r = LateTimeFormulae.Sigma(cPlus, sPlus)
        if arity(flux) == 2:
            f = flux(cPlus, cMinus)
        else:
            f = flux(cPlus, cMinus, sPlus)
        dcPlus = -f / (1 - zeta) + alpha_c * Da * r / (1 - sPlus)
        dcMinus = f / zeta
        dsPlus = -epsilon * alpha_s * Da * sPlus *(1 - cPlus)
        return [dcPlus, dcMinus, dsPlus]
    
    @staticmethod
    def ode_rhs_zeta_constraint(
        cPlus,
        sPlus,
        Da,
        epsilon,
        zeta0,
        sr, 
        cr,
        alpha_c,
        alpha_s,
        flux,
    ):
        """
        `dc⁺(t) / dt = -F(c⁺, c⁻) / (1 - ζ₀) + Da Σ(s⁺, c⁺)` / (1 - s⁺)`
        
        `ds⁺(t) / dt = -εDa Σ(s⁺, c⁺)`
        """
        cMinus = LateTimeFormulae.cMinus_from_zeta_constaint(cPlus, sPlus, epsilon, zeta0, sr, cr)
        r = LateTimeFormulae.Sigma(cPlus, sPlus)
        f = flux(cPlus, cMinus)
        dcPlus = -f / (1 - zeta0) + alpha_c * Da * r / (1 - sPlus)
        dsPlus = -epsilon * alpha_s * Da * sPlus *(1 - cPlus)
        return [dcPlus, dsPlus]
    
    @staticmethod
    def cMinus_from_zeta_constaint(cPlus, sPlus, epsilon, zeta0, sr, cr):
        m0_per_vol = LateTimeFormulae.m0(epsilon, zeta0, sr, cr, 1)
        cMinus = m0_per_vol
        cMinus += -(1 - zeta0) * (1 - sPlus) * cPlus 
        cMinus += -(1 - zeta0) * sPlus / epsilon
        cMinus *= (1 / zeta0)
        return cMinus
    
    @staticmethod
    def ode_rhs_cPlus_constraint(
        cPlus,
        cMinus,
        sPlus,
        epsilon,
        zeta0,
        sr, 
        cr,
        flux,
    ):
        zeta = LateTimeFormulae.zeta(cPlus, cMinus, sPlus, epsilon, zeta0, sr, cr)
        f = flux(cPlus, cMinus)
        dcMinus = f / zeta
        return [dcMinus]
    
    @staticmethod
    def Sigma(
        cPlus,
        sPlus,
    ):
        """
        `Σ(s⁺, c⁺) = s⁺(1 - c⁺)`
        """
        return sPlus * (1 - cPlus)

    @staticmethod
    def zeta(
        cPlus,
        cMinus,
        sPlus,
        epsilon,
        zeta0,
        sr,
        cr,
    ):
        """
        `ζ(c⁺, c⁻, s⁺) = ... / ...` horizontally-averaged interface height
        """ 
        m0_per_vol = LateTimeFormulae.m0(epsilon, zeta0, sr, cr, 1)
        numerator = m0_per_vol - sPlus / epsilon - (1 - sPlus) * cPlus
        denom = cMinus - sPlus/epsilon - (1 - sPlus) * cPlus
        return numerator / denom
    
    @staticmethod
    def mD(
        cPlus,
        cMinus,
        sPlus,
        epsilon,
        zeta0,
        sr,
        cr,
        vol: float,
    ):
        """
        `mᴰ(c⁺, c⁻, s⁺) = ...` dissolved mass
        """
        zeta = LateTimeFormulae.zeta(cPlus, cMinus, sPlus, epsilon, zeta0, sr, cr)
        m_per_vol = (1 - zeta) * (1 - sPlus) * cPlus + zeta * cMinus
        return m_per_vol * vol

    @staticmethod
    def mC(
        cPlus,
        cMinus,
        sPlus,
        epsilon,
        zeta0,
        sr,
        cr,
        vol: float,
    ):
        """
        `mᶜ(c⁺, c⁻, s⁺) = ...` capillary-trapped mass
        """
        zeta = LateTimeFormulae.zeta(cPlus, cMinus, sPlus, epsilon, zeta0, sr, cr)
        m_per_vol = (1 - zeta) * sPlus / epsilon 
        return m_per_vol * vol
    
    @staticmethod
    def m_sum(
        cPlus,
        cMinus,
        sPlus,
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
            LateTimeFormulae.mD(cPlus, cMinus, sPlus, epsilon, zeta0, sr, cr, vol),
            LateTimeFormulae.mC(cPlus, cMinus, sPlus, epsilon, zeta0, sr, cr, vol),
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
S_PLUS = TypeVar('S_PLUS')
FLUX = TypeVar('FLUX')
@dataclass
class LateTimeModel(Model):
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
    alpha_cs: float | tuple[float, float]
    flux: Callable[[C_PLUS, C_MINUS], FLUX] | Callable[[C_PLUS, C_MINUS, S_PLUS], FLUX]
    constraint: Literal['zeta', 'cPlus'] | None = None
    ics: tuple[float, float, float] | None = None
    method: str = 'RK45'

    def __post_init__(self):
        if self.ics is None:
            self.ics = (self.cr, 0.0, self.sr)
        self._ode_solution = self.call(LateTimeFormulae.ode_solve)

    @property
    def cPlus(self) -> np.ndarray:
        """
        `c⁺(t)` upper concentration
        """
        if self.constraint is None:
            return self._ode_solution.y[0]
        elif self.constraint == 'zeta':
            return self._ode_solution.y[0]
        elif self.constraint == 'cPlus':
            return self.ics[0]
        else:
            raise ValueError(self.constraint)
    
    @property
    def cMinus(self) -> np.ndarray:
        """
        `c⁻(t)`lower concentration
        """
        if self.constraint is None:
            return self._ode_solution.y[1]
        elif self.constraint == 'zeta':
            return self.call(LateTimeFormulae.cMinus_from_zeta_constaint)
        elif self.constraint == 'cPlus':
            self._ode_solution.y[0]
        else:
            raise ValueError(self.constraint)
    
    @property
    def sPlus(self) -> np.ndarray:
        """
        s⁺(t)` upper saturation
        """
        if self.constraint is None:
            return self._ode_solution.y[2]
        elif self.constraint == 'zeta':
            return self._ode_solution.y[1]
        elif self.constraint == 'cPlus':
            return self.ics[2]
        else:
            raise ValueError(self.constraint)
    
    @property
    def zeta(self) -> float | np.ndarray:
        """
        `ζ(t)` horizontally-averaged interface height 
        """
        if self.constraint:
            return self.zeta0
        else:
            return self.call(LateTimeFormulae.zeta)
    
    @property
    def mC(self) -> np.ndarray:
        """
        `mᶜ(t)` capillary-trapped mass
        """
        return self.call(LateTimeFormulae.mC)
    
    @property
    def mD(self) -> np.ndarray:
        """
        `mᴰ(t)` dissolved mass 
        """
        return self.call(LateTimeFormulae.mD)
    
    @property
    def f(self) -> np.ndarray:
        """
        `F(t) = F(c⁺(t), c⁻(t))` flux 
        """
        if arity(self.flux) == 2:
            return self.flux(self.cPlus, self.cMinus)
        else:
            return self.flux(self.cPlus, self.cMinus, self.sPlus)


TIME = TypeVar('TIME')
C_PLUS = TypeVar('C_PLUS')
C_MINUS = TypeVar('C_MINUS')
FLUX = TypeVar('FLUX')
def late_model_penalties(
    attr: str,
    model_factory: Callable[[TIME, Callable[[C_PLUS, C_MINUS], FLUX]], LateTimeModel],
    flux_factory: Callable[..., Callable[[C_PLUS, C_MINUS], FLUX]],
    flux_params: Iterable[tuple[float, ...] | float],
    dns_series: Iterable[float],
    dns_time_series: Iterable[float],
    penalty_metric: Callable[[Iterable[float], Iterable[float]], float] | None = None,
) -> list[float]:

    if penalty_metric is None:
        penalty_metric = lambda mdl, dns: np.sqrt(np.mean(mdl - dns))

    pens = []
    for prm in flux_params:
        if not isinstance(prm, tuple):
            prm = (prm, )
        model = model_factory(dns_time_series, flux_factory(prm))
        model_series = getattr(model, attr)
        pens.append(penalty_metric(model_series, dns_series))

    return pens


def find_optimal_late_model(
    attr: str,
    model_factory: Callable[[TIME, Callable[[C_PLUS, C_MINUS], FLUX]], LateTimeModel],
    flux_factory: Callable[..., Callable[[C_PLUS, C_MINUS], FLUX]],
    flux_params: Iterable[tuple[float, ...] | float],
    dns_series: Iterable[float],
    dns_time_series: Iterable[float],
    penalty_metric: Callable | None = None,
) -> tuple[float, ...] | float:
    penalties = late_model_penalties(
        attr,
        model_factory, 
        flux_factory, 
        flux_params, 
        dns_series, 
        dns_time_series, 
        penalty_metric,
    )
    index = penalties.index(min(penalties))
    return flux_params[index]
