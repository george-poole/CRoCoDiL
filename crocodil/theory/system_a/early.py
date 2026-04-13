from typing import Callable, Iterable, Any
from functools import lru_cache, cached_property, partial
from dataclasses import dataclass, field

import numpy as np
import scipy.special as sp
from scipy.integrate import quad
from scipy.optimize import fsolve

from .model import Model


class EarlyTimeExactFormulae:

    @staticmethod
    def c(
        y: float, 
        t: float, 
        Lmbda: float,
        Di: float,
        zeta0: float,
        eigenvalues,
        coefficients,
    ) -> float:
        """
        `c(y, t) = 1 - Σ ...`
        """
        # eigenvalues = EarlyTimeExactFormulae.eigenvalues(
        #     Lmbda, zeta0, eigen_guesses, **kwargs
        # )
        # if eigen_min is not None:
        #     eigenvalues = [i for i in eigenvalues if i >= eigen_min]
        # if eigen_max is not None:
        #     eigenvalues = [i for i in eigenvalues if i <= eigen_max]

        c = 1.0
        for Ln, Cn in zip(eigenvalues, coefficients):
            Yn = EarlyTimeExactFormulae.Yn(y, Ln, Lmbda, zeta0)
            c -= Cn * np.exp(-Di * Ln * t) * Yn
        return c

    @staticmethod
    def Yn(
        y: float, 
        Ln: float, 
        Lmbda: float,
        zeta0: float, 
    ) -> float:
        if y > zeta0:
            Bn = EarlyTimeExactFormulae.Bn_plus(Ln, Lmbda, zeta0)
        else:
            Bn = EarlyTimeExactFormulae.Bn_minus(Ln, Lmbda, zeta0)

        if Ln > Lmbda:
            if y > zeta0:
                cfunc = np.cos
                sfunc = np.sin
                arg = Ln - Lmbda
            else:
                cfunc = np.cos
                sfunc = np.sin
                arg = Ln
        elif Ln < Lmbda and Ln > 0:
            if y > zeta0:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = Lmbda - Ln
            else:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = Ln
        else:
            if y > zeta0:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = Lmbda + np.abs(Ln)
            else:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = np.abs(Ln)

        arg = (y - zeta0) * np.sqrt(arg)
        return cfunc(arg) + Bn * sfunc(arg)
    
    @staticmethod
    def Bn_plus(Ln, Lmbda, zeta0):
        if Ln > Lmbda:
            arg = Ln - Lmbda
            func = np.tan
        elif Ln < Lmbda and Ln > 0:
            arg = Lmbda - Ln
            func = np.tanh
        else:
            arg = Lmbda + np.abs(Ln)
            func = -np.tanh

        arg = (1 - zeta0) * np.sqrt(arg)
        return func(arg)      

    @staticmethod
    def Bn_minus(Ln, Lmbda, zeta0):
        if Ln > Lmbda:
            arg = Ln
            func = lambda x: -np.tan(x)
        elif Ln < Lmbda and Ln > 0:
            arg = Ln
            func = lambda x: np.tan(x)
        else:
            arg = np.abs(Ln)
            func = np.tanh

        arg = zeta0 * np.sqrt(arg)
        return func(arg)
    
    @staticmethod
    def Cn(
        Ln: float,
        Lmbda: float,
        zeta0: float,  
        s0: Callable[[float], float],
        c0: Callable[[float], float],
        Ly: float | Iterable[float] ,
    ) -> float:
        Yn = partial(EarlyTimeExactFormulae.Yn, Ln=Ln, Lmbda=Lmbda, zeta0=zeta0)
        numer = lambda y: (1 - s0(y)) * (1 - c0(y)) * Yn(y)
        denom = lambda y: (1 - s0(y)) * Yn(y) ** 2
        if isinstance(Ly, float):
            Ly = (0.0, Ly)
        return quad(numer, min(Ly), max(Ly))[0] / quad(denom,  min(Ly), max(Ly))[0]

    @staticmethod
    def eigenvalues(
        Lmbda: float,
        zeta0: float,
        guesses: Iterable[float],
        **kwargs: Any,
    ) -> list[float]:
        return list(EarlyTimeExactFormulae._eigenvalues(Lmbda, zeta0, guesses, **kwargs))
        
    @staticmethod
    @lru_cache
    def _eigenvalues(
        Lmbda: float,
        zeta0: float,
        guesses: Iterable[float],
        **kwargs: Any,
    ) -> np.ndarray:
        
        root_func = partial(
            EarlyTimeExactFormulae.eigenvalue_characteristic, Lmbda=Lmbda, zeta0=zeta0,
        ) 
        eigens = []
        
        for guess in guesses:
            roots, _, ier, _ = fsolve(root_func, guess, full_output=True, **kwargs)
            if ier == 1:
                if not np.any(np.isclose(eigens, roots[0])):
                    eigens.append(roots[0])
                continue

        return np.sort(eigens)

    @staticmethod
    def eigenvalue_characteristic(
        Ln: float,
        Lmbda: float,
        zeta0: float
    ) -> float:
        if Ln > Lmbda:
            nfunc = np.tan
            dfunc = np.tan
            arg = Ln - Lmbda
            sgn = -1
        elif Ln < Lmbda and Ln > 0:
            nfunc = np.tanh
            dfunc = np.tan
            arg =  Lmbda - Ln
            sgn = 1
        else:
            nfunc = np.tanh
            dfunc = np.tanh
            arg = Lmbda + np.abs(Ln)
            sgn = -1


        lhs = np.sqrt(np.abs(Ln))
        rhs = sgn * np.sqrt(arg) * nfunc((1 - zeta0) * np.sqrt(arg))
        rhs *= 1 / dfunc(zeta0 * np.sqrt(np.abs(Ln)))

        return lhs - rhs
    
    @staticmethod
    def cLowerApprox():
        ...

    @staticmethod
    def Lmbda(Di, Ki, sr):
        """
        `Λ = Ki·sr / Di·(1 - sr)`
        """
        return (1 / Di) * Ki * sr / (1 - sr)

    @staticmethod
    def Pi(Ki, sr):
        """
        `Π = Di · Λ = Ki·sr / (1 - sr)`
        """
        return Ki * sr / (1 - sr)


@dataclass
class EarlyTimeExactModel(Model):
    """
    `ε = 0` \\
    `s₀(y) = sr·H(y - ζ₀)` \\
    `c₀(y) = cr·H(y - ζ₀)`  
    """
    t: Iterable[float]
    y: Iterable[float]
    Di: float
    Ki: float
    zeta0: float
    sr: float
    cr: float
    eigen_guesses: Iterable[float] = range(1, 100, 1),
    eigen_min: float | None = None
    eigen_max: float | None = None
    kwargs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        self.Ly = (min(self.y), max(self.y))
        self._c: dict[float, np.ndarray] = {}
    
    @property
    def Lmbda(self):
        return EarlyTimeExactFormulae.Lmbda(self.Di, self.Ki, self.sr)

    @property
    def Pi(self):
        return EarlyTimeExactFormulae.Pi(self.Ki, self.sr)


    @cached_property
    def c_series(self) -> list[np.ndarray]:
        return [self.c(ti) for ti in self.t]
    
    def c(
        self,
        t: float,
    ) -> np.ndarray:
        if t in self._c:
            return self._c[t]
        else:
            c_formula = partial(
                EarlyTimeExactFormulae.c, 
                Lmbda=self.Lmbda, Di=self.Di, zeta0=self.zeta0, 
                eigenvalues=self.eigenvalues,
                coefficients=self.coefficients,
                **self.kwargs,
            )
            # c_formula = self.create_partial(EarlyTimeExactFormulae.c)
            c_arr = np.array([c_formula(yi, t)for yi in self.y])
            self._c[t] = c_arr
            return c_arr
        
    def eigenfunction(self, n: int) -> np.ndarray:
        Ln = self.eigenvalues[n]
        return [EarlyTimeExactFormulae.Yn(i, Ln, self.Lmbda, self.zeta0) for i in self.y]
    
    @cached_property
    def eigenvalues(self) -> np.ndarray:
        return EarlyTimeExactFormulae.eigenvalues(
            self.Lmbda, self.zeta0, self.eigen_guesses,
        )
    
    @cached_property
    def coefficients(self) -> list[float]:
        s0 = lambda y: self.sr if y > self.zeta0 else 0
        c0 = lambda y: self.cr if y > self.zeta0 else 0
        # Cn = self.create_partial(EarlyTimeExactFormulae.Cn, s0=s0, c0=c0)
        return [EarlyTimeExactFormulae.Cn(i, self.Lmbda, self.zeta0, s0=s0, c0=c0, Ly=self.Ly) for i in self.eigenvalues]


class EarlyTimeSimilarityFormulae:

    Lmbda = EarlyTimeExactFormulae.Lmbda
    Pi = EarlyTimeExactFormulae.Pi
    
    @staticmethod
    def cPlus(
        t,
        Ki,
        sr,
        cr,
    ):
        Pi = EarlyTimeSimilarityFormulae.Pi(Ki, sr)
        return 1 - (1 - cr) * np.exp(-t * Pi)
    
    @staticmethod
    def cZeta(
        t,
        Di,
        Ki,
        zeta0,
        sr,
        cr,
        infty: bool,
    ):
        Pi = EarlyTimeSimilarityFormulae.Pi(Ki, sr)
        const = np.sqrt(np.pi * Pi)
        if not infty:
            Lmbda = EarlyTimeSimilarityFormulae.Lmbda(Di, Ki, sr)
            const *= np.tanh((1 - zeta0) * np.sqrt(Lmbda))
        return (np.sqrt(t) * const + cr) / (np.sqrt(t) * const + 1)

    @staticmethod
    def cUpper(
        t,
        y,
        Di,
        Ki,
        zeta0,
        sr,
        cr,
        infty: bool,
    ):
        cZeta = EarlyTimeSimilarityFormulae.cZeta(t, Di, Ki, zeta0, sr, cr, infty)
        Lmbda = EarlyTimeSimilarityFormulae.Lmbda(Di, Ki, sr)
        if infty:
            return 1 - (1 - cZeta) * np.exp(-(y - zeta0) * np.sqrt(Lmbda))
        else:
            num = np.cosh((1 - y) * np.sqrt(Lmbda))
            denom = np.cosh((1 - zeta0) * np.sqrt(Lmbda))
            return 1 - (1 - cZeta) * num / denom

    @staticmethod
    def cLower(
        t,
        y,
        Di,
        Ki,
        zeta0,
        sr,
        cr,
        infty: bool,
    ):
        cZeta = EarlyTimeSimilarityFormulae.cZeta(t, Di, Ki, zeta0, sr, cr, infty)
        return cZeta + (cZeta - cr) * sp.erf((y - zeta0) / (2 * np.sqrt(Di * t)))

    
@dataclass
class EarlyTimeSimilarityModel(Model):
    t: Iterable[float]
    y: Iterable[float]
    Di: float
    Ki: float
    zeta0: float
    sr: float
    cr: float
    infty: bool = False
    tol: float = 1e-10

    def __post_init__(self):
        self.t = [i if i > self.tol else self.tol for i in self.t]

    @cached_property
    def c(self) -> list[np.ndarray]:
        """
        `c(y,t)`
        """
        cUpper = self.create_partial(EarlyTimeSimilarityFormulae.cUpper)
        cLower = self.create_partial(EarlyTimeSimilarityFormulae.cLower)

        c_series = []
        for ti in self.t:
            ci = np.array(
                [cUpper(t=ti, y=yi) if yi > self.zeta0 else cLower(t=ti, y=yi) for yi in self.y]
            )
            c_series.append(ci)

        return c_series
    
    @cached_property
    def cZeta(self) -> list[float]:
        c = self.create_partial(EarlyTimeSimilarityFormulae.cZeta)
        return [c(t=ti)for ti in self.t]
    
    @cached_property
    def cPlus(self) -> list[float]:
        """
        `c(y > ζ, t << 1/Π)` 
        """
        c = self.create_partial(EarlyTimeSimilarityFormulae.cPlus)
        return [c(t=ti)for ti in self.t]
    
    @property
    def Pi(self):
        return EarlyTimeSimilarityFormulae.Pi(self.Ki, self.sr)
    
    @property
    def Lmbda(self):
        return EarlyTimeSimilarityFormulae.Lmbda(self.Di, self.Ki, self.sr)
