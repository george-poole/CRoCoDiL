from typing import Callable, Iterable
from functools import lru_cache, cached_property
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve


class ExprSeparableSolution:

    @staticmethod
    def c(
        y: float, 
        t: float, 
        Lmbda: float,
        Ra: float,
        h0: float,
        c0: Callable,
        s0: Callable,  
        cutoff: tuple[float, float] | None = None
    ) -> float:
        """
        `c(y, t)`
        """
        eigenvalues = ExprSeparableSolution.eigenvalues(Lmbda, h0, cutoff)
        c = 1
        for lmbda in eigenvalues:
            Cn = ExprSeparableSolution.Cn(c0, s0, h0, Lmbda, lmbda)
            Yn = ExprSeparableSolution.Yn(y, h0, Lmbda, lmbda)
            c -= Cn * np.exp(-lmbda * t / Ra) * Yn(y, h0, Lmbda, lmbda)
        return c

    @staticmethod
    def c_minus():
        """
        `Da -> ∞`
        """
        ...

    @staticmethod
    def Yn(
        y: float, 
        lmbda: float, 
        Lmbda: float,
        h0: float, 
    ) -> float:
        if y > h0:
            Bn = ExprSeparableSolution.Bn_plus(lmbda, Lmbda, h0)
        else:
            Bn = ExprSeparableSolution.Bn_minus(lmbda, Lmbda, h0)

        if lmbda > Lmbda:
            if y > h0:
                cfunc = np.cos
                sfunc = np.sin
                arg = lmbda - Lmbda
            else:
                cfunc = np.cos
                sfunc = np.sin
                arg = lmbda
        elif lmbda < Lmbda and lmbda > 0:
            if y > h0:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = Lmbda - lmbda
            else:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = lmbda
        else:
            if y > h0:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = Lmbda + np.abs(lmbda)
            else:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = np.abs(lmbda)

        arg = (y - h0) * np.sqrt(arg)
        return cfunc(arg) + Bn * sfunc(arg)
    
    @staticmethod
    def Bn_plus(lmbda, Lmbda, h0):
        if lmbda > Lmbda:
            arg = lmbda - Lmbda
            func = np.tan
        elif lmbda < Lmbda and lmbda > 0:
            arg = Lmbda - lmbda
            func = np.tanh
        else:
            arg = Lmbda + np.abs(lmbda)
            func = -np.tanh

        arg = (1 - h0) * np.sqrt(arg)
        return func(arg)      

    @staticmethod
    def Bn_minus(lmbda, Lmbda, h0):
        if lmbda > Lmbda:
            arg = lmbda
            func = -np.tan
        elif lmbda < Lmbda and lmbda > 0:
            arg = lmbda
            func = -np.tan
        else:
            arg = np.abs(lmbda)
            func = np.tanh

        arg = h0 * np.sqrt(arg)
        return func(arg)
    
    @staticmethod
    def Cn(
        lmbda: float,
        Lmbda: float,
        h0: float,  
        s0: Callable,
        c0: Callable,
    ) -> float:
        Yn = ExprSeparableSolution.Yn
        numer = lambda y: (1 - s0(y)) * (1 - c0(y)) * Yn(y, h0, Lmbda, lmbda)
        denom = lambda y: (1 - s0(y)) * Yn(y, h0, Lmbda, lmbda) ** 2
        y_interval = (0.0, 1.0)
        return quad(numer, *y_interval)[0] / quad(denom, *y_interval)[0]

    @staticmethod
    def eigenvalues(
        Lmbda: float,
        h0: float,
        cutoff: tuple[float, float] | None = None,
    ) -> np.ndarray:
        return ExprSeparableSolution._eigenvalues(Lmbda, h0, cutoff)
        
    @staticmethod
    @lru_cache
    def _eigenvalues(
        Lmbda: float,
        h0: float,
        cutoff: tuple[float, float] | None = None,
    ) -> np.ndarray:
        func = lambda x: ExprSeparableSolution.characteristic(x, Lmbda, h0)
        roots = fsolve(func, full_output=True)
        return roots

    @staticmethod
    def characteristic(
        lmbda: float,
        Lmbda: float,
        h0: float
    ) -> float:
        if lmbda > Lmbda:
            nfunc = np.tan
            dfunc = np.tan
            arg = lmbda - Lmbda
            coeff = -np.sqrt(arg)
        elif lmbda < Lmbda and lmbda > 0:
            nfunc = np.tanh
            dfunc = np.tan
            arg =  Lmbda - lmbda
            coeff = np.sqrt(arg)
        else:
            nfunc = np.tanh
            dfunc = np.tanh
            arg = Lmbda + np.abs(lmbda)
            coeff = -np.sqrt(arg)

        arg = (1 - h0) * np.sqrt(arg)
        return np.sqrt(np.abs(lmbda)) - coeff * nfunc(arg) / dfunc(h0 * np.sqrt(np.abs(lmbda)))



@dataclass
class SeparableSolution:
    """
    separation of variables solution
    `ε = 0` \\
    `s₀(y) = sr·H(y - h₀)` \\
    `c₀(y) = cr·H(y - h₀)`  
    """
    t: Iterable[float]
    y: Iterable[float]
    Ra: float
    Da: float
    h0: float
    sr: float
    cr: float

    expr = ExprSeparableSolution
    
    @property
    def Lmbda(self):
        """
        `Λ = Ra·Da·sr / (1 - sr)`
        """
        return self.Ra * self.Da * self.sr / (1 - self.sr)

    @cached_property
    def c(
        self,
    ) -> list[np.ndarray]:
        c0 = ...
        s0 = ...
    
        _c = []
        for t in self.t:
            _ct = []
            for y in self.y:
                _cyt = self.expr.c(y, t, self.Lmbda, self.Ra, self.h0, c0, s0)
                _ct.append(_cyt)
            _c.append(np.array(_ct))

        return _c


@dataclass
class SimilaritySolution:
    t: Iterable[float]
    y: Iterable[float]
    Ra: float
    Da: float
    h0: float
    sr: float
    cr: float

    def c():
        ...

    @property
    def c_early_time(
        self,
    ) -> list[float]:
        """
        `Da -> ∞` and `t << O(1/Da)` 
        """
        
        return [1 - np.exp(-self.Da * self.sr * t / (1 - self.sr)) for t in self.t]
