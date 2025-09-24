from typing import Callable, Iterable
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar


class ExprSVS:

    @staticmethod
    def cb(
        y: float, 
        t: float, 
        Ra: float, 
        h0: float, 
        Lmbda: float,
    ) -> float:
        eigenvalues = ExprSVS.eigenvalues()
        c = 1
        for lmbda in eigenvalues:
            Cn = ExprSVS.Cn
            Yn = ExprSVS.Yn(y, h0, Lmbda, )
            c -= Cn * np.exp(-lmbda * t / Ra) * Yn(y, h0, Lmbda, lmbda)

        return c

    @staticmethod
    def Yn(
        y: float, 
        h0: float, 
        C: float, 
        lmbda: float,
    ) -> float:
        if y > h0:
            Bn = ExprSVS.Bn_plus(lmbda, C, h0)
        else:
            Bn = ExprSVS.Bn_minus(lmbda, C, h0)

        if lmbda > C:
            if y > h0:
                cfunc = np.cos
                sfunc = np.sin
                arg = lmbda - C
            else:
                cfunc = np.cos
                sfunc = np.sin
                arg = lmbda
        elif lmbda < C and lmbda > 0:
            if y > h0:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = C - lmbda
            else:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = lmbda
        else:
            if y > h0:
                cfunc = np.cosh
                sfunc = np.sinh
                arg = C + np.abs(lmbda)
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
        c0: Callable,
        s0: Callable,
        h0: float, 
        c: float, 
        lmbda: float,
    ) -> float:
        Yn = ExprSVS.Yn
        numer = lambda y: (1 - s0(y)) * (1 - c0(y)) * Yn(y, h0, c, lmbda)
        denom = lambda y: (1 - s0(y)) * Yn(y, h0, c, lmbda) ** 2
        y_interval = (0.0, 1.0)
        return quad(numer, *y_interval)[0] / quad(denom, *y_interval)[0]

    @staticmethod
    def eigenvalues(
    ) -> list[float]:
        ...

    @staticmethod
    def 


@dataclass
class SVS:
    """
    separation of variables solution
    `ε -> 0` and initial condition `s₀(y) = sr·H(y - h₀)` 
    """
    t: Iterable[float]
    y: np.ndarray
    Ra: float
    Da: float
    sr: float
    
    @property
    def C(self):
        return self.Ra * self.Da * self.sr / (1 - self.sr)

    @property
    def c(self) -> list[np.ndarray]:
        series = []
        for t in self.t:
            ct = ExprSVS.cb()
        return 



        


class DaInftyEpsZeroBaseState:
    """
    `Da -> ∞`, `ε -> 0` and initial condition `s₀(y) = sr·H(y - h₀)` 
    """
    
    def __init__(self, t, y, Ra, Da, Sr, epsilon):
        super.__init__()

    @staticmethod
    def concentration_upper(t, Da, Sr):
        if not isinstance(t, np.ndarray):
            t = np.array(t)
        return 1 - np.exp(-Da * Sr * t / (1 - Sr))
    

class DaInftyDiffusiveModel:

    ...
    
