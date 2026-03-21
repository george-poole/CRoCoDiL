from typing import Callable, Any, ParamSpec, TypeVar
from inspect import signature
from functools import partial


P = ParamSpec('P')
R = TypeVar('R')
class Model:

    def _get_kws(
        self,
        clbl: Callable,
    ) -> dict[str, Any]:
        return {k: getattr(self, k) for k in signature(clbl).parameters}
    
    def call(
        self,
        clbl: Callable[..., R],
        *args: Any,
        **kws: Any,
    ) -> R:
        _kws = self._get_kws(clbl)
        _kws.update(kws)
        return clbl(*args, **_kws)
    
    def create_partial(
        self,
        clbl: Callable[P, R],
        **kws: Any,
    ) -> Callable[..., R]:
        _kws = self._get_kws(clbl)
        _kws.update(kws)
        return partial(clbl, **_kws)