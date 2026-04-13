from typing import Callable, Any, ParamSpec, TypeVar, Iterable
from inspect import signature
from functools import partial


P = ParamSpec('P')
R = TypeVar('R')
class Model:

    def _get_kws(
        self,
        clbl: Callable,
    ) -> dict[str, Any]:
        return {k: getattr(self, k) for k in signature(clbl).parameters} #TODO strict vs non-strict?
    
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
        exclude: Iterable[str] = (),
        **kws: Any,
    ) -> Callable[P, R]:
        _kws = self._get_kws(clbl)
        _kws = {k: v for k, v in _kws.items() if not k in exclude}
        _kws.update(kws)
        return partial(clbl, **_kws)