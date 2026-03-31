from typing import Any, Iterable
import numpy as np

pass_anything = lambda *_, **__: pass_anything

when_condition = lambda u, t, op: np.array(t)[np.where(op(np.array(u)))[0]]

when_geq = lambda u, t, val: when_condition(u, t, lambda _: _ >= val)

when_geq_first = lambda u, t, val, sentinel: (
    when_geq(u, t, val)[0] if len(when_geq(u, t, val)) else sentinel
)

def include_prm(prm: Any, targets: Iterable[Any] | None) -> bool:
    if targets is None:
        return True
    return prm in targets