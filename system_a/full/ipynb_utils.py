import numpy as np

from lucifex.fem import GridFunction

pass_anything = lambda *_, **__: pass_anything

when_condition = lambda u, t, op: np.array(t)[np.where(op(np.array(u)))[0]]

when_geq = lambda u, t, val: when_condition(u, t, lambda _: _ >= val)

when_geq_first = lambda u, t, val: (
    when_geq(u, t, val)[0] if len(when_geq(u, t, val)) else np.inf
)


def correct_interpolation_error(
    s: GridFunction,
    
):
    ...