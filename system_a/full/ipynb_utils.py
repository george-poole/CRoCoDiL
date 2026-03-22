import numpy as np

as_int_if_poss = lambda p: (
    int(p) if isinstance(p, float) and float.is_integer(p) else p

)

when_condition = lambda u, t, op: np.array(t)[np.where(op(np.array(u)))[0]]
when_geq = lambda u, t, val: when_condition(u, t, lambda _: _ >= val)
when_geq_first = lambda u, t, val: (
    when_geq(u, t, val)[0] if len(when_geq(u, t, val)) else np.inf
)