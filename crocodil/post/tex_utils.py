class TeX:
    X = '$x$'
    """
    x-coordinate
    """
    Y = '$y$'
    """
    y-coordinate
    """
    T = '$t$'
    """
    time
    """
    MIN_X = staticmethod(lambda u: f'$\\min_{{\mathbf{{x}}}}({u})$')
    """
    `minₓ(u)`
    """
    MAX_X = staticmethod(lambda u: f'$\\max_{{\mathbf{{x}}}}({u})$')
    """
    `maxₓ(u)`
    """
    ABS_MIN_X = staticmethod(lambda u: f'$\\min_{{\mathbf{{x}}}}|\mathbf{{{u}}}|$')
    """
    `minₓ|𝐮|`
    """
    ABS_MAX_X = staticmethod(lambda u: f'$\\max_{{\mathbf{{x}}}}|\mathbf{{{u}}}|$')
    """
    `maxₓ|𝐮|`
    """
    BRAKET = lambda a: f'$\langle {a}\\rangle$'
    """
    `⟨u⟩`
    """
