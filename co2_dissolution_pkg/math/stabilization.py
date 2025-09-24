from dolfinx.mesh import Mesh
from dolfinx.fem import Constant, Function
from ufl import (as_vector, conditional, as_vector, sqrt, conditional, lt, tanh)
from ufl.geometry import CellDiameter, MinCellEdgeLength, Circumradius
from ufl.core.expr import Expr

from lucifex.fdm import DT, FiniteDifference
from lucifex.fdm.ufl_operators import inner, grad


def supg_velocity(
    phi: Expr | Function | Constant,
    u: Function, 
    Ra: Constant,
    D_adv: FiniteDifference | tuple[FiniteDifference, FiniteDifference],
    D_diff: FiniteDifference,
):
    DiscretizationError = RuntimeError(
        'SUPG stabilization requires advection term to be implicit in concentration',
    )
    match D_adv:
        case D_adv_u, D_adv_c:
            if D_adv_c.is_explicit:
                raise DiscretizationError
            u_coeff = D_adv_u(u) * D_adv_c.explicit_coeff
        case D_adv:
            if D_adv.is_explicit:
                raise DiscretizationError
            u_coeff = u[1] * D_adv.explicit_coeff

    u_eff = (1/phi) * u_coeff \
                - (1/Ra) * (1/phi) * grad(phi) * D_diff.explicit_coeff
    
    return u_eff


def supg_diffusivity(
    Ra: Constant,
    D_diff: FiniteDifference,
):
    return (1/Ra) * D_diff.explicit_coeff


def supg_reaction(
    dt: Constant,
    phi: Expr | Function | Constant,
    Da: Constant,
    D_reac: FiniteDifference | tuple[FiniteDifference, ...],
    r_factor: float = 0.1, # order of magnitude estimate
):
    # NOTE assumes reaction such that `R(s,c) = R(s)R(c)` and `R(s,0) = R(s)`
    # r_func, r_args = r
    # s = r_args[0]
    # match D_reac:
    #     case D_reac_s, D_reac_c:
    #         r_coeff = r(D_reac_s[False](s), 0) * D_reac_c.explicit_coeff 
    #     case D_reac:
    #         r_coeff = r_func(s[1], 0) * D_reac.explicit_coeff 

    if isinstance(D_reac, tuple):
        c_trial_coeff = r_factor * D_reac[1].explicit_coeff 
    else:
        c_trial_coeff = r_factor * D_reac.explicit_coeff 

    r_eff = (1/dt) * DT.explicit_coeff + Da * (1/phi) * c_trial_coeff
    return r_eff


def supg_tau(
    method: str,
    mesh: Mesh,
    u_eff: Expr,
    d_eff: Expr,
    r_eff: Expr = 0,
    cell_size: type[MinCellEdgeLength] | type[Circumradius] = CellDiameter,
) -> Expr:
    
    h = cell_size(mesh)

    if method == 'codina':
        return tau_codina(u_eff, d_eff, r_eff, h)
    elif method == 'shakib':
        return tau_shakib(u_eff, d_eff, r_eff, h)
    elif method == 'upwind':
        return tau_upwind(u_eff, d_eff, h)
    elif method == 'upwind_xy':
        # TODO how to get hx, hy in ufl?
        return tau_upwind_xy(u_eff, d_eff, h, h)
    else:
        raise RuntimeError(f'Invalid SUPG choice {method}.')


def tau_codina(U_eff, D_eff, R_eff, h) -> Expr:
    u = sqrt(inner(U_eff, U_eff))
    return ((2 * u / h) + (4 * D_eff / h**2) + R_eff) ** (-1) 


def tau_shakib(U_eff, D_eff, R_eff, h) -> Expr:
    u = sqrt(inner(U_eff, U_eff))
    return ((2 * u / h)**2 + 9 * (4 * D_eff / h**2)**2 + R_eff**2) ** (-0.5) 


def tau_upwind(U_eff, D_eff, h) -> Expr:
    u = sqrt(inner(U_eff, U_eff))
    return (0.5 * h / u) * xi(Pe(u, h, D_eff))


def tau_upwind_xy(U_eff, D_eff, hx, hy) -> Expr:
    ex = as_vector((1, 0))
    ey = as_vector((0, 1))
    Ux = inner(ex, U_eff)
    Uy = inner(ey, U_eff)
    ux = sqrt(inner(Ux, Ux))
    uy = sqrt(inner(Uy, Uy))
    xi_x = xi(Pe(ux, hx, D_eff))
    xi_y = xi(Pe(uy, hy, D_eff))
    return 0.5 * (xi_x * ux * hx + xi_y * uy * hy) / inner(U_eff, U_eff)


def Pe(u, h, D):
    return 0.5 * u * h / D


def xi(Pe):
    return 1/tanh(Pe) - Pe


def xi_approx(Pe):
    return conditional(lt(Pe, 3), Pe / 3, 1)


# def edge_stabilization():
#     """
#     NOTE incompatible with quadrilateral cells
#     """
#     gamma = Constant(c.function_space.mesh, edge)
#     h = FacetArea(c.function_space.mesh) 
#     n = FacetNormal(c.function_space.mesh)
#     dim = c.function_space.mesh.geometry.dim
#     F_edge = (0.5 * gamma * h("+")**dim) * jump(inner(n, grad(v))) * jump(inner(n, grad(trialfunction(c)))) * dS
#     forms.append(F_edge)