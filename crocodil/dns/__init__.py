from .generic import dns_generic
from .benchmark import (
    dns_darcy_evolving, dns_darcy_rayleigh_benard, dns_darcy_rayleigh_taylor,
    dns_darcy_thermosolutal,
)
from .utils import (
    CONVECTION_REACTION_SCALINGS, vertical_flux, critical_wavelength,
    critical_dx, critical_Nx,
    heaviside, rectangle_mesh_closure, mass_capillary_trapped, mass_dissolved
)