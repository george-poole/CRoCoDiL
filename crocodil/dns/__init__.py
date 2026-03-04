from .generic import dns_generic
from .classic import (
    dns_darcy_evolving, dns_darcy_rayleigh_benard, dns_darcy_rayleigh_taylor,
    dns_darcy_thermosolutal,
)
from .diagnostic import mass_capillary, mass_dissolved, vertical_flux
from .theory import CONVECTION_REACTION_SCALINGS, CONVECTION_CONSTANTS