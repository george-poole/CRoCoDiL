from lucifex.fdm import AB1
from lucifex.solver import (
    projection_solver, 
    interpolation_solver, bvp_solver, ivp_solver, ibvp_solver, 
    OptionsPETSc, BoundaryValueProblem, ProjectionProblem, InterpolationProblem)


from .transport import (darcy_streamfunction, streamfunction_velocity, darcy_incompressible, 
                  advection_diffusion_cg, advection_diffusion_dg, 
                  advection_diffusion_reaction_dg, advection_diffusion_reaction_cg,
                  evolution_expression, evolution)


def streamfunction_method(
    petsc: OptionsPETSc | tuple,
) -> bool:
    if isinstance(petsc, tuple):
        return True
    return False