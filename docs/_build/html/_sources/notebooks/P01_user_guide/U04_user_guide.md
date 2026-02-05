# User guide

## Installation

...

## Creating convection-reaction systems

User-defined constitutive relations, boundary conditions and initial conditions may be prescribed to investigate novel models in CRoCoDiL. Solutal or thermal transport may be 'switched off' by setting the solutal or thermal dispersion relation to `None`. Porosity evolution may be 'switched off' by setting the evolution number to `None`, or equivalently $\varepsilon=0$.

## Best practices for advection-dominated systems

### Beware sharp gradients
...

### Resolving boundary layers

### Timestep contraints

CFL, reactive, diffusive

### Finite difference discretizations

### SUPG stabilization

pros: computational cost, cons: ambiguous $\tau$

### DG formulation

pros: sharp interfaces, cons: computational cost