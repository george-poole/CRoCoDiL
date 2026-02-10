# Best practices

## Mesh resolution

critical wavelength

$$\lambda_{\text{c}} = \frac{90 L_y}{Ra}$$

$$2\min_{\textbf{x}}h(\textbf{x}) \leq\lambda_{\text{c}}$$

$h=\Delta x = \frac{L_x}{N_x}$ for a uniform Cartesian mesh.

$$
\begin{align*}
2\Delta x & \leq\lambda_{\text{c}} \\
\frac{2L_x}{N_x} & \leq\frac{90 L_y}{Ra} \\
\frac{\mathcal{A}Ra}{45} & \leq N_x \\
\end{align*}
$$


### Boundary layers

Dirichlet boundary conditions on the concentration or temperature can lead to boundary layers across which the concentration or temperature adjusts from its value in the bulk of the domain to its Dirichlet value.


## Adaptive timestep

$$
\begin{align*}
\Delta t &\leq \min\{\Delta t_{\textbf{u}}, \Delta t_{\mathsf{D}} ,\Delta t_{\Sigma}\} \\
\Delta t_{\textbf{u}} &= \min_{\textbf{x}}\left(\frac{h}{|\textbf{u}|}\right)  \\
\Delta t_{\mathsf{D}} & = \min_{\textbf{x}}\left(\frac{h^2}{|\mathsf{D}|}\right) \\
\Delta t_{\Sigma} &= \min_{\textbf{x}}\left(\frac{1}{|\Sigma|}\right) \\
\end{align*}
$$

$$
\Delta t = \min\{C_{\textbf{u}}\Delta t_{\textbf{u}}, C_{\mathsf{D}}\Delta t_{\mathsf{D}}, C_{\Sigma}\Delta t_{\Sigma}\}
$$



## Time discretization


| Finite difference operator | Recommended | Description |
| -------- | ------- | ------- | 
| $\mathcal{D}_\phi$ | $\text{AB}_1$ | forward Euler |  
| $\mathcal{D}_{\textbf{u},\theta}$ | $\text{AB}_2\circ\text{CN}$ | second-order Adams-Bashforth on $\textbf{u}$ and Crank-Nicolson on $\theta$ | 
| $\mathcal{D}_{\textbf{u},c}$ | $\text{AB}_2\circ\text{CN}$ | second-order Adams-Bashforth on $\textbf{u}$ and Crank-Nicolson on $c$ | 
| $\mathcal{D}_{\mathsf{G}, \theta}$ | $\text{AB}_1\circ\text{CN}$ | forward Euler on $\mathsf{G}$ and Crank-Nicolson on $\theta$ | 
| $\mathcal{D}_{\mathsf{D}, c}$ | $\text{AB}_1\circ\text{CN}$ | forward Euler on $\mathsf{D}$ and Crank-Nicolson on $c$ |
| $\mathcal{D}_{R,c}$ | $\text{AB}_1\circ\text{AM}_1$ | forward Euler on $R(s)$ and backward Euler on $c$|  
| $\mathcal{D}_{J}$ | $\text{AB}_1$ | forward Euler on $J(s)$ |  
| $\mathcal{D}_{\Sigma, s}$ | $\text{AB}_1\circ\text{AM}_1$ | forward Euler on $s$ and backward Euler on $c$ | 

For mass conservation, ideally choose finite difference operators such that

$$\mathcal{D}_{R,c}(Rc) + \mathcal{D}_J(J) = \mathcal{D}_{\Sigma, s}(\Sigma)$$


## Stabilization of the advection-diffusion-reaction equation

| Finite element method | Pros | Cons |
| -------- | ------- | ------- | 
| unstabilized continuous Galerkin | simplest formulation, least computational cost  | unstable to small diffusivities and large velocities | 
| SUPG-stabilized continuous Galerkin | simple formulation, low computational cost | can introduce artificial and crosswind diffusivities, ambiguous choice of parameter | 
| discontinuous Galerkin | stable to sharp gradients and discontinuities  |  intricate formulation, higher computational cost |

## Linear algebra options

Significant computional speedups can be achieved with a suitable choice of Krylov subspace  (`ksp_type`) and preconditioner (`pc_type`), informed by the properties of matrix operator resulting from the finite element discretization.

### Darcy equations

#### Mixed formulation

If no natural boundary conditions on the pressure are prescibed ($\partial\Omega_{\text{N}}=\varnothing\iff\partial\Omega_{\text{E}}=\partial\Omega$), then the linear algebra options need to be configured to handle the nullspace airising from the invariance of the equations under the addition of a constant to the pressure.

``` python
pc_factor_mat_solver_type = 'mumps'
```

#### Streamfunction formulation

The matrix is (in theory) symmetric. Interpolation, as opposed to projection, of the streamfunction gradients to obtain the velocity is the recommended method.

```python
psi_petsc = OptionsPETSc('cg', 'hypre')
u_petsc = None
flow_petsc = (psi_petsc, u_petsc)
```

### Advection-diffusion-reaction equation

The matrix is not symmetric

```python
c_petsc = OptionsPETSc('gmres', 'ilu')
theta_petsc = OptionsPETSc('gmres', 'ilu')
```