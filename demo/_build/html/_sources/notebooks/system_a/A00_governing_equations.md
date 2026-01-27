# Governing equations

## Dimensional governing equations

$\textbf{x}=(x,y)\in\Omega\sub\mathbb{R}^2$ or $\textbf{x}=(x,y,z)\in\Omega\sub\mathbb{R}^3$ spanned by unit vectors $\{\textbf{e}_x, \textbf{e}_y\}$ or $\{\textbf{e}_x, \textbf{e}_y, \textbf{e}_z\}$

$$
\begin{align}
\phi\frac{\partial c}{\partial t} + \textbf{u}\cdot\nabla c &= \nabla\cdot(\mathsf{D}(\phi, \textbf{u})\cdot\nabla c) + R(s, c, \theta) & \\
\phi\frac{\partial\theta}{\partial t} + \textbf{u}\cdot\nabla\theta &= \nabla\cdot(\mathsf{G}(\phi, \textbf{u})\cdot\nabla\theta) & \\
\nabla\cdot\textbf{u} &= 0 & \\
\textbf{u}&=-\frac{\mathsf{K}(\phi)}{\mu(c, \theta)}\cdot(\nabla p - \,\rho(c, \theta) g\,\textbf{e}_g) \\
\varrho\varphi(\textbf{x})\frac{\partial s}{\partial t} &= -R(s,c,\theta)
\end{align}
$$

$\phi(\textbf{x}, t)=\varphi(\textbf{x})(1-s(\textbf{x},t))$ is the effective porosity. The Boussinesq approximation has been used to neglect variations in density expect in the gravitational term and assume that $c\ll\rho_{\text{ref}}$ for some reference density $\rho_{\text{ref}}$. The porosity is assumed to be slowly-evolving and so incompressibility is also approximated.

## Non-dimensionalization

### Scalings

| Quantity | $\vert\textbf{x}\vert$ | $\vert\textbf{u}\vert$ | $t$ | $c$ | $\theta$ | $\rho g$ | $p$ | $\psi$ |
| -------- | ------- | ------- | ------- | ------- | ------- |  ------- |  ------- |  ------- | 
| **Scaling** | $\mathcal{L}$ | $\mathcal{U}$ |$\mathcal{T}$ | $\Delta c$ | $\Delta\theta$ | $g \Delta\rho$ | $\mu_{\text{ref}}\,\mathcal{U}\mathcal{L}/K_{\text{ref}}$ | $\mathcal{U}\mathcal{L}$ |

| $\mu$ | $\phi$ | $K$ | $\vert\mathsf{D}\vert$ | $\vert\mathsf{G}\vert$ | $R$ | $Q$ |
| ------- | ------- | ------- | ------- | ------- |  ------- | ------- |  
| $\mu_{\text{ref}}$ | $\phi_{\text{ref}}$ |$K_{\text{ref}}$ | $D_{\text{ref}}$ | $G_{\text{ref}}$ | $\Delta R$ | $\Delta Q$ |

### Dimensionless numbers

#### Abstract dimensionless numbers

$$
\begin{align*}
Ad&=\frac{\mathcal{U}\mathcal{T}}{\phi_{\text{ref}}\mathcal{L}} \\
Pe&=\frac{\phi_{\text{ref}}\mathcal{L}^2}{D_{\text{ref}}\mathcal{T}} \\
Ki&=\frac{\mathcal{T}\Delta R}{\phi_{\text{ref}}\Delta c} \\
Bu&=\frac{K_{\text{ref}}\,g\Delta\rho}{\mu_{\text{ref}}\,\mathcal{U}} \\
Xl&=\frac{\mathcal{L}}{\mathcal{L}_\Omega}
\end{align*} \\
$$

A particular choice of non-dimensionalization with $\mathcal{L}$, $\mathcal{U}$ and $\mathcal{T}$ chosen according to the details of the specified system (e.g. boundary conditions, dominant transport processes) will map the abstract dimensionless numbers $\{Ad, Pe, Ki, Bu, Xl\}$ onto combinations of the physical dimensionless numbers introduced below.

| Name | $\mathcal{L}$ | $\mathcal{U}$ |$ \mathcal{T}$ | $\{Ad, Pe, Ki, Bu, Xl\}$ | Examples |
| -------- | ------- | ------- | ------- | ------- | ------- |
| $\text{advective}$ | $\mathcal{L}_\Omega$  |  $K_{\text{ref}}\,g\Delta\rho/\mu_{\text{ref}}$  | $\phi_{\text{ref}}\mathcal{L}/\mathcal{U}$ | $\{1, Ra, Da, 1, 1\}$| [Hewitt et al. (2012)](https://link.aps.org/doi/10.1103/PhysRevLett.108.224503) |
| $\text{advective-diffusive}$ | $D_{\text{ref}}/\mathcal{U}$  |  $K_{\text{ref}}\,g\Delta\rho/\mu_{\text{ref}}$  | $\phi_{\text{ref}}\mathcal{L}/\mathcal{U}$ | $\{1, 1, Da/Ra, 1, Ra\}$| [Slim (2014)](https://www.cambridge.org/core/product/identifier/S0022112013006733/type/journal_article) |
| $\text{diffusive}$ | $\mathcal{L}_\Omega$  |  $D_{\text{ref}}/\mathcal{L}$  | $\phi_{\text{ref}}\mathcal{L}/\mathcal{U}$ | $\{1, 1, RaDa, Ra, 1\}$| [Ritchie \& Pritchard  (2011)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/natural-convection-and-the-evolution-of-a-reactive-porous-medium/71E5FB557F61CB9125E5B4E4EE9D828F) |
| $\text{reactive}$ | $\sqrt{D_{\text{ref}}\mathcal{T}/\phi_{\text{ref}}}$  |  $\phi_{\text{ref}}\mathcal{L}/\mathcal{T}$  | $\phi_{\text{ref}}\Delta c/\Delta R$  | $\{1, 1, 1, \sqrt{Ra/Da}, \sqrt{RaDa}\}$| [Kabbadj et al. (2025)](https://nlpc.ulb.be/pdf/25.Kabbadj_MATRIX.pdf) |


#### Physical dimensionless numbers

Rayleigh number (defined with respect to the transport of $c$ and domain length scale)
$$Ra=\frac{\mathcal{L}_\Omega K_{\text{ref}}g\Delta\rho}{\mu_{\text{ref}}D_{\text{ref}}}=\underbrace{\frac{K_{\text{ref}}\,g\Delta\rho}{\mu_{\text{ref}}}}_{\text{convective speed}} \big/ \underbrace{\frac{D_{\text{ref}}}{\mathcal{L}_\Omega}}_{\text{diffusive speed}}$$

$Da$ is the Damköhler number (defined with respect to the transport of $c$ and domain length scale)

$$Da=\frac{\mathcal{L}_\Omega \mu_{\text{ref}}\,\Delta R}{K_{\text{ref}}\,g\Delta\rho\Delta c} = \underbrace{\frac{\Delta R}{\Delta c}}_{\text{reaction rate}} \big/ \underbrace{\frac{K_{\text{ref}}\,g\Delta\rho}{\mathcal{L}_\Omega \mu_{\text{ref}}}}_{\text{convection rate}}$$

evolution number

$$\varepsilon = \frac{\Delta c}{\varrho}\ll1$$

Lewis number for the ratio of thermal to solutal diffusivity

$$Le=\frac{G_{\text{ref}}}{D_{\text{ref}}}$$


## Non-dimensional governing equations

### Mixed formulation

In strong form, the non-dimensionalized initial boundary value problem on the domain $\Omega\sub\mathbb{R}^d$ is
$$
\begin{align*}
&\text{Find} \\
&\text{$c(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$, } \\
&\text{$\theta(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$, } \\
&\text{$s(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$, } \\
&\text{$\textbf{u}(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}^d$ and $p(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$} \\
& \text{such that} \\
&\begin{cases}
\phi\frac{\partial c}{\partial t} + Ad\,\textbf{u}\cdot\nabla c = \frac{1}{Pe}\nabla\cdot(\mathsf{D}(\phi, \textbf{u})\cdot\nabla c) + KiR(s, c, \theta) & \\
\phi\frac{\partial\theta}{\partial t} + Ad\,\textbf{u}\cdot\nabla\theta = \frac{1}{LePe}\nabla\cdot(\mathsf{G}(\phi, \textbf{u})\cdot\nabla\theta) & \\
\nabla\cdot\textbf{u} = 0 & \\
\textbf{u}=-\frac{\mathsf{K}(\phi)}{\mu(c, \theta)}\cdot(\nabla p + Bu\,\rho(c, \theta) \,\textbf{e}_g) \\
\varphi(\textbf{x})\frac{\partial s}{\partial t} = -\varepsilon Ki R(s,c,\theta) & \forall(\textbf{x}, t)\in\Omega\times[0,\infty) \\
c(\textbf{x},t=0)=c_0 & \forall\textbf{x}\in\Omega \\
\theta(\textbf{x},t=0)=\theta_0 & \forall\textbf{x}\in\Omega \\
s(\textbf{x},t=0)=s_0 & \forall\textbf{x}\in\Omega \\
c=c_{\text{D}} & \forall(\textbf{x}, t)\in\partial\Omega_{\text{D}, c} \times [0,\infty] \\
\textbf{n}\cdot(\mathsf{D}\cdot\nabla c) = c_{\text{N}} & \forall(\textbf{x}, t)\in\partial\Omega_{\text{N}, c}
\times [0,\infty]~,~\partial\Omega_{\text{N}, c}=\partial\Omega/\partial\Omega_{\text{D}, c} \\
\theta=\theta_{\text{D}} & \forall (\textbf{x}, t)\in\partial\Omega_{\text{D}, \theta} \times [0,\infty] \\
\textbf{n}\cdot(\mathsf{G}\cdot\nabla \theta) = \theta_{\text{N}} & \forall(\textbf{x}, t)\in\partial\Omega_{\text{N}, \theta}
\times [0,\infty]~,~\partial\Omega_{\text{N}, \theta}=\partial\Omega/\partial\Omega_{\text{D}, \theta} \\
\textbf{n}\cdot\textbf{u} = u_{\text{E}} & \forall(\textbf{x}, t)\in\partial\Omega_{\text{E}} \times [0,\infty] \\
p = p_{\text{N}} & \forall(\textbf{x}, t)\in\partial\Omega_{\text{N}}\times [0,\infty]~,~\partial\Omega_{\text{N}}=\partial\Omega/\partial\Omega_{\text{E}}
\end{cases}~.
\end{align*}
$$

### Streamfunction formulation

In $d=2$ with $\partial\Omega_{\text{E}}=\partial\Omega\iff\Omega_{\text{N}}=\varnothing$, the streamfunction $\psi$ is solved for instead of the velocity $\textbf{u}$ and pressure $p$.

$$
\begin{align*}
&\text{Find} \\
&\text{$c(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$, } \\
&\text{$\theta(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$, } \\
&\text{$s(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$, } \\
&\text{$\psi(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$ and $\textbf{u}(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}^d$} \\
& \text{such that} \\
&\begin{cases}
\phi\frac{\partial c}{\partial t} + Ad\,\textbf{u}\cdot\nabla c = \frac{1}{Pe}\nabla\cdot(\mathsf{D}(\phi, \textbf{u})\cdot\nabla c) + KiR(s, c, \theta) & \\
\phi\frac{\partial\theta}{\partial t} + Ad\,\textbf{u}\cdot\nabla\theta = \frac{1}{LePe}\nabla\cdot(\mathsf{G}(\phi, \textbf{u})\cdot\nabla\theta) & \\
\nabla\cdot\bigg(\frac{\mu\mathsf{K}^{\mathsf{T}}\cdot\nabla\psi}{\det\mathsf{K}}\bigg) = -\frac{\partial(\rho\,\textbf{e}_g\cdot\textbf{e}_y)}{\partial x} + \frac{\partial(\rho\,\textbf{e}_g\cdot\textbf{e}_x)}{\partial y} \\
\textbf{u} = (-\frac{\partial\psi}{\partial y}, \frac{\partial\psi}{\partial x}) & \forall(\textbf{x}, t)\in\Omega \times [0,\infty] \\
c(\textbf{x},t=0)=c_0 & \forall\textbf{x}\in\Omega \\
\theta(\textbf{x},t=0)=\theta_0 & \forall\textbf{x}\in\Omega \\
s(\textbf{x},t=0)=s_0 & \forall\textbf{x}\in\Omega \\
c=c_{\text{D}} & \forall(\textbf{x}, t)\in\partial\Omega_{\text{D}, c} \times [0,\infty] \\
\textbf{n}\cdot(\mathsf{D}\cdot\nabla c) = c_{\text{N}} & \forall(\textbf{x}, t)\in\partial\Omega_{\text{N}, c}
\times [0,\infty]~,~\partial\Omega_{\text{N}, c}=\partial\Omega/\partial\Omega_{\text{D}, c} \\
\theta=\theta_{\text{D}} & \forall (\textbf{x}, t)\in\partial\Omega_{\text{D}, \theta} \times [0,\infty] \\
\textbf{n}\cdot(\mathsf{G}\cdot\nabla \theta) = \theta_{\text{N}} & \forall(\textbf{x}, t)\in\partial\Omega_{\text{N}, \theta}
\times [0,\infty]~,~\partial\Omega_{\text{N}, \theta}=\partial\Omega/\partial\Omega_{\text{D}, \theta} \\
\psi = \psi_{\text{D}} & \forall(\textbf{x}, t)\in\partial\Omega\times [0,\infty]
\end{cases}~.
\end{align*}
$$


## Constitutive relations

The above equations are closed by the constitutive relations
$$
\begin{cases}
\textbf{e}_g[~=\text{$-(\textbf{e}_z$ if $\Omega\sub\mathbb{R}^3$ else $\textbf{e}_y$)~]} & \text{gravitational unit vector} \\
\varphi(\textbf{x})[~=1] & \text{rock porosity}\\
\mathsf{D}(\phi, \textbf{u})[~=\mathsf{I}~] & \text{solutal dispersion} \\
\mathsf{G}(\phi, \textbf{u})[~=\mathsf{I}~]  & \text{thermal dispersion} \\
K(\phi)[~=\phi^2\mathsf{I}~] & \text{permeability}\\
\mu(c, \theta)[~=1~] & \text{fluid viscosity } \\
\rho(c, \theta)[~=c-\theta~] & \text{fluid density}\\
R(s, c, \theta)[~=s(1-c)\theta~] & \text{reaction rate} \\
\end{cases}
$$

which assume their default values denoted inside square brackets if not prescribed. 

## Boundary conditions

The partial differential equations are constrained by the boundary conditions

$$
\begin{cases}
\partial\Omega_{\text{D},\theta}[~=\varnothing~]~,~c\vert_{\partial\Omega_{\text{D},\theta}}\equiv\theta_{\text{D}} & \text{thermal Dirichlet} \\
\partial\Omega_{\text{N},\theta}[~=\partial\Omega~]~,~c\vert_{\partial\Omega_{\text{N},\theta}}\equiv\theta_{\text{N}}=[~0~] & \text{thermal Neumann} \\
\partial\Omega_{\text{D},c}[~=\varnothing~]~,~c\vert_{\partial\Omega_{\text{D},c}}\equiv c_{\text{D}} & \text{solutal Dirichlet} \\
\partial\Omega_{\text{N},c}[~=\partial\Omega~]~,~c\vert_{\partial\Omega_{\text{N},c}}\equiv c_{\text{N}}[~=0~] & \text{solutal Neumann} \\
\partial\Omega_{\text{E}}[~=\partial\Omega~]~,~(\textbf{n}\cdot\textbf{u})\vert_{\partial\Omega_{\text{E}}} \equiv u_{\text{E}}=[~=0~] & \text{essential} \\
\partial\Omega_{\text{N}}[~=\varnothing~]~,~p\vert_{\partial\Omega_{\text{N}}} \equiv p_{\text{N}} & \text{natural} \\
~\psi\vert_{\partial\Omega}\equiv\psi_{\text{D}}[~=0~] & \text{streamfunction Dirichlet} 
\end{cases}
$$

which assume their default values denoted inside square brackets if not prescribed. 