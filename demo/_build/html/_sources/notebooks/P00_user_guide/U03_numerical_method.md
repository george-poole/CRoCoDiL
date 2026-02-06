# Numerical method

## Linearization


Decompose the reaction into a term linear in $c$ and a source term.
$$
\Sigma(s,c,\theta) = R(s,\theta)c + J(s,\theta)
$$

Each governing equation is linear with respect to its own variable, so solving the equations sequentially in a suitable order and with a suitable choice of discretization in time linearizes the otherwise nonlinear advection, dispersion and reaction terms.

## Discretization in time

Finite difference method

$$
\begin{align}
\mathcal{D}_\phi(\phi)\frac{\theta^{n+1} - \theta^n}{\Delta t^n} + Ad\,\mathcal{D}_{\textbf{u},\theta}(\textbf{u}\cdot\nabla\theta) &= \frac{1}{LePe}\nabla\cdot(\mathcal{D}_{\mathsf{G}}(\mathsf{G})\cdot\nabla\mathcal{D}_{\theta}(\theta)) \\
\mathcal{D}_\phi(\phi)\frac{c^{n+1} - c^n}{\Delta t^n} + Ad\,\mathcal{D}_{\textbf{u},c}(\textbf{u}\cdot\nabla c) &= \frac{1}{Pe}\nabla\cdot(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\mathcal{D}_{c}(c)) + Ki\left(\mathcal{D}_{R,c}(Rc) + \mathcal{D}_{J}(J)\right) \\
\nabla\cdot\textbf{u}^n &= 0 \\
\textbf{u}^n &= - \frac{\mathsf{K}^n}{\mu^n}(\nabla p^n + Bu\,\rho^n\textbf{e}_g) \\
\varphi\frac{s^{n+1} - s^n}{\Delta t^n} &= -\varepsilon Ki\,\mathcal{D}_{s,c,\theta}(\Sigma)
\end{align}
$$

## Adaptive timestep

$$
\Delta t^n=\min\{\Delta t_{\text{max}}, \max\{\Delta t_{\text{min}}, \Delta t^n_{\text{CFLR}}\}\}
$$

where

$$
\begin{align*}
\Delta t^n_{\text{CFLR}} &= \min\{C_{\text{CFL}}\Delta t^n_{\text{CFL}}, C_R\Delta t^n_R\} \\
\Delta t^n_{\text{CFL}} &= \min_{\textbf{x}}\bigg(\frac{h(\textbf{x})}{Ad\,|\textbf{u}^n|}\bigg)  \\
\Delta t^n_{\text{R}} &= \frac{1}{\max_{\textbf{x}}(Ki\,|R^n|)}
\end{align*}
$$

Default values are $C_{\text{CFL}}=0.5$ and $C_{\text{R}}=0.1$ for safety.

## Discretization in space

finite element method

mesh $\Omega\approx\cup_{\mathcal{K}}\mathcal{K}$

trial and test function spaces $V, \hat{V}$

multiply by a test function in $\hat{V}$ and integrate over the domain $\Omega$, integrating by parts where necessary to apply boundary conditions.

weak vs strong enforcement

Obtain $\mathsf{A}\cdot\textbf{u}=\textbf{b}$ where the structure of the matrix $\mathsf{A}$ should be used to make a smart choice of the preconditioner matrix $\mathsf{P}$ used in a Krlov subspace method.

### Mixed formulation


$$
V_{\textbf{u}} = \hat{V}_{\textbf{u}} =
$$

$$
V_p = \hat{V}_p =
$$


### Flow equations

#### Streamfunction formulation

continuous Lagrange finite element space of degree $k$
$$
V_\psi = \hat{V}_\psi = \{v\in C^0(\Omega)~:~v|_{\partial\Omega}=\psi_{\text{D}}~,~v|_{\mathcal{K}}\in\text{P}_k(\mathcal{K})~\forall\mathcal{K}\in\cup_{\mathcal{K}}\mathcal{K}\}
$$

$$
\begin{align*}
&\text{Find $\psi^n\in V_\psi$ such that } \\
&F_\psi(\psi^n, v)\equiv-\int_\Omega\text{d}\Omega~\nabla v\cdot\bigg(\frac{\mu^n\,(\mathsf{K}^n)^{\mathsf{T}}\cdot\nabla\psi^n}{\det(\mathsf{K}^n)}\bigg)+v\frac{\partial(\rho^n\textbf{e}_g\cdot\textbf{e}_y)}{\partial x} - v\frac{\partial(\rho^n\textbf{e}_g\cdot\textbf{e}_x)}{\partial y}=0 ~~~\forall v\in\hat{V}_\psi
\end{align*}
$$

Matrix is symmetric, so algebraic multigrid preconditioning is effective.

### Transport equations

#### Continuous Galerkin formulation

continuous Lagrange finite element space of degree $k$
$$
V_c = \hat{V}_c = \{v\in C^0(\Omega)~:~v|_{\partial{\Omega}_{\text{D}, c}}=c_{\text{D}}~,~v|_{\mathcal{K}}\in\text{P}_k(\mathcal{K})~\forall\mathcal{K}\in\cup_{\mathcal{K}}\mathcal{K}\}
$$
$$
V_\theta = \hat{V}_\theta  = \{v\in C^0(\Omega)~:~v|_{\partial{\Omega}_{\text{D}, \theta }}=\theta_{\text{D}}~,~v|_{\mathcal{K}}\in\text{P}_k(\mathcal{K})~\forall\mathcal{K}\in\cup_{\mathcal{K}}\mathcal{K}\}
$$

$$
\begin{align*}
\text{Find $c^{n+1}\in V_c$ such that } \\
F_c(c^{n+1}, v)&\equiv
\int_\Omega\text{d}\Omega~v\frac{c^{n+1}-c^n}{\Delta t^n} \\
&\quad +Ad \int_\Omega\text{d}\Omega~v\frac{\mathcal{D}_{\textbf{u}, c}(\textbf{u}\cdot\nabla c)}{\mathcal{D}_\phi(\phi)} \\
&\quad +\frac{1}{Pe}\int_\Omega\text{d}\Omega~\nabla\bigg(\frac{v}{\mathcal{D}_\phi(\phi)}\bigg)\cdot(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla \mathcal{D}_c(c)) \\
&\quad -Ki\int_\Omega\text{d}\Omega~v\frac{\mathcal{D}_{R,c}(Rc) \\
+ \mathcal{D}_{J}(J)}{\mathcal{D}_\phi(\phi)} \\
&\quad -\frac{1}{Pe}\int_{\partial\Omega_{\text{N},c}}\text{d}\Gamma~\frac{vc_{\text{N}}}{\mathcal{D}_\phi(\phi)} \\
&\quad + F^{\text{SUPG}}_c(c^{n+1}, v) \\
& =0 \qquad \forall v\in\hat{V}_c
\end{align*}
$$

$$
\begin{align*}
\text{Find $c^{n+1}\in V_\theta$ such that } \\
F_\theta(\theta^{n+1}, v)&\equiv
\int_\Omega\text{d}\Omega~v\frac{\theta^{n+1}-\theta^n}{\Delta t^n} \\
&\quad +Ad\int_\Omega\text{d}\Omega~v\frac{\mathcal{D}_{\textbf{u}, \theta}(\textbf{u}\cdot\nabla \theta)}{\mathcal{D}_\phi(\phi)} \\
&\quad +\frac{1}{LePe}\int_\Omega\text{d}\Omega~\nabla\bigg(\frac{v}{\mathcal{D}_\phi(\phi)}\bigg)\cdot(\mathcal{D}_{\mathsf{G}}(\mathsf{G})\cdot\nabla \mathcal{D}_\theta(\theta)) \\
&\quad -\frac{1}{LePe}\int_{\partial\Omega_{\text{N}, \theta}}\text{d}\Gamma~\frac{v\theta_{\text{N}}}{\mathcal{D}_\phi(\phi)} \\
&\quad +F^{\text{SUPG}}_\theta(\theta^{n+1}, v) \\
&=0 \qquad \forall v\in\hat{V}_\theta
\end{align*}
$$

unstable if advection-dominated flow, add SUPG stabilization term

$$F^{\text{SUPG}}_c(c^{n+1}, v) = \tau^n_c(\nabla v\cdot \textbf{u}^{n}_{\text{eff}, c})\mathcal{R}_c^{n+1}$$
$$F^{\text{SUPG}}_\theta(\theta^{n+1}, v) = \tau^n_\theta(\nabla v\cdot \textbf{u}^{n}_{\text{eff},\theta})\mathcal{R}_\theta^{n+1}$$

residual

$$
\mathcal{R}_c^{n+1}=\frac{c^{n+1}-c^n}{\Delta t^n} + \frac{\mathcal{D}_{\textbf{u},c}(\textbf{u}\cdot\nabla c)}{\mathcal{D}_\phi(\phi)} - \frac{1}{Pe}\frac{\nabla\cdot(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla \mathcal{D}_c(c))}{\mathcal{D}_\phi(\phi)} - Ki\frac{\mathcal{D}_{R,c}(Rc) + \mathcal{D}_{J}(J)}{\mathcal{D}_\phi(\phi)}
$$

$$
\mathcal{R}_\theta^{n+1}=
$$

stabilization parameter

$$
\tau^n=\tau^n(h, \textbf{u}_{\text{eff}}, \mathsf{D}_{\text{eff}}, R_{\text{eff}})=\begin{cases}
(2|\textbf{u}_{\text{eff}}| / h  +  4D_{\text{eff}} / h^2  -  R_{\text{eff}})^{-1} & \text{Codina} \\
((2|\textbf{u}_{\text{eff}}| / h)^2  +  9(4D_{\text{eff}} / h^2)^2  +  R_{\text{eff}}^2)^{-1/2} & \text{Shakib} \\
((2 / \Delta t)^2 + (2|\textbf{u}_{\text{eff}}| / h)^2  +  (2D_{\text{eff}} / h^2)^2  +  R_{\text{eff}}^2)^{-1/2} & \text{transient} \\
(h/2|\textbf{u}_{\text{eff}}|)(\text{coth}(|\textbf{u}_{\text{eff}}|h/2D_{\text{eff}}) - 2D_{\text{eff}}/|\textbf{u}_{\text{eff}}|h) & \text{coth} \\
0 & \text{none}
\end{cases}
$$


effective quantities with reference to steady state problem such that with $\mathsf{D}\approx\text{tr}(\mathsf{D})/\text{dim}(\mathsf{D})\mathsf{I}$

$$
\begin{align*}
&\frac{c^{n+1} - c^n}{\Delta t^n} + Ad\frac{\mathcal{D}_{\textbf{u},c}(\textbf{u}\cdot\nabla c)}{\mathcal{D}_\phi(\phi)} - \frac{1}{Pe}\frac{\nabla\cdot(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\mathcal{D}_{c}(c))}{\mathcal{D}_\phi(\phi)} - Ki\frac{\mathcal{D}_{R,c}(Rc) + \mathcal{D}_{J}(J)}{\mathcal{D}_\phi(\phi)} \\
&=\textbf{u}_{\text{eff},c}\cdot\nabla c^{n+1} - \nabla\cdot(D_{\text{eff}}\nabla c^{n+1}) - R_{\text{eff}}c^{n+1} - J_{\text{eff}}
\end{align*}
$$

Due to the advection term, matrix is asymmetric. generalized minimal residual method (GMRES) effective. 

#### Disontinuous Galerkin formulation

In contrast to the continuous formulation, Dirichlet boundary conditions are not strongly enforced in the definition of the discontinuous function space, but are instead weakly enforced by terms in the variational formulation with penalty parameters $\alpha_{\text{DG}}$ and $\gamma_{\text{DG}}$.

$$
V_c = \hat{V}_c = \{v\in L^2(\Omega)~:~v|_{\partial{\Omega}_{\text{D}, c}}~v|_{\mathcal{K}}\in\text{P}_k(\mathcal{K})~\forall\mathcal{K}\in\cup_{\mathcal{K}}\mathcal{K}\}
$$
$$
V_\theta = \hat{V}_\theta = \{v\in L^2(\Omega)~:~v|_{\partial{\Omega}_{\text{D},\theta}}~v|_{\mathcal{K}}\in\text{P}_k(\mathcal{K})~\forall\mathcal{K}\in\cup_{\mathcal{K}}\mathcal{K}\}
$$

jump and average operators

$$
\begin{align*}
\{v\} &= \\
\{\textbf{v}\} &= \\
\llbracket v\rrbracket &= \\
\llbracket \textbf{v}\rrbracket &= 
\end{align*}
$$

$$
\begin{align*}
\text{Find $c^{n+1}\in V_c$ such that } \\
F_c(c^{n+1}, v)&\equiv
\int_\Omega\text{d}\Omega~v\frac{c^{n+1}-c^n}{\Delta t^n} \\
&\quad -Ad \int_\Omega\text{d}\Omega~\,\nabla\left(\frac{v}{\mathcal{D}_\phi(\phi)}\right)\cdot\nabla\mathcal{D}_{\textbf{u},c}(\textbf{u}c) \\
&\quad + Ad\int_{\mathcal{F}/\partial\Omega}\text{d}\Gamma~2\llbracket v\rrbracket\textbf{n}\cdot\begin{cases}
\{\mathcal{D}_{\textbf{u},c}(\textbf{u}c)\} & \textbf{n}\cdot\textbf{u} > 0 \\
0 & \text{otherwise} \\
\end{cases} \\
&\quad + Ad\int_{\partial\Omega}\text{d}\Gamma~\,\frac{v}{\mathcal{D}_\phi(\phi)}\begin{cases}
0 & \textbf{n}\cdot\textbf{u} >0 \\
c_{\text{N}}\,\textbf{n}\cdot\textbf{u} & \text{otherwise} \\
\end{cases} \\
&\quad +\frac{1}{Pe}\int_\Omega\text{d}\Omega~\nabla\left(\frac{v}{\mathcal{D}_\phi(\phi)}\right)\cdot\nabla(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\mathcal{D}_c(c)) \\
&\quad - \frac{1}{Pe}\int_{\mathcal{F}/\partial\Omega}\text{d}\Gamma~\left\{\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\left(\frac{v}{\mathcal{D}_\phi(\phi)}\right)\right\}\cdot\llbracket\mathcal{D}_c(c)\rrbracket\textbf{n}+\left\{\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\mathcal{D}_c(c)\right\}\cdot\llbracket \frac{v}{\mathcal{D}_\phi(\phi)}\rrbracket\textbf{n} \\
&\quad +\frac{1}{Pe}\int_{\mathcal{F}/\partial\Omega}\text{d}\Gamma~ \frac{\alpha_{\text{DG}}}{h}\llbracket \frac{v}{\mathcal{D}_\phi(\phi)}\rrbracket\cdot\llbracket\mathcal{D}_c(c)\rrbracket\\
&\quad -\frac{1}{Pe}\quad \int_{\partial\Omega_{\text{D},c}}\text{d}\Gamma~\left(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\left(\frac{v}{\mathcal{D}_\phi(\phi)}\right)\right)\cdot(\mathcal{D}_c(c)-c_{\text{D}})\textbf{n}\\

&\quad -\frac{1}{Pe}\quad \int_{\partial\Omega_{\text{D},c}}\text{d}\Gamma~(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\mathcal{D}_c(c))\cdot\left(\frac{v}{\mathcal{D}_\phi(\phi)}\textbf{n}\right)\\

&\quad +\frac{1}{Pe}\quad \int_{\partial\Omega_{\text{D},c}}\text{d}\Gamma~\frac{\gamma_{\text{DG}}}{h} v(c-c_{\text{D}})\\
&\quad -\frac{1}{Pe}\quad \int_{\partial\Omega_{\text{N},c}}\text{d}\Gamma~\frac{vc_{\text{N}}}{\mathcal{D}_\phi(\phi)}\\
&\quad -Ki\int_\Omega\text{d}\Omega~v\,\frac{\mathcal{D}_{R,c}(Rc)}{\mathcal{D}_\phi(\phi)} \\
&=0 \qquad \forall v\in\hat{V}_c 
\end{align*}
$$


## Algorithm

1. Initialize namespace $\mathcal{Q}=\{t, \Delta t, s\}$
1. Initialize $n\gets0$
1. Initialize $t^n\gets0$
1. **if** $\exists\psi$ **then**
1. &emsp; Update namespace $\mathcal{Q} \gets \mathcal{Q}\cup\{\psi, \textbf{u}\}$
1. **else**
1. &emsp; Update namespace $\mathcal{Q} \gets \mathcal{Q}\cup\{\textbf{u}, p\}$
1. **if** $\exists\theta$ **then**
1. &emsp; Update namespace $\mathcal{Q} \gets \mathcal{Q}\cup\{\theta\}$
1. &emsp; Initialize $\theta^n \gets \theta_0$
1. **if** $\exists c$ **then**
1. &emsp; Update namespace $\mathcal{Q} \gets \mathcal{Q}\cup\{c\}$
1. &emsp; Initialize $c^n \gets c_0$
1. **while** $t^n<t_{\text{stop}}$ **and** $n<n_{\text{stop}}$ **do**
1. &emsp; **if** $\exists\psi$ **then**
1. &emsp;&emsp; Update $\psi^n$ by solving ...
1. &emsp;&emsp; Update $\textbf{u}^n$ by projecting or interpolating ...
1. &emsp; **else**
1. &emsp;&emsp; Update $\textbf{u}^n, p^n$ by solving ...
1. &emsp; **if** $n<n_{\text{init}}$ **then**
1. &emsp;&emsp; Update $\Delta t^n\gets\Delta t_{\text{init}}$
1. &emsp; **else**
1. &emsp;&emsp; Update $\Delta t^n\gets\Delta t^n_{\text{CFLR}}$
1. &emsp; Update $t^{n+1}\gets t^n + \Delta t^n$
1. &emsp; **if** $\exists\theta$ **then**
1. &emsp;&emsp; Update $\theta^{n+1}$ by solving ... 
1. &emsp; **if** $\exists c$ **then**
1. &emsp;&emsp; Update $c^{n+1}$ by solving ...
1. &emsp; **if** $\varepsilon\neq 0$ **then**
1. &emsp;&emsp; Update $s^{n+1}$ by solving ...
1. &emsp; Update $q^{n-j}\gets q^{n-j+1}~\forall q\in\mathcal{Q}$ and $\forall j\geq 0$ 
1. &emsp; Update $n\gets n+1$
1. **end while**
