# Numerical method

## Linearization


Decompose the reaction into a term linear in $c$ and a source term.
$$
\Sigma(s,c,\theta) = R(s,\theta)c + J(s,\theta)
$$

Each governing equation is linear with respect to its own variable, so solving the equations sequentially in suitable order and with a suitable choice of discretization in time linearizes the otherwise nonlinear advection, dispersion and reaction terms.

## Discretization in time

Finite difference method

$$
\begin{align}
\mathcal{D}_\phi(\phi)\frac{\theta^{n+1} - \theta^n}{\Delta t^n} + Ad\,\mathcal{D}_{A,\theta}(\textbf{u}\cdot\nabla\theta) &= \frac{1}{LePe}\nabla\cdot(\mathcal{D}_{\phi}(\mathsf{G})\cdot\nabla\mathcal{D}_{D,\theta}(\theta)) \\
\mathcal{D}_\phi(\phi)\frac{c^{n+1} - c^n}{\Delta t^n} + Ad\,\mathcal{D}_{A,c}(\textbf{u}\cdot\nabla c) &= \frac{1}{Pe}\nabla\cdot(\mathcal{D}_{\phi}(\mathsf{D})\cdot\nabla\mathcal{D}_{D,c}(c)) + Ki\,\left(\mathcal{D}_{R,c}(Rc) + \mathcal{D}_{J,c}(J)\right) \\
\nabla\cdot\textbf{u}^n &= 0 \\
\textbf{u}^n &= - \frac{\mathsf{K}^n}{\mu^n}(\nabla p^n + Bu\rho^n\textbf{e}_g) \\
\varphi\frac{s^{n+1} - s^n}{\Delta t^n} &= -\varepsilon Da\left(\mathcal{D}_{R,s}(Rc) + \mathcal{D}_{J,s}(J)\right)
\end{align}
$$

| $\mathcal{D}$ operator | Default | Description |
| -------- | ------- | ------- | 
| $\mathcal{D}_\phi$ | `AB(1)` | forward Euler | 
| $\mathcal{D}_{A,\theta}$ | `AB(2) @ CN` | second-order Adams-Bashforth on $\textbf{u}$ and Crank-Nicolson on $\theta$ | 
| $\mathcal{D}_{A,c}$ | `AB(2) @ CN` | second-order Adams-Bashforth on $\textbf{u}$ and Crank-Nicolson on $c$ | 
| $\mathcal{D}_{D,\theta}$ | `CN` | Crank-Nicolson | 
| $\mathcal{D}_{D,c}$ | `CN` | Crank-Nicolson |  
| $\mathcal{D}_{R,c}$ | `AM(2) @ AM(1)` | second-order Adams-Bashforth on $s$ and backward Euler on $c$|  
| $\mathcal{D}_{R,s}$ | `AM(2) @ AM(1)` | second-order Adams-Bashforth on $s$ and backward Euler on $c$|  
| $\mathcal{D}_{J,c}$ | `AM(2) @ AM(1)` | second-order Adams-Bashforth on $s$ and backward Euler on $c$|  
| $\mathcal{D}_{J,s}$ | `AM(2) @ AM(1)` | second-order Adams-Bashforth on $s$ and backward Euler on $c$|  

Ideally set $\mathcal{D}_{R,c}=\mathcal{D}_{R,s}$ to assist mass conservation.

## Adaptive timestep

$$
\Delta t^n=\min\{\Delta t_{\text{max}}, \max\{\Delta t_{\text{min}}, \Delta t^n_{\text{CFLR}}\}\}
$$

where

$$
\begin{align*}
\Delta t^n_{\text{CFLR}} &= \min\{C_{\text{CFL}}\Delta t^n_{\text{CFL}}, C_R\Delta t^n_R\} \\
\Delta t^n_{\text{CFL}} &= \min_{\textbf{x}}\bigg(\frac{\ell(\textbf{x})}{Ad\,|\textbf{u}^n|}\bigg)  \\
\Delta t^n_{\text{R}} &= \frac{1}{\max_{\textbf{x}}(Ki\,|R^n|)}
\end{align*}
$$

Default values are $C_{\text{CFL}}=0.5$ and $C_{\text{R}}=0.1$ for safety.

## Discretization in space

Finite element method

Derive weak form of each equation and solve sequentially.
