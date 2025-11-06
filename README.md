# LUCiFEx CO<sub>2</sub> Package

## Installation (macOS)

See `https://github.com/george-poole/LUCiFEx` to install the `lucifex` package.

`git clone https://github.com/george-poole/CO2DissolutionPackage.git`

## Governing equations

This package solves a system of PDEs modelling flow in a porous medium coupled to the thermosolutal transport of dissolved CO<sub>2</sub> and porosity evolution due to the dissolution of capillary-trapped CO<sub>2</sub>.

In strong form, the non-dimensionalized initial boundary value problem is
$$
\begin{align*}
&\text{Find} \\
&\text{$c(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$, } \\
&\text{$\theta(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$, } \\
&\text{$s(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$, } \\
&\text{$\textbf{u}(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}^d$ and $p(\textbf{x}, t): \Omega\times[0, \infty) \to \mathbb{R}$} \\
&\text{such that} \\
&\begin{cases}
\phi\frac{\partial c}{\partial t} + Ad\,\textbf{u}\cdot\nabla c = \frac{1}{Pe}\nabla\cdot(\mathsf{D}(\phi, \textbf{u})\cdot\nabla c) + KiR(s, c, \theta) & \\
\phi\frac{\partial\theta}{\partial t} + Ad\,\textbf{u}\cdot\nabla\theta = \frac{1}{LePe}\nabla\cdot(\mathsf{G}(\phi, \textbf{u})\cdot\nabla\theta) & \\
\nabla\cdot\textbf{u} = 0 & \\
\textbf{u}=-\frac{\mathsf{K}(\phi)}{\mu(c, \theta)}\cdot(\nabla p + Bu\,\rho(c, \theta) g\,\textbf{e}_g) \\
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

where $\phi(\textbf{x}, t)=\varphi(\textbf{x})(1 - s(\textbf{x}, t))$ is the effective porosity. Alternatively in 2D with $\partial\Omega_{\text{E}}=\partial\Omega\iff\Omega_{\text{N}}=\varnothing$, the streamfunction formulation solving for $\psi(\textbf{x}, t): \Omega\times[0,\infty]\to\mathbb{R}$ instead of $\textbf{u}$ and $p$ contains

$$
\begin{cases}
\textbf{u} = (-\frac{\partial\psi}{\partial y}, \frac{\partial\psi}{\partial x}) & \\
\nabla\cdot\bigg(\frac{\mu\mathsf{K}^{\mathsf{T}}\cdot\nabla\psi}{\det\mathsf{K}}\bigg) &= -\frac{\partial(\rho\,\textbf{e}_g\cdot\textbf{e}_y)}{\partial x} + \frac{\partial(\rho\,\textbf{e}_g\cdot\textbf{e}_x)}{\partial y} & \forall(\textbf{x}, t)\in\Omega \times [0,\infty] \\
\psi=\psi_{\text{D}} & \forall(\textbf{x}, t)\in\partial\Omega \times [0,\infty] \\
\end{cases}
$$

The following constitutive relations
$$
\begin{cases}
\textbf{e}_g & \text{gravitational unit vector} & \text{$-\textbf{e}_y$ in 2D or $-\textbf{e}_z$ in 3D}  \\
\varphi(\textbf{x}) & \text{inert rock formation's porosity} & 1\\
\mathsf{D}(\phi, \textbf{u}) & \text{solutal dispersion} & \mathsf{I} \\
\mathsf{G}(\phi, \textbf{u}) & \text{thermal dispersion} & \mathsf{I} \\
K(\phi) & \text{permeability} & \phi^2\mathsf{I}\\
\mu(c, \theta) & \text{fluid viscosity} & 1\\
\rho(c, \theta) & \text{fluid density} & c-\theta\\
R(s, c, \theta) & \text{reaction rate} & s(1-c)\theta \\
\end{cases}
$$

are to be prescribed. Default values denoted in the rightmost column model a horizontal, isotropic, homogeneous porous medium containing an isoviscous fluid. The default boundary conditions, defined on the entire boundary, are
$$
\begin{cases}
c_{\text{N}}\vert_{\partial\Omega}=0 & \text{no solutal flux} \\
\theta_{\text{N}}\vert_{\partial\Omega}=0 & \text{no thermal flux} \\
u_{\text{E}}\vert_{\partial\Omega}=0 \iff\psi_{\text{D}}\vert_{\partial\Omega}=0 & \text{no fluid penetration}
\end{cases}
$$

With a length scale $\mathcal{L}$, velocity scale $\mathcal{U}$, time scale $\mathcal{T}$ and other scales denoted by $\Delta\,\cdot$ or $\cdot~_{\text{ref}}$, the dimensionless numbers appearing in the governing equations are
$$
\begin{align*}
Ad&=\frac{\mathcal{U}\mathcal{T}}{\phi_{\text{ref}}\mathcal{L}} \\
Pe&=\frac{\phi_{\text{ref}}\mathcal{L}^2}{D_{\text{ref}}\mathcal{T}} \\
Ki&=\frac{\mathcal{T}\Delta R}{\phi_{\text{ref}}\Delta c} \\
Bu&=\frac{K_{\text{ref}}\,g\Delta\rho}{\mu_{\text{ref}}\,\mathcal{U}} \\
\end{align*}
$$

$Le$ is the Lewis number for the ratio of thermal to solutal diffusivity; $Ra$ is the Rayleigh number (defined with respect to the transport of $c$ and domain length scale) is
$$Ra=\frac{\mathcal{L}_\Omega K_{\text{ref}}g\Delta\rho}{\mu_{\text{ref}}D_{\text{ref}}}=\underbrace{\frac{K_{\text{ref}}\,g\Delta\rho}{\mu_{\text{ref}}}}_{\text{convective speed}} \big/ \underbrace{\frac{D_{\text{ref}}}{\mathcal{L}_\Omega}}_{\text{diffusive speed}}$$

and $Da$ is the Damköhler number (defined with respect to the transport of $c$ and domain length scale) is

$$Da=\frac{\mathcal{L}_\Omega \mu_{\text{ref}}\,\Delta R}{K_{\text{ref}}\,g\Delta\rho\Delta c} = \underbrace{\frac{\Delta R}{\Delta c}}_{\text{reaction rate}} \big/ \underbrace{\frac{K_{\text{ref}}\,g\Delta\rho}{\mathcal{L}_\Omega \mu_{\text{ref}}}}_{\text{convection rate}}$$

A particular choice of non-dimensionalization with $\mathcal{L}$, $\mathcal{U}$ and $\mathcal{T}$ chosen according to the details of the problem (e.g. boundary conditions, dominant transport processes) will map $\{Ad, Pe, Ki, Bu, Xl\}$ onto combinations of the Rayleigh, Damköhler and Lewis numbers.

| $\mathcal{L}$ | $\mathcal{U}$ |$ \mathcal{T}$ | $\{Ad, Pe, Ki, Bu, Xl\}$ | Examples | 
| -------- | ------- | ------- | ------- | ------- |
| $\mathcal{L}_\Omega$  |  $K_{\text{ref}}\,g\Delta\rho/\mu_{\text{ref}}$  | $\phi_{\text{ref}}\mathcal{L}/\mathcal{U}$ | $\{1, Ra, Da, 1, 1\}$| [Hewitt et al. (2012)](https://link.aps.org/doi/10.1103/PhysRevLett.108.224503) |
| $D_{\text{ref}}/\mathcal{U}$  |  $K_{\text{ref}}\,g\Delta\rho/\mu_{\text{ref}}$  | $\phi_{\text{ref}}\mathcal{L}/\mathcal{U}$ | $\{1, 1, Da/Ra, 1, Ra\}$| [Slim (2014)](https://www.cambridge.org/core/product/identifier/S0022112013006733/type/journal_article) | 
| $\mathcal{L}_\Omega$  |  $D_{\text{ref}}/\mathcal{L}$  | $\phi_{\text{ref}}\mathcal{L}/\mathcal{U}$ | $\{1, 1, RaDa, Ra, 1\}$| [Ritchie \& Pritchard  (2011)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/natural-convection-and-the-evolution-of-a-reactive-porous-medium/71E5FB557F61CB9125E5B4E4EE9D828F) | 
| $\sqrt{D_{\text{ref}}\mathcal{T}/\phi_{\text{ref}}}$  |  $\phi_{\text{ref}}\mathcal{L}/\mathcal{T}$  | $\phi_{\text{ref}}\Delta c/\Delta R$  | $\{1, 1, 1, \sqrt{Ra/Da}, \sqrt{RaDa}\}$| [Kabbadj et al. (2025)](https://nlpc.ulb.be/pdf/25.Kabbadj_MATRIX.pdf) |


The above mappings between $\{Ra, Da\}$ and $\{Ad, Pe, Ki, Bu, Xl\}$ are implemented by the enumeration class `co2_pkg.sim.ScalingType` to avoid hard-coding simulations with a fixed non-dimensionlization.


## Discretization

TODO


## Models

User-defined choices of constitutive relations, boundary conditions and initial conditions may be prescribed to investigate novel models. Solutal and thermal transport may be 'switched off' by setting the solutal or thermal dispersion relation to `None`. Setting the density ratio parameter to `None` or $\varepsilon=0$ switches off the porosity evolution. 

$$
\text{Model A~}
\begin{cases}
\mathcal{L}, \mathcal{U}, \mathcal{T}= \\
\Omega=[0, L_x] \times [0, L_y] \\ 
\rho(c)=c \\
r(s,c)s(1-c) \\
s_0(y)=s_r\text{H}(y-h_0) \\
c_0(x,y)=c_r\text{H}(y-h_0) + \mathcal{N}(x,y)
\end{cases}
$$

$$
\text{Model B~}
\begin{cases}
\mathcal{L}, \mathcal{U}, \mathcal{T}= \\
\Omega=[0, L_x] \times [0, L_y] \\ 
\rho(c)=c - \gamma\theta \\
r(s,c, \theta)=s(1+\delta\theta-c) \\
s_0=s_r \\
c_0(x,y)=f(y) + \mathcal{N}(x,y)\\
\theta_0(x,y)=1 - y + \mathcal{N}(x, y) \\
\theta_{\text{D}}(x,y=0)=1 \\
\theta_{\text{D}}(x, y=L_y)=0 \\
\theta_{\text{N}}(x=0,y)=0 \\
\theta_{\text{N}}(x=L_x, y)=0
\end{cases}
$$ 

Classic models such as Rayleigh-Taylor or Rayleigh-Bénard convection are also recoverable in this framework.


$$
\text{Rayleigh-Taylor~}
\begin{cases}
\mathcal{L}, \mathcal{U}, \mathcal{T}= \\
\Omega=[0, L_x] \times [0, 1] \\ 
\rho(\theta)=-\theta \\
\theta_0(x,y)=\text{H}(\tfrac{1}{2} - y) \\
\end{cases}
$$

$$
\text{Rayleigh-Bénard~}
\begin{cases}
\mathcal{L}, \mathcal{U}, \mathcal{T}= \\
\Omega=[0, L_x] \times [0, 1] \\ 
\rho(\theta)=-\theta \\
\theta_0(x,y)=1-y + \mathcal{N}(x,y) \\
\theta_{\text{D}}(x,y=0)=1 \\
\theta_{\text{D}}(x, y=1)=0 \\
\theta_{\text{N}}(x=0,y)=0 \\
\theta_{\text{N}}(x=L_x, y)=0
\end{cases}
$$

## Simulations

### Creating a custom simulation

```python
# file `custom_simulation.py` 
from lucifex.sim import configure_simulation, integrate_from_cli
from co2_dissolution_pkg.sim import abstract_simulation

@configure_simulation(
    write_step=...,
    # optional arguments to configure I/O
)
def custom_simulation(*args, **kwargs):
    # code specifying choices of constitutive relations, 
    # boundary conditions and initial conditions
    return abstract_simulation(*args, **kwargs)

if __name__ == "__main__":
    integrate_from_cli(custom_simulation)
```

Examples of simulations can be found in the module `co2_dissolution_pkg.sim`. 

### Running simulations

#### From the command line

In the terminal, <br>
`python custom_simulation.py --help` to list all argument names and default values.

To run a simulation with the argument `X` set to $X=X_0$ 
<br>
`python custom_custom_simulation.py --X X0`<br>

To run multiple simulations in parallel, the command line utility GNU `parallel` is recommended.

To run `N_PROC` simultaneous processes with the arguments `X` and `Y` taking values $(X, Y)\in\{(X_0, Y_0), (X_0, Y_1), (X_1, Y_0), (X_1, Y_1)\}$
<br>
`parallel -j N_PROC "python custom_simulation.py --X {1} --Y {2}" ::: X0 X1 ::: Y0 Y1` <br>

To run `N_PROC` simultaneous processes with the arguments `X` and `Y` taking values $(X, Y)\in\{(X_0, Y_0), (X_1, Y_1), (X_2, Y_2)\}$
<br>
`parallel -j N_PROC --link "python custom_simulation.py --X {1} --Y {2}" ::: X0 X1 X2 ::: Y0 Y1 Y2` <br>

Further command line utilities:
* `caffeinate` e.g. `caffeinate -d -i -s -t <SECONDS> <COMMAND>` to prevent sleeping
* `nohup` e.g. `nohup <COMMAND> & disown` to run without interruption
* `htop` and `kill` for process mangement
* (Ctrl + Z) followed by `bg` to move process to background <br>

Timeseries data will be written to disk in accordance with the `write_step` argument given at the command line, or if unspecified with its value configured in the decorator function `configure_simulation`. 

Postprocessing can the be carried out in a separate script or notebook by making use of the `lucifex.io` and `lucifex.viz` modules.

#### In an iPython notebook

For example,
```python
# file `custom_simulation.ipynb`
from lucifex.sim import integrate
from custom_simulation import custom_simulation

simulation = custom_simulation(store_step=...)(*args, **kwargs)
integrate(simulation, n_stop=..., t_stop=...)
```

Provided that `store_step` is not `None`, timeseries data is kept in memory by the simulation object once the `integrate` routine has finished. A quantity stored under the name `'q'` can be accessed as `simulation['q']` for postprocessing in the iPython notebook.
