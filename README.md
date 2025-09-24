# CO<sub>2</sub> Dissolution Package

## Installation (macOS)

See `https://github.com/george-poole/LUCiFEx` to install the `lucifex` package.

`git clone https://github.com/george-poole/CO2Dissolution.git`

## Governing equations

This package solves a non-dimensionalized system of PDEs describing flow in a porous medium coupled to the advection-diffusion-reaction of dissolved CO<sub>2</sub>, advection-diffusion of heat and and dissolution of capillary-trapped CO<sub>2</sub>. For $(\textbf{x}, t)\in\Omega\times[0, \infty)$, the governing equations are

$$
\begin{align}
\phi\frac{\partial\theta}{\partial t} + \textbf{u}\cdot\nabla\theta &= \frac{1}{Rb}\nabla\cdot(\mathsf{G}\cdot\nabla\theta)~~~,\\
\phi\frac{\partial c}{\partial t} + \textbf{u}\cdot\nabla c &= \frac{1}{Ra}\nabla\cdot(\mathsf{D}\cdot\nabla c) + Da\,R~~~, \\
\nabla\cdot\textbf{u} &= 0 \\
\textbf{u} &= -\frac{\mathsf{K}}{\mu}\cdot(\nabla p + \rho\,\textbf{e}_g) \\
\varphi\frac{\partial s}{\partial t}&=-\varepsilon Da\,R
\end{align}
$$

with constitutive relations for the inert rock formation's porosity $\varphi(\textbf{x})$, permeability $K(\phi)$, solutal dispersion $\mathsf{D}(\phi, \textbf{u})$, thermal dispersion $\mathsf{G}(\phi, \textbf{u})$, fluid density $\mu(c, \theta)$, fluid viscosity $\rho(c, \theta)$ and reaction rate $R(s, c, \theta)$. The porous medium's effective porosity is $\phi = \varphi(1 - s)$ for a saturation $s$ of capillary-trapped CO<sub>2</sub>. 

On $\partial\Omega$ a combination of Dirichlet $(\text{D})$ and Neumann $(\text{N})$ boundary conditions apply for $\theta$ and $c$, or essential $(\text{E})$ and natural $(\text{N})$ boundary conditions for $\textbf{u}$ and $p$.

$$
\begin{align}
\theta &= \theta_{\text{D}}\quad\text{for~}\textbf{x}\in\partial\Omega_{\text{D},\theta} \\
\textbf{n}\cdot(\mathsf{G}\cdot\nabla\theta) &= \theta_{\text{N}}\quad\text{for~}\textbf{x}\in\partial\Omega_{\text{N},\theta}=\partial\Omega/\partial\Omega_{\text{D},\theta} \\
c &= c_{\text{D}}\quad\text{for~}\textbf{x}\in\partial\Omega_{\text{D},c}\\
\textbf{n}\cdot(\mathsf{D}\cdot\nabla c) &= c_{\text{N}}\quad\text{for~}\textbf{x}\in\partial\Omega_{\text{N},c}=\partial\Omega/\partial\Omega_{\text{D},c} \\
\textbf{n}\cdot\textbf{u} &= u_{\text{E}}\quad\text{for~}\textbf{x}\in\partial\Omega_{\text{E}} \\
p &= p_{\text{N}}\quad\text{for~}\textbf{x}\in\partial\Omega_{\text{N}}=\partial\Omega/\partial\Omega_{\text{E}}
\end{align}
$$

Initial conditions are 

$$
\begin{align}
\theta(\textbf{x}, t=0) &= \theta_0(\textbf{x}) \\
c(\textbf{x}, t=0) &= c_0(\textbf{x}) \\
s(\textbf{x}, t=0) &= s_0(\textbf{x})
\end{align}
$$

The streamfunction formulation in 2D

$$
\begin{align}
\textbf{u} &= 
\begin{pmatrix}
-\frac{\partial\psi}{\partial y} \\
\frac{\partial\psi}{\partial x}
\end{pmatrix} \\
\nabla\cdot\bigg(\frac{\mu\mathsf{K}^{\mathsf{T}}\cdot\nabla\psi}{\det\mathsf{K}}\bigg) &= \cos\beta\frac{\partial\rho}{\partial x}-\sin\beta\frac{\partial\rho}{\partial y}
\end{align}
$$

is an alternative to equations $(3)$ and $(4)$, given that the unit vector pointing in the direction of gravity is $\textbf{e}_g=-\sin\beta\textbf{e}_x -\cos\beta\textbf{e}_y$ for a domain inclined at an angle $\beta$ to the horizontal.

Unless otherwise specified, constitutive relations, boundary conditions and initial conditions assume their defaults

* gravity unit vector $\textbf{e}_g=-\textbf{e}_y$ in 2D or $\textbf{e}_g=-\textbf{e}_z$ in 3D
* rock porosity $\varphi=1$
* permeability $K(\phi)=\phi^2$
* solutal dispersion $D(\phi)=\phi$
* solutal dispersion $G(\phi)=\phi$
* fluid density $\rho=c-\theta$
* fluid viscosity $\mu=1$
* no-flux temperature boundary condition $\textbf{n}\cdot(\mathsf{G}\cdot\nabla\theta)=0$ for $\textbf{x}\in\partial\Omega$
* no-flux concentration boundary condition $\textbf{n}\cdot(\mathsf{D}\cdot\nabla c)=0$ for $\textbf{x}\in\partial\Omega$
* no-flux velocity boundary condition $\textbf{n}\cdot\textbf{u}=0$ for $\textbf{x}\in\partial\Omega$

to model an impermeable horizontal domain containing an isotropic, homogeneous porous medium and isoviscous fluid with a linear density.

Solutal or thermal convection can be 'turned off' by setting the solutal Rayleigh number $Ra$ or thermal Rayleigh number $Rb$ to zero. Likewise the reaction and/or dissolution can be 'turned off' by setting the Damkohler number $Da$ and/or density ratio $\varepsilon$ to zero.

Specific choices of constitutive relations, boundary conditions and initial conditions may be user-defined to create novel models. Classic models such as Rayleigh-Taylor or Rayleigh-Benard convection can also be recovered from appropriate choices.

See `co2_dissolution_pkg.math` module for full details of the finite element formulations.

### Novel Models

#### Solutal Convection-Reaction

* rectangular domain $\Omega=[0, L_x] \times [0, L_y]$
* fluid density $\rho(c)=c$ 
* reaction rate $r(s,c)=s(1-c)$
* initial saturation $s(x,y,t=0)=s_r\text{H}(y-h_0)$
* initial concentration with noise $c(x,y,t=0)=c_r\text{H}(y-h_0) + \mathcal{N}(x,y)$
* parameterised by $(Ra, Da, \varepsilon, h_0, s_r, c_r)$

#### Thermosolutal Convection-Reaction

* rectangular domain $\Omega=[0, L_x] \times [0, L_y]$
* fluid density $\rho(c, \theta)=c - \gamma\theta$ 
* reaction rate $r(s,c, \theta)=s(1+\delta\theta-c)$
* temperature boundary conditions $$\begin{align*}
\textbf{n}\cdot\nabla\theta&=0\quad\text{on~} x=0,L_x \\
\theta&=1\quad\text{on~} y=0 \\
\theta&=0\quad\text{on~} y=L_y
\end{align*}$$
* initial saturation $s(x,y,t=0)=s_r$
* initial concentration $c(x,y,t=0)=c_0(y)$
* initial temperature with noise $\theta(x,y,t=0)=1 - y + \mathcal{N}(x, y)$
* parameterised by $(Ra, Rb, Da, \varepsilon, \gamma, \delta, s_r)$

### Classic Models

#### Rayleigh-Taylor Convection

* rectangular domain $\Omega=[0, L_x] \times [0, L_y]$
* fluid density $\rho(\theta)=-\theta$ 
* initial temperature $\theta(x,y,t=0)=\text{H}(h_0 - y)$
* parameterised by $(Rb, h_0)$

#### Rayleigh-Benard Convection

* rectangular domain $\Omega=[0, L_x] \times [0, L_y]$
* fluid density $\rho(\theta)=-\theta$ 
* temperature boundary conditions $$\begin{align*}
\textbf{n}\cdot\nabla\theta&=0\quad\text{on~} x=0,L_x \\
\theta&=1\quad\text{on~} y=0 \\
\theta&=0\quad\text{on~} y=L_y
\end{align*}$$
* initial temperature with noise $\theta(x,y,t=0)=1 - y + \mathcal{N}(x, y)$
* parameterised by $Rb$

## Simulations

### Creating a custom simulation

To define a simulation for a model within the framework of governing equations $(1)$-$(5)$,

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
