# Convective CO<sub>2</sub> Dissolution Package

### Installation (macOS)

See `https://github.com/george-poole/LUCiFEx` to install the `lucifex` package.

`git clone https://github.com/george-poole/CO2Dissolution.git`

### Running Simulations from the Command Line

In the terminal, <br>
`python simulate.py --help` to list all argument names and default values

To run a single simulation <br>
e.g. `python python simulate .py --write_step 1 --n_stop 10 --Ra 420 `<br>
for a simulation with $Ra=420$ integrating 10 steps and writing data at every step

To run multiple simulations in parallel, use GNU `parallel`.

e.g. `parallel -j 2 "python simulate.py --Ra {1} --Da {2}" ::: 200.0 400.0 ::: 5.0 10.0` <br>
for 2 simultaneous simulations with $(Ra, Da)\in\{(200, 5), (200, 10), (400, 5), (400, 10)\}$

e.g. `parallel -j 3 --link "python simulate.py --Ra {1} --Nx {2}" ::: 250 500 750 ::: 100, 150, 200` <br>
for 3 simultaneous simulations with $(Ra, N_x)\in\{(250, 100), (500, 150), (750, 200)\}$

See also `nohup`, `htop`, `bg`, `disown`, `jobs`, `kill` and `caffeinate` (MacOS) for process management.

`<COMMAND>` (Ctrl + Z) then `bg` to move process to background <br>
`<COMMAND> &` to begin process in background <br>

`nohup <COMMAND> & disown`

`caffeinate -d -i -s -t <SECONDS> <COMMAND>`

### Running Simulations in iPython notebooks

In a `.ipynb` notebook, for example
```python
from lucifex.sim import integrate
from simulate import carbon_dissolution

simulation = carbon_dissolution(store_step=1)(Ra=420)
integrate(simulation, *args, **kwargs)
```

Timeseries data is available for postprocessing once `integrate` has finished running.
```python
from lucifex.viz import plot_line
umax = simulation['umax']
fig, ax = plot_line((umax.time_series, umax.value_series))
```

### Exporting Results

To copy a simulation directory containing data and figures <br>
`rsync -r <SOURCE> <DESTINATION>`

or to include the cumbersome data files <br>
`rsync -r --exclude "*.h5" --exclude "*.xdmf" <SOURCE> <DESTINATION>`


### Postprocessing Interactively

`...`

### Postprocessing Non-Interactively

`...`

## Governing Equations


See `formulae.py` for full details of the finite element (in space) and finite difference (in time) discretizations. 

This package solves the non-dimensionalized system of PDEs describing flow in a porous medium

$$\begin{equation}
\tag{$1$}
\begin{matrix}
\nabla\cdot\textbf{u} = 0 \\
\textbf{u} = -\frac{\mathsf{K}}{\mu}\cdot(\nabla p + \rho\,\textbf{e}_g)
\end{matrix}\iff\nabla\cdot\bigg(\frac{\mu\mathsf{K}^{\mathsf{T}}\cdot\nabla\psi}{\det\mathsf{K}}\bigg) = \cos\beta\frac{\partial\rho}{\partial x}-\sin\beta\frac{\partial\rho}{\partial y}~~~,
\end{equation}$$

coupled to the advection-diffusion-reaction of dissolved CO<sub>2</sub> coupled to 
$$\begin{equation}
\tag{$2$}
\phi\frac{\partial c}{\partial t} + \textbf{u}\cdot\nabla c = \frac{1}{Ra}\nabla\cdot(\mathsf{D}\cdot\nabla c) + Da\,R~~~,
\end{equation}
$$

advection-diffusion of heat
$$\begin{equation}
\tag{$2$}
\phi\frac{\partial\theta}{\partial t} + \textbf{u}\cdot\nabla\theta = \frac{1}{Rb}\nabla\cdot(\mathsf{G}\cdot\nabla\theta)~~~,
\end{equation}
$$

and dissolution of capillary-trapped CO<sub>2</sub>
$$\begin{equation}
\tag{$3$}
\varphi\frac{\partial s}{\partial t}=-\varepsilon Da\,R
\end{equation}$$

in an abitrary domain $\Omega$ with prescribed relations for the inert rock porosity $\varphi(\textbf{x})$, permeability $K(\phi)$, solutal dispersion $\mathsf{D}(\phi, \textbf{u})$, thermal dispersion $\mathsf{G}(\phi, \textbf{u})$, fluid density $\mu(c, \theta)$, fluid viscosity $\rho(c, \theta)$ and reaction rate $R(s, c, \theta)$, and furthermore subject to arbitrary initial conditions for $s, c, \theta$ and boundary conditions on $\partial\Omega$ for $c, \theta$ and either $\textbf{u}, p$ or $\psi$. The effective porosity in the medium is $\phi = \varphi(1 - s)$ for a saturation $s$ of capillary-trapped CO<sub>2</sub>. The unit vector pointing in the direction of gravity is $\textbf{e}_g=(\sin\beta, -\cos\beta)$ for a domain inclined at an angle $\beta$ to the horizontal.

Default settings are for an isotropic, homogeneous, horizontal porous medium with isoviscous fluid such that

* rock porosity $\varphi=1$
* permeability $K(\phi)=\phi^2$
* solutal dispersion $D(\phi)=\phi$
* solutal dispersion $G(\phi)=\phi$
* fluid viscosity $\mu=1$ 
* inclination angle $\beta=0$

Solutal or thermal convection can be 'turned off' by setting the solutal Rayleigh number $Ra$ or thermal Rayleigh number $Rb$ to zero. Likewise the reaction and/or dissolution can be 'turned off' by setting the Damkohler number $Da$ and/or density ratio $\varepsilon$ to zero.

### 2D Solutal Case

* rectangular domain $\Omega=[0, L_x] \times [0, L_y]$
* fluid density $\rho(c)=c$ 
* reaction rate $r(s,c)=s(1-c)$
* no-flux boundary condition for the flow $\textbf{n}\cdot\textbf{u}=0\Leftrightarrow\psi=0$ on $\partial\Omega$ <br>
* no-flux boundary condition for the solute $\textbf{n}\cdot\nabla c=0$ on $\partial\Omega$
* initial saturation $s(x,y,t=0)=s_r\text{H}(y-h_0)$
* initial concentration with noise $c(x,y,t=0)=c_r\text{H}(y-h_0) + \mathcal{N}(x,y)$
* physical parameter space $(Ra, Da, \varepsilon, h_0, s_r, c_r)$

### 2D Thermosolutal Case

* rectangular domain $\Omega=[0, L_x] \times [0, L_y]$
* fluid density $\rho(c, \theta)=c - \gamma\theta$ 
* reaction rate $r(s,c, \theta)=s(1+\delta\theta-c)$
* no-flux boundary condition for the flow $\textbf{n}\cdot\textbf{u}=0\Leftrightarrow\psi=0$ on $\partial\Omega$ <br>
* no-flux boundary condition for the solute $\textbf{n}\cdot\nabla c=0$ on $\partial\Omega$
* lateral no-flux boundary conditions for the heat $\frac{\partial\theta}{\partial y}=0$ on $x=0$ and $x=L_x$
* fixed temperature boundary conditions for the heat $\theta=1$ on $y=0$ and $\theta=0$ on $y=L_y$
* initial saturation $s(x,y,t=0)=s_r$
* initial concentration $c(x,y,t=0)=c_0(y)$
* initial temperature with noise $\theta(x,y,t=0)=1 - y + \mathcal{N}(x, y)$
* physical parameter space $(Ra, Rb, Da, \varepsilon, \gamma, \delta, s_r)$