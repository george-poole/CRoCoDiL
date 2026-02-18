# Introduction

## Capillary-trapped CO<sub>2</sub>

![Figure](./figures/pore_scale.png)

$$
\phi + \phi_{\text{CO}_2} + \phi_{\text{rock}} = 1 
$$

$$
\phi(\textbf{x}, t)=\varphi(\textbf{x})(1-s(\textbf{x},t)) 
$$

$$
\begin{align*}
\frac{\partial\varphi}{\partial t}=0&\implies\frac{\partial\phi}{\partial t}=-\varphi\frac{\partial s}{\partial t} \\
s=0 &\iff \phi = \varphi \\
s=1 &\iff \phi = 0 \\
\end{align*}
$$

## Mathematical notation

### Flow & transport

| Symbol | Description |
| -------- | ------- |
| $\textbf{e}_x, \textbf{e}_y, \textbf{e}_z$ | Cartesian unit vectors | 
| $\textbf{x}=(x, y, z) = x\textbf{e}_x + y\textbf{e}_y + z\textbf{e}_z$ | spatial coordinates |
| $t$ | time |
| $s$ | saturation of capillary-trapped CO<sub>2</sub> | 
| $c$ | concentration of dissolved CO<sub>2</sub> | 
| $\theta$ | temperature | 
| $\textbf{u}=(u_x, u_y, u_z)$ | fluid velocity |
| $p$| pressure |
| $\psi$| streamfunction |
| $\rho$ | fluid density |
| $\mu$ | fluid viscosity |
| $g$ | gravity constant |
| $\,{\textbf{e}}_g$ | gravity unit vector |
| $\varphi$ | rock porosity |
| $\phi$ | effective porosity |
| $\mathsf{K}$ | permeability |
| $\mathsf{D}$ | solutal dispersion |
| $\mathsf{G}$ | thermal dispersion |
| $\varepsilon$ | ratio of CO<sub>2</sub> concentration scale to single-phase CO2 density |
| $\varrho$ | density of single-phase CO<sub>2</sub> |
| $\mathcal{L}$ | length scale |
| $\mathcal{U}$ | velocity scale |
| $\mathcal{T}$ | time scale |

### Finite element method

| Symbol | Description |
| -------- | ------- |
| $\Omega\subset\mathbb{R}^d$ | domain |
| $d$ | number of spatial dimensions |
| $\partial\Omega$ | domain boundary |
| $\text{d}\Omega$ | integration measure over the cells | 
| $\text{d}\Gamma$ | integration measure over the cell facets | 
| $\textbf{n}$ | outward normal unit vector |
| $\mathscr{T}$ | tesselation of the domain | 
| $\mathcal{K}$ | mesh cell | 
| $h$ | mesh cell size | 
| $\mathcal{F}$ | set of cell facets |
| $\mathcal{V}$ | set of cell vertices |
| | | 
| $L^2(\Omega)$ | Lebesgue space |
| $H^1(\Omega)$ | Sobolev space |
| $\text{P}_k$ | continuous Lagrange polynomial of degree $k$ |
| $\text{BDM}_k$ | Brezzi-Douglas-Marini polynomial of degree $k$ |
| $\{\cdot\}$ | cell facet jump operator
| $\left[\!\left[\cdot\right]\!\right]$ | cell facet average operator |

### Finite difference method

| Symbol | Description |
| -------- | ------- |
| $\Delta t$ | timestep |
| $t^n$ | $n^{\text{th}}$ time level |
| $u^n=u(t^n)$| time-dependent quantity at the $n^{\text{th}}$ time level |

### Operators