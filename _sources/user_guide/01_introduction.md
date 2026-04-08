# Introduction


## Capillary-trapped CO<sub>2</sub>

![Figure](./figures/pore_scale_unreactive.png)

$$
\phi + \phi_{\text{CO}_2} + \phi_{\text{rock}} = 1 
$$

$$
\phi(\textbf{x}, t)=\varphi(\textbf{x})(1-s(\textbf{x},t)) 
$$

$$
\begin{align*}
s=0 &\iff \phi = \varphi \\
s=1 &\iff \phi = 0 \\
\frac{\partial\varphi}{\partial t}=0&\implies\frac{\partial\phi}{\partial t}=-\varphi\frac{\partial s}{\partial t} \\
\end{align*}
$$

## Notation

In addition to the notation defined in the [LUCiFEx notebooks](https://george-poole.github.io/LUCiFEx/notation.html).

| Symbol | Description |
| -------- | ------- |
| $s$ | saturation of capillary-trapped CO<sub>2</sub> | 
| $c$ | concentration of dissolved CO<sub>2</sub> | 
| $\theta$ | temperature | 
| $\textbf{u}$ | fluid velocity |
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