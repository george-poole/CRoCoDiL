# Introduction

$$
\begin{cases}
\Omega=[0, \mathcal{A}X] \times [0, X] & \text{rectangular domain} \\
s_0(y)=s_r + \mathcal{N}_s(x,y) & \text{perturbed initial saturation} \\
\theta_0(y)=1-y+ \mathcal{N}_{\theta}(x,y) & \text{perturbed initial temperature} \\
c_0(x,y)=c_{\text{b}}(y)+ \mathcal{N}_{c}(x,y) & \text{perturbed initial concentration} \\
\theta_{\text{D}}(x,y=0)=1 & \text{hot lower boundary} \\
\theta_{\text{D}}(x, y=L_y)=0 & \text{cold upper boundary} \\
\theta_{\text{N}}(x=0,y)=0 & \text{no-thermal-flux on left boundary} \\
\theta_{\text{N}}(x=L_x, y)=0 & \text{no-thermal-flux on right boundary} \\
c_{\text{N}}\vert_{\partial\Omega}=0 & \text{no-solutal-flux on entire boundary} \\
u_{\text{E}}\vert_{\partial\Omega}=0\iff\psi_{\text{D}}\vert_{\partial\Omega}=0 & \text{no-penetration on entire boundary} \\
\varphi = 1 & \text{constant rock porosity} \\
\mathsf{K} = \phi\mathsf{I} & \text{isotropic permeability} \\
\mathsf{D} = \phi\mathsf{I} & \text{isotropic solutal dispersion}\\
\mathsf{G} = \phi\mathsf{I} & \text{isotropic thermal dispersion}\\
\mu = 1 & \text{constant viscosity} \\
\rho(c)=c -  \gamma\theta & \text{linear density} \\
\Sigma(c, s)=s(1+\delta\theta-c) & \text{linear reaction} \\
\textbf{e}_g = -\textbf{e}_y & \text{vertically downward gravity} \\ 
\end{cases}
$$ 

![System B at t=0](./figures/system_b_Omega0.png)

![System B at t>0](./figures/system_b_Omega.png)

## Concentration base state

If the evolution of $s$ is negligible, then

$$
Di\frac{\text{d}}{\text{d}y}\left((1-s_r)\frac{\text{d}c_{\text{b}}}{\text{d}y}\right) + Ki\,s_r(1+\delta-\delta y - c_{\text{b}})= 0
$$

The boundary conditions $\frac{\text{d}c_{\text{b}}}{\text{d}y}\vert_{y=0}=\frac{\text{d}c_{\text{b}}}{\text{d}y}\vert_{y=X}=0$ give the base state

$$
c_{\text{b}}(y) = 1 + \delta(1-y) + \frac{\delta}{\sqrt{\Lambda}}\frac{\sinh\left(\sqrt{\Lambda}(y-X/2)\right)}{\cosh\left(\sqrt{\Lambda}X/2\right)}
$$

where $\Lambda=\frac{Ki\,s_r}{Di(1-s_r)}$. The reaction term evaluates to 

$$
\Sigma(c_{\text{b}}, s_r) = \frac{\delta}{\sqrt{\Lambda}}\frac{\sinh\left(\sqrt{\Lambda}(y-X/2)\right)}{\cosh\left(\sqrt{\Lambda}X/2\right)}
$$

so the saturation evolution is of order $\mathcal{O}(\varepsilon Da\delta/\sqrt{\Lambda})$, justifying its negligibility.