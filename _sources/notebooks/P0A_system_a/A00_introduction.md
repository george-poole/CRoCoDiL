# Introduction

$$
\mathbb{S}\begin{cases}
\Omega=[0, \mathcal{A}X] \times [0, X] & \text{rectangular domain} \\
\textbf{e}_g = -\textbf{e}_y & \text{vertically downward gravity} \\ 
\varphi = 1 & \text{constant rock porosity} \\
\mathsf{K} = \phi\mathsf{I} & \text{isotropic permeability} \\
\mathsf{D} = \phi\mathsf{I} & \text{isotropic dispersion}\\
\mu = 1 & \text{constant viscosity} \\
\rho(c)=c & \text{linear density} \\
\Sigma(c, s)=s(1-c) \iff R(s)=-J(s)=-s & \text{linear reaction} \\
s_0(y)=s_r\text{H}(y-\mathcal{A}\zeta_0) + \mathcal{N}_s(x,y) & \text{perturbed initial saturation} \\
c_0(x,y)=c_r\text{H}(y-\mathcal{A}\zeta_0) + \mathcal{N}_c(x,y) & \text{perturbed initial concentration} \\
c_{\text{N}}\vert_{\partial\Omega}=0 & \text{no-flux on entire boundary} \\
u_{\text{E}}\vert_{\partial\Omega}=0\iff\psi_{\text{D}}\vert_{\partial\Omega}=0 & \text{no-penetration on entire boundary} \\
\end{cases}
$$

![System A initial conditions](./figures/Omega0.png)

![System A time evolution](./figures/Omega(t).png)

## Notation

| Symbol | Description |
| -------- | ------- |
| $s_r$ | initial residual saturation |
| $c_r$ | initial residual concentration |
| $\zeta_0$ | initial front height |
| $\zeta$ | front height |
| $\Omega^+$ | upper subdomain containing capillary-trapped  CO<sub>2</sub> |
| $\Omega^-$ | lower subdomain without capillary-trapped  CO<sub>2</sub> |
| $\mathcal{A}$ | domain aspect ratio | 