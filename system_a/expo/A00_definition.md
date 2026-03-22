# Definition

$$
\mathbb{S}_{s,c,\psi}\begin{cases}
\Omega=[0, \mathcal{A}X] \times [0, X] & \text{rectangular domain} \\
s_0(y)=s_r\text{H}(y-\zeta_0X) + \mathcal{N}_s(x,y) & \text{perturbed initial saturation} \\
c_0(x,y)=c_r\text{H}(y-\zeta_0X) + \mathcal{N}_c(x,y) & \text{perturbed initial concentration} \\
c_{\text{N}}\vert_{\partial\Omega}=0 & \text{no-flux on entire boundary} \\
u_{\text{E}}\vert_{\partial\Omega}=0\iff\psi_{\text{D}}\vert_{\partial\Omega}=0 & \text{no-penetration on entire boundary} \\
\varphi = 1 & \text{constant rock porosity} \\
\mathsf{K} = \phi^2\mathsf{I} & \text{isotropic permeability} \\
\mathsf{D} = \phi\mathsf{I} & \text{isotropic dispersion}\\
\mu = 1 & \text{constant viscosity} \\
\rho(c)=c & \text{linear density} \\
\Sigma(c, s)=s(1-c) & \text{linear reaction} \\
\iff R(s)=-s\quad J(s)=s \\
\textbf{e}_g = -\textbf{e}_y & \text{vertically downward gravity} \\ 
\end{cases}
$$

![System A at t=0](./figures/system_a_Omega_0.png)

![System A at t>0](./figures/system_a_Omega_t.png)

## Notation

| Symbol | Description |
| -------- | ------- |
| $\mathcal{A}$ | domain aspect ratio | 
| $s_r$ | initial residual saturation |
| $c_r$ | initial residual concentration |
| $\zeta_0$ | initial front height |
| $\zeta$ | front height |
| $\Omega^+$ | upper subdomain containing capillary-trapped  CO<sub>2</sub> |
| $\Omega^-$ | lower subdomain without capillary-trapped  CO<sub>2</sub> |
| $\langle\cdot\rangle_x$ | horizontal average | 
| $\cdot^+=\langle\cdot\rangle_{\Omega^+}$ | upper subdomain average | 
| $\cdot^-=\langle\cdot\rangle_{\Omega^-}$ | upper subdomain average | 