# Definition

$$
\mathbb{S}_{s,c,\theta,\psi}\begin{cases}
\Omega=[0, \mathcal{A}X] \times [0, X] & \text{rectangular domain} \\
s_0(y)=s_r\text{H}(y-\zeta_0X) + \mathcal{N}_s(x,y) & \text{perturbed initial saturation} \\
c_0(x,y)=c_r\text{H}(y-\zeta_0X) + \mathcal{N}_c(x,y) & \text{perturbed initial concentration} \\
\theta_0(x,y)=\pm\text{H}(y-\zeta_0X) + \mathcal{N}_\theta(x,y) & \text{perturbed initial temperature} \\
c_{\text{N}}\vert_{\partial\Omega}=0 & \text{no-solutal-flux on entire boundary} \\
\theta_{\text{N}}\vert_{\partial\Omega}=0 & \text{no-thermal-flux on entire boundary} \\
u_{\text{E}}\vert_{\partial\Omega}=0\iff\psi_{\text{D}}\vert_{\partial\Omega}=0 & \text{no-penetration on entire boundary} \\
\varphi = 1 & \text{constant rock porosity} \\
\mathsf{K} = \phi^2\mathsf{I} & \text{isotropic permeability} \\
\mathsf{D} = \phi\mathsf{I} & \text{isotropic solutal dispersion}\\
\mathsf{G} = \phi\mathsf{I} & \text{isotropic thermal dispersion}\\
\mu = 1 & \text{constant viscosity} \\
\rho(c,\theta)=c - \gamma\theta & \text{linear density} \\
\Sigma(c, \theta, s)=s(1+\delta\theta)(1-c) & \text{linear reaction}  \\
\iff R(\theta, s)=-J(\theta, s)=-s(1+\delta\theta) \\
\textbf{e}_g = -\textbf{e}_y & \text{vertically downward gravity} \\ 
\end{cases}
$$ 

![System B at t=0](./figures/system_b_Omega_0.png)

## Notation

| Symbol | Description |
| -------- | ------- |
| $Ra^{\mathrm{solutal}}=Ra$ | solutal Rayleigh number | 
| $Ra^{\mathrm{thermal}}=\gamma Ra/Le$ | thermal Rayleigh number | 