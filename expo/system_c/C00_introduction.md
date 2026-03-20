# Theory

$$
\mathbb{S}_{s,c,\psi}\begin{cases}
\Omega=[0, \mathcal{A}X] \times [0, X] & \text{rectangular domain} \\
s_0(y)=s_r\text{H}(y-\zeta_0(x)) + \mathcal{N}_s(x,y) & \text{perturbed initial saturation} \\
c_0(x,y)=c_r\text{H}(y-\zeta_0(x)) + \mathcal{N}_c(x,y) & \text{perturbed initial concentration} \\
c_{\text{N}}\vert_{\partial\Omega}=0 & \text{no-flux on entire boundary} \\
u_{\text{E}}\vert_{\partial\Omega}=0\iff\psi_{\text{D}}\vert_{\partial\Omega}=0 & \text{no-penetration on entire boundary} \\
\varphi = 1 & \text{constant rock porosity} \\
\mathsf{K} = \phi^2\begin{pmatrix}
\cos^2\vartheta +\kappa\sin^2\vartheta & (1 -  \kappa)\cos\vartheta\sin\vartheta \\
(1 -  \kappa)\cos\vartheta\sin\vartheta & \cos^2\vartheta +\kappa\sin^2\vartheta\\
 \end{pmatrix} & \text{anisotropic cross-bedded permeability} \\ 
\mathsf{D} = \phi\mathsf{I} & \text{isotropic dispersion}\\
\mu = 1 & \text{constant viscosity} \\
\rho(c)=c & \text{linear density} \\
\Sigma(c, s)=s(1-c) & \text{linear reaction}  \\
\iff R(s)=-J(s)=-s \\
\textbf{e}_g = -\sin\beta\textbf{e}_x -\cos\beta\textbf{e}_y & \text{inclined plane} \\ 
\end{cases}
$$ 

![System C at t=0](./figures/system_c_Omega_0.png)

## Notation

| Symbol | Description |
| -------- | ------- |
| $\beta$ | inclination angle | 
| $\kappa$ | permeability cross-bedding anisotropy | 
| $\vartheta$ | permeability cross-bedding angle | 
| $\zeta_0(x) = (\mathcal{A}X - x)\tan\beta + \zeta_0^{\mathrm{min}}X$| initial front height |