## Models






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
