# Novel convection-reaction systems in CRoCoDiL


User-defined constitutive relations, boundary conditions and initial conditions may be prescribed to investigate novel models in CRoCoDiL. Solutal or thermal transport may be 'switched off' by setting the solutal or thermal dispersion relation to `None`. Porosity evolution may be 'switched off' by setting the evolution number to `None`, or equivalently $\varepsilon=0$.

The creatively-named systems A, B and C are of interest to ongoing research.

## System A

Solutal convection-reaction in a closed, horizontal rectangle.

$$
\text{System A~}
\begin{cases}
\{\mathcal{L}, \mathcal{U}, \mathcal{T}\}_{\text{advective}} \\
\Omega=[0, L_x] \times [0, 1] \\ 
\rho(c)=c \\
R(s,c)=s(1-c) \\
s_0(y)=s_r\text{H}(y-h_0) \\
c_0(x,y)=c_r\text{H}(y-h_0) + \mathcal{N}(x,y)
\end{cases}
$$

## System B

Thermosolutal convection-reaction in a closed, horizontal rectangle with a heated lower boundary and a cooled upper boundary.

$$
\text{System B~}
\begin{cases}
\{\mathcal{L}, \mathcal{U}, \mathcal{T}\}_{\text{advective}} \\
\Omega=[0, L_x] \times [0, 1] \\ 
\rho(c)=c - \gamma\theta \\
R(s,c, \theta)=s(1+\delta\theta-c) \\
s_0=s_r \\
c_0(x,y)=f(y) + \mathcal{N}(x,y)\\
\theta_0(x,y)=1 - y + \mathcal{N}(x, y) \\
\theta_{\text{D}}(x,y=0)=1 \\
\theta_{\text{D}}(x, y=L_y)=0 \\
\theta_{\text{N}}(x=0,y)=0 \\
\theta_{\text{N}}(x=L_x, y)=0
\end{cases}
$$ 

## System C

Solutal convection-reaction in a closed, inclined rectangle.