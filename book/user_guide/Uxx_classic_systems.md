# Classic convection systems in CRoCoDiL

For demonstration and benchmarking purposes, some classic convection systems are implemented in CRoCoDiL.

## Rayleigh-Bénard

$$
\text{Rayleigh-Bénard~}
\begin{cases}
\{\mathcal{L}, \mathcal{U}, \mathcal{T}\}_{\text{advective}} \\
\Omega=[0, L_x] \times [0, 1] \\ 
\rho(\theta)=-\theta \\
\theta_0(x,y)=1-y + \mathcal{N}(x,y) \\
\theta_{\text{D}}(x,y=0)=1 \\
\theta_{\text{D}}(x, y=1)=0 \\
\theta_{\text{N}}(x=0,y)=0 \\
\theta_{\text{N}}(x=L_x, y)=0
\end{cases}
$$

## Rayleigh-Taylor 

$$
\text{Rayleigh-Taylor~}
\begin{cases}
\{\mathcal{L}, \mathcal{U}, \mathcal{T}\}_{\text{advective}} \\
\Omega=[0, L_x] \times [0, 1] \\ 
\rho(\theta)=-\theta \\
\theta_0(x,y)=\text{H}(\tfrac{1}{2} - y) \\
\end{cases}
$$