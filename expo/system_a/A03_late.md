# Late-time modelling

Subdomain averaging of the advection-diffusion-reaction and evolution equations obtains

$$
\begin{align*}
\frac{\text{d}c^+}{\text{d}t} &= -\frac{F(t)}{1-\zeta(t)} + \alpha{c} Da\frac{s^+(1-c^+)}{1-s^+} \\
\frac{\text{d}c^-}{\text{d}t} &= \frac{F(t)}{\zeta(t)} \\
\frac{\text{d}s^+}{\text{d}t} &= -\varepsilon \alpha_{s} Da s^+(1-c^+)
\end{align*}
$$

where due to the positive correlation of $s$ and $c$ in $\Omega^+$ the inequalities

$$
\begin{align*}
\frac{1}{\text{vol}(\Omega^+)}\int_{\Omega^+}\frac{s(1-c)}{1-s}~\text{d}\Omega &\leq \frac{s^+(1-c^+)}{1-s^+} \\
\frac{1}{\text{vol}(\Omega^+)}\int_{\Omega^+}s(1-c)~\text{d}\Omega &\leq s^+(1-c^+) \\
\end{align*}
$$

may be assumed to suppose that

$$
\begin{align*}
\int_{\Omega^+}\frac{s(1-c)}{1-s}~\text{d}\Omega &= \alpha_c\frac{s^+(1-c^+)}{1-s^+} \\
\int_{\Omega^+}s(1-c)~\text{d}\Omega &= \alpha_s s^+(1-c^+) \\
\end{align*}
$$

for some $\alpha_c, \alpha_s \leq 1$.

From total mass conservation, 

$$\underbrace{(1-\zeta_0)\left(\frac{s_r}{\varepsilon}+(1-s_r)c_r\right)}_{=m_0/\text{vol}(\Omega)}=(1-\zeta)(1-s^+)c^+ + \zeta c^- + (1-\zeta)\frac{s^+}{\varepsilon}$$

so

$$\zeta=\frac{\frac{m_0}{\text{vol}(\Omega)}-\frac{s^+}{\varepsilon}-(1-s^+)c^+}{c^- -\frac{s^+}{\varepsilon} - (1-s^+)c^+}$$

Empirically, let

$$F(c^+, c^-)=\alpha_1(c^+ - c^-) + \alpha_2(c^+ - c^-)^2$$

to close the ordinary differential equations, which may be solved numerically for $t>0$ with initial conditions  

$$
\begin{align*}
c^+(t=0)&=c_r \\
c^-(t=0)&=0 \\
s^+(t=0) &= s_r
\end{align*}
$$

or for $t>t_*$ with initial conditions  

$$
\begin{align*}
c^+(t=t_*) &=c^+_* \\
c^-(t=t_*)&=c^-_* \\
s^+(t=t_*) &= s^+_*
\end{align*}
$$

deduced from the direct numerical simulation data at $t=t_*$.

## Constrained equations

### Constant height constraint

If $\zeta=\zeta_0$ for all times, then instead $c^-$ can be eliminated via total mass conservation.

### Constant concentration constraint

If $c^+=1$ for all times, then $s^+$ is constant and only $c^-$ needs to be solved for.