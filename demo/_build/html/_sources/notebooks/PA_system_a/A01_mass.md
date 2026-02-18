# Dissolved and capillary-trapped mass

## Conserved mass

From the initial conditions, the initial masses are

$$
\begin{align*}
m_C^{0}&=\frac{1}{\varepsilon}\int_\Omega \varphi s_0~\text{d}\Omega \\
m_D^{0}&=\int_\Omega \varphi(1-s_0)c_0~\text{d}\Omega
\end{align*}
$$

If $\varphi=1$, $s_0 = s_r\text{H}(y-\zeta_0X) + \mathcal{N}_s$ and $c_0 = c_r\text{H}(y-\zeta_0X)+ \mathcal{N}_c$ on the domain $\Omega = [0, \mathcal{A}X] \times [0, X]$, then

$$
\begin{align*}
m_C^0 &= \frac{1}{\varepsilon}\mathcal{A}X^2(1-\zeta_0)\left(s_r + \mathcal{O}(\mathcal{N}_s)\right) + \frac{1}{\varepsilon}\mathcal{A}X^2\zeta_0\mathcal{O}(\mathcal{N}_s)\\
m_D^0 &= \mathcal{A}X^2(1-\zeta_0)\left((1-s_r - \mathcal{O}(\mathcal{N}_s))(c_r + \mathcal{O}(\mathcal{N}_c))\right) + \mathcal{A}X^2\zeta_0(1-\mathcal{O}(\mathcal{N}_s))\mathcal{O}(\mathcal{N}_c)
\end{align*}
$$

If $\mathcal{N}_s=\mathcal{N}_c=0$, then

$$
\begin{align*}
m_C^0 &= \frac{1}{\varepsilon}\mathcal{A}X^2(1-\zeta_0)s_r\\
m_D^0 &= \mathcal{A}X^2(1-\zeta_0)(1-s_r)c_r \\
\end{align*}
$$

so the total conserved mass is

$$
m = \mathcal{A}X^2(1-\zeta_0)\left(\frac{s_r}{\varepsilon}+(1-s_r)c_r\right)
$$

## Complete dissolution

$c\to c_\infty(\textbf{x})$ and $s\to0$ as $t\to\infty$

$$
\begin{align*}
m_C^{\infty} & = 0 \\
m_D^{\infty} & = \int_\Omega \varphi c_\infty~\text{d}\Omega
\end{align*}
$$

If $\varphi=1$, then $m_D^{\infty}=\text{vol}(\Omega)\langle c_\infty\rangle$ so by mass conservation on $\Omega = [0, \mathcal{A}X] \times [0, X]$ 

$$(1-\zeta_0)\left(\frac{s_r}{\varepsilon}+(1-s_r)c_r\right)=\langle c_\infty\rangle$$

Due to the physical constraint $\langle c_\infty\rangle\leq 1$ the condition on complete dissolution is obtained as

$$s_r\leq \frac{\varepsilon\left(\frac{1}{1-\zeta_0}-c_r\right)}{1-\varepsilon c_r}\equiv s_r^*$$

## Incomplete dissolution

$c\to 1$ and $s\to s_\infty(\textbf{x})$ as $t\to\infty$

$$
\begin{align*}
m_C^{\infty} & = \frac{1}{\varepsilon}\int_\Omega \varphi s_\infty~\text{d}\Omega \\
m_D^{\infty} & = \int_\Omega \varphi(1 - s_\infty)~\text{d}\Omega = \text{vol}(\Omega)\langle\varphi\rangle -\varepsilon m_C^{\infty}
\end{align*}
$$

By mass conservation $m=m_C^{\infty}+m_D^{\infty}$ so

$$
\begin{align*}
m_C^{\infty} &= \frac{m - \text{vol}(\Omega)\langle\varphi\rangle}{1-\varepsilon} \\
m_D^{\infty} &= \frac{m - \frac{\text{vol}(\Omega)\langle\varphi\rangle}{\varepsilon}}{1-\frac{1}{\varepsilon}}
\end{align*}
$$

With $\varphi=1$ and $\Omega = [0, \mathcal{A}X] \times [0, X]$, $\langle\varphi\rangle=1$ and $\text{vol}(\Omega)=\mathcal{A}X^2$.

Using the previous expression for the total conserved mass and noting the physical constraint $m_C^{\infty}\geq0$ obtains the same expression for the critical residual saturation $s_r^*$.