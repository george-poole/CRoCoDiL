# Diagnostic equations

## Dimensional equations

| Name | Definition |
| -------- | ------- | 
| mass of capillary-trapped CO<sub>2</sub> | $m_C=\int_{\Omega}\varphi\varrho s~\text{d}\Omega$ | $\mathcal{L}$ | $\varphi_{\text{ref}}\mathcal{L}^d\Delta c$ |
| mass of dissolved CO<sub>2</sub> | $m_D=\int_{\Omega}\phi(\rho-\rho_{\mathrm{ref}})~\text{d}\Omega$ |
| total mass of CO<sub>2</sub> | $m=m_C + m_D$ |
| dissolution rate | $-\frac{\text{d}m_C}{\text{d}t}$ |
| dissolution flux | $F_{\Gamma}^m=-\frac{1}{\mathrm{len}(\Gamma)}\frac{\text{d}m_C}{\text{d}t}$ |
| solutal flux | $F_{\Gamma}=\frac{1}{\mathrm{len}(\Gamma)}\int_\Gamma(\textbf{u}c-\mathsf{D}\cdot\nabla c)\cdot\textbf{n}~\text{d}\Gamma = F_{\Gamma}^{\textbf{u}} + F_{\Gamma}^{\mathsf{D}}$ |
| thermal flux | $Q_{\Gamma}=\frac{1}{\mathrm{len}(\Gamma)}\int_\Gamma(\textbf{u}\theta-\mathsf{G}\cdot\nabla\theta)\cdot\textbf{n}~\text{d}\Gamma = Q_{\Gamma}^{\textbf{u}} + Q_{\Gamma}^{\mathsf{G}}$ |
| maximum speed | $\max_{\textbf{x}}\vert\textbf{u}\vert$ |
| root-mean-square speed | $\text{rms}(\textbf{u})=\Vert\textbf{u}\cdot\textbf{u}\Vert_{L_2(\Omega)}=\left(\int_\Omega\textbf{u}\cdot\textbf{u}~\text{d}\Omega\right)^{1/2}$ |
| velocity divergence norm | $\text{divnorm}(\textbf{u})=\Vert\nabla\cdot\textbf{u}\Vert_{L_2(\Omega)}=\left(\int_\Omega(\nabla\cdot\textbf{u})^2~\text{d}\Omega\right)^{1/2}$ |
| spatial average | $\langle c\rangle_{\Omega^\prime\subseteq\Omega} = \frac{1}{\text{vol}(\Omega^\prime)}\int_{\Omega^\prime} c~\text{d}\Omega$ |
| correction | $\mathcal{C}(c)=\max(c_{\text{min}}, \min(c, c_{\text{max}})) - c$ |


## Non-dimensionalization

| Quantity | $\text{d}\Omega$ | $\text{d}\Gamma$ | $\varphi, \phi$ | $m_C, m_D$ | $F$ | $Q$ |
| -------- | ------- | ------- | ------- | ------- | ------- | ------- | 
| **Scaling** | $\mathcal{L}^d$ | $\mathcal{L}$ | $\varphi_{\text{ref}}$ | $\varphi_{\text{ref}}\mathcal{L}^d\Delta c$ | ... | ... |


## Non-dimensional equations

| Name | Definition |
| -------- | ------- | 
| mass of capillary-trapped CO<sub>2</sub> | $m_C=\frac{1}{\varepsilon}\int_{\Omega} \varphi s~\text{d}\Omega$ |
| mass of dissolved CO<sub>2</sub> | $m_D=\int_{\Omega}\phi\rho~\text{d}\Omega$ |
| solutal flux | $F=...$ |
| thermal flux | $Q=...$ |