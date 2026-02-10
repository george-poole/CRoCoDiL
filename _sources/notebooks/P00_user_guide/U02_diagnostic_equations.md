# Diagnostic equations

## Dimensional equations

| Name | Definition |
| -------- | ------- | 
| mass of capillary-trapped CO<sub>2</sub> | $m_C=\int_{\Omega}\varphi\varrho s~\text{d}\Omega$ | $\mathcal{L}$ | $\varphi_{\text{ref}}\mathcal{L}^d\Delta c$ |
| mass of dissolved CO<sub>2</sub> | $m_D=\int_{\Omega}\phi(\rho(c,\theta)-\rho(0,\theta))~\text{d}\Omega$ |
| total mass of CO<sub>2</sub> | $m=m_C + m_D$ |
| solutal flux | $F=\frac{1}{\vert\Gamma\vert}\int_\Gamma(\textbf{u}c-\nabla\cdot(\mathsf{D}\cdot\nabla c))\cdot\textbf{n}~\text{d}\Gamma = F_{\textbf{u}} + F_{\mathsf{D}}$ |
| thermal flux | $Q=\frac{1}{\vert\Gamma\vert}\int_\Gamma(\textbf{u}\theta-\nabla\cdot(\mathsf{G}\cdot\nabla \theta))\cdot\textbf{n}~\text{d}\Gamma = Q_{\textbf{u}} + Q_{\mathsf{G}}$ |
| maximum speed | $\max_{\textbf{x}}\vert\textbf{u}\vert$ |
| root-mean-square speed | $\text{rms}(\textbf{u})=\Vert\textbf{u}\cdot\textbf{u}\Vert_{L_2(\Omega)}=\left(\int_\Omega\textbf{u}\cdot\textbf{u}~\text{d}\Omega\right)^{1/2}$ |
| velocity divergence norm | $\text{divnorm}(\textbf{u})=\Vert\nabla\cdot\textbf{u}\Vert_{L_2(\Omega)}=\left(\int_\Omega(\nabla\cdot\textbf{u})~\text{d}\Omega\right)^{1/2}$ |
| solutal correction | $\mathcal{C}(c)=\max(c_{\text{min}}, \min(c, c_{\text{max}})) - c$ |
| thermal correction | $\mathcal{C}(\theta)=\max(\theta_{\text{min}}, \min(\theta, \theta_{\text{max}})) - \theta$ |


## Non-dimensionalization

| Quantity | $\text{d}\Omega$ | $\text{d}\Gamma$ | $m_C, m_D$ | $F$ | $Q$ |
| -------- | ------- | ------- | ------- | ------- | ------- | 
| **Scaling** | $\mathcal{L}^d$ | $\mathcal{L}$ | $\varphi_{\text{ref}}\mathcal{L}^d\Delta c$ | ... | ... |


## Non-dimensional equations

| Name | Definition |
| -------- | ------- | 
| mass of capillary-trapped CO<sub>2</sub> | $m_C=\frac{1}{\varepsilon}\int_{\Omega} s~\text{d}\Omega$ |
| mass of dissolved CO<sub>2</sub> | $m_D=\int_{\Omega} \phi c~\text{d}\Omega$ |
| solutal flux | $F=...$ |
| thermal flux | $Q=...$ |