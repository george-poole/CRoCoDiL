# Finite element method

## Mesh

$$\Omega\approx\cup_{\mathcal{K}\in\mathscr{T}}\mathcal{K}\subset\mathbb{R}^d$$

## Function spaces

### Mixed formulation

$$
V_{\textbf{u}} = \{\textbf{v}\ \in H(\text{div},\Omega)~:~(\textbf{n}\cdot\textbf{v})|_{\partial\Omega_{\text{E}}}=u_{\text{E}}~,~\textbf{v}|_{\mathcal{K}}\in\left[\text{BDM}_1(\mathcal{K})\right]^d ~\forall\mathcal{K}\in\mathscr{T}\}
$$

$$
V_p = \{v\in L^2(\Omega)~:~v|_{\mathcal{K}}\in\text{P}_0(\mathcal{K})~\forall\mathcal{K}\in\mathscr{T}\}
$$

### Streamfunction formulation

$$
V_\psi = \{v\in C^0(\Omega)~:~v|_{\partial{\Omega}_{\text{D}}}=\psi_{\text{D}}~,~v|_{\mathcal{K}}\in\text{P}_2(\mathcal{K})~\forall\mathcal{K}\in\mathscr{T}\}
$$

$$
V_{\textbf{u}} = \{\textbf{v}\in \left[C^0(\Omega)\right]^2~:~\textbf{v}|_{\partial{\Omega}_{\text{D}}}=\textbf{u}_{\text{E}}~,~\textbf{v}|_{\mathcal{K}}\in\left[\text{P}_1(\mathcal{K})\right]^2~\forall\mathcal{K}\in\mathscr{T}\}
$$

### Continuous Galerkin

$$
V_w^{\text{DG}} = \{v\in C^0(\Omega)~:~v|_{\partial{\Omega}_{\text{D}, w}}=w_{\text{D}}~,~v|_{\mathcal{K}}\in\text{P}_1(\mathcal{K})~\forall\mathcal{K}\in\mathscr{T}\}~,~w\in\{c,\theta\}
$$
### Discontinuous Galerkin

$$
V_w^{\text{DG}} = \{v\in L^2(\Omega)~:~v|_{\mathcal{K}}\in\text{P}_1(\mathcal{K})~\forall\mathcal{K}\in\mathscr{T}\}~,~w\in\{c,\theta\}
$$

## Weak forms

### Mixed formulation

$$
\begin{align*}
&\text{Find} \\
&(\textbf{u}^n, p^n)\in V_\textbf{u}\times V_p \\
&\theta^{n+1}\in V_\theta \\
&c^{n+1}\in V_c \\
&s^{n+1}\in V_s \\
&\text{such that} \\
&\mathbb{F} 
\begin{cases}
F_{\textbf{u},p}(\textbf{u}^n, p^n, \textbf{v}, q)=0 \quad\forall(\textbf{v}, q)\in V_{\textbf{u}} \times V_p \\
F_\theta(\theta^{n+1}, v) = 0\quad\forall v\in V_\theta \\ 
F_c(c^{n+1}, v) = 0\quad\forall v\in V_c \\ 
F_s(s^{n+1}, v) = 0\quad\forall v\in V_s \\ 
\end{cases}
\end{align*}
$$

### Streamfunction formulation

$$
\begin{align*}
&\text{Find} \\
&\psi^n\in V_\psi \\
&\textbf{u}^n\in V_\textbf{u} \\
&\theta^{n+1}\in V_\theta \\
&c^{n+1}\in V_c \\
&s^{n+1}\in V_s \\
&\text{such that} \\
&\mathbb{F} 
\begin{cases}
F_\psi(\psi^n, v) = 0\quad\forall v\in V_\psi \\
F_{\textbf{u}}(\textbf{u}^n, \textbf{v}) = 0\quad\forall \textbf{v}\in V_{\textbf{u}} \\ 
F_\theta(\theta^{n+1}, v) = 0\quad\forall v\in V_\theta \\ 
F_c(c^{n+1}, v) = 0\quad\forall v\in V_c \\ 
F_s(s^{n+1}, v) = 0\quad\forall v\in V_s \\ 
\end{cases}
\end{align*}
$$

### Variational forms

$$
\begin{align*}
F_{\textbf{u},p}(\textbf{u}^n, p^n, \textbf{v}, q) &= \dots \\
F_\psi(\psi^n, v) & =-\int_\Omega\text{d}\Omega~\nabla v\cdot\bigg(\frac{\mu^n\,(\mathsf{K}^n)^{\mathsf{T}}\cdot\nabla\psi^n}{\det(\mathsf{K}^n)}\bigg)+v\frac{\partial(\rho^n\textbf{e}_g\cdot\textbf{e}_y)}{\partial x} - v\frac{\partial(\rho^n\textbf{e}_g\cdot\textbf{e}_x)}{\partial y}=0 \\
F_{\textbf{u}}(\textbf{u}^n, \textbf{v}) &= \int_\Omega\text{d}\Omega~\textbf{v}\cdot\textbf{u}^n-\textbf{v}\cdot\begin{pmatrix}
-\frac{\partial\psi^n}{\partial y} \\
\frac{\partial\psi^n}{\partial x}
\end{pmatrix} \\
F_c(c^{n+1}, v) &\in F^{j}(\theta^{n+1}, v; Ad\,\textbf{u}, Di\mathsf{D}, Ki R, Ki J)~,~j\in\{\text{CG}, \text{DG}\}\\
F_\theta(\theta^{n+1}, v) &\in F^{j}(\theta^{n+1}, v; Ad\,\textbf{u}, \tfrac{Di}{Le}\mathsf{G}, 0, 0)~,~j\in\{\text{CG}, \text{DG}\}\\
F^{\text{CG}}(w^{n+1}, v; \textbf{u}, \mathsf{D}, R, J) &= \int_\Omega\text{d}\Omega~v\frac{w^{n+1}-w^n}{\Delta t^n} \\
&\quad +\int_\Omega\text{d}\Omega~v\frac{\mathcal{D}_{\textbf{u}, w}(\textbf{u}\cdot\nabla w)}{\mathcal{D}_\phi(\phi)} \\
&\quad +\int_\Omega\text{d}\Omega~\nabla\bigg(\frac{v}{\mathcal{D}_\phi(\phi)}\bigg)\cdot(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla \mathcal{D}_w(w)) \\
&\quad -\int_\Omega\text{d}\Omega~v\frac{\mathcal{D}_{R,w}(Rw) \\
+ \mathcal{D}_{J}(J)}{\mathcal{D}_\phi(\phi)} \\
&\quad - \int_{\partial\Omega_{\text{N},c}}\text{d}\Gamma~\frac{vw_{\text{N}}}{\mathcal{D}_\phi(\phi)} \\
&\quad + F^{\text{SUPG}}(w^{n+1}, v) \\
F^{\text{SUPG}}(w^{n+1}, v; \textbf{u}, \mathsf{D}, R, J) &= \tau^n(h, \textbf{u}_{\text{eff}}, \mathsf{D}_{\text{eff}}, R_{\text{eff}})(\nabla v\cdot \textbf{u}^{n}_{\text{eff}, w})\mathcal{R}(w^{n+1}, \textbf{u}, \mathsf{D}, R, J)  \\

\tau^n(h, \textbf{u}_{\text{eff}}, \mathsf{D}_{\text{eff}}, R_{\text{eff}}) &= \begin{cases}
(2|\textbf{u}_{\text{eff}}| / h  +  4D_{\text{eff}} / h^2  -  R_{\text{eff}})^{-1} & \text{Codina} \\
((2|\textbf{u}_{\text{eff}}| / h)^2  +  9(4D_{\text{eff}} / h^2)^2  +  R_{\text{eff}}^2)^{-1/2} & \text{Shakib} \\
((2 / \Delta t)^2 + (2|\textbf{u}_{\text{eff}}| / h)^2  +  (2D_{\text{eff}} / h^2)^2  +  R_{\text{eff}}^2)^{-1/2} & \text{transient} \\
(h/2|\textbf{u}_{\text{eff}}|)(\text{coth}(|\textbf{u}_{\text{eff}}|h/2D_{\text{eff}}) - 2D_{\text{eff}}/|\textbf{u}_{\text{eff}}|h) & \text{coth} \\
0 & \text{unstabilized} 
\end{cases} \\

\mathcal{R}(w^{n+1}; \textbf{u}, \mathsf{D}, R, J) &= \frac{w^{n+1}-w^n}{\Delta t^n} + \frac{\mathcal{D}_{\textbf{u},w}(\textbf{u}\cdot\nabla w)}{\mathcal{D}_\phi(\phi)}\\
&\quad -\frac{\nabla\cdot(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla \mathcal{D}_w(w))}{\mathcal{D}_\phi(\phi)} \\
&\quad  -\frac{\mathcal{D}_{R,w}(Rw) + \mathcal{D}_{J}(J)}{\mathcal{D}_\phi(\phi)} \\
F^{\text{DG}}(w^{n+1}, v; \textbf{u}, \mathsf{D}, R, J) &= \dots

% F_c(c^{n+1}, v)&\equiv
% \int_\Omega\text{d}\Omega~v\frac{c^{n+1}-c^n}{\Delta t^n} \\
% &\quad -Ad \int_\Omega\text{d}\Omega~\,\nabla\left(\frac{v}{\mathcal{D}_\phi(\phi)}\right)\cdot\nabla\mathcal{D}_{\textbf{u},c}(\textbf{u}c) \\
% &\quad + Ad\int_{\mathcal{F}/\partial\Omega}\text{d}\Gamma~2\llbracket v\rrbracket\textbf{n}\cdot\begin{cases}
% \{\mathcal{D}_{\textbf{u},c}(\textbf{u}c)\} & \textbf{n}\cdot\textbf{u} > 0 \\
% 0 & \text{otherwise} \\
% \end{cases} \\
% &\quad + Ad\int_{\partial\Omega}\text{d}\Gamma~\,\frac{v}{\mathcal{D}_\phi(\phi)}\begin{cases}
% 0 & \textbf{n}\cdot\textbf{u} >0 \\
% c_{\text{N}}\,\textbf{n}\cdot\textbf{u} & \text{otherwise} \\
% \end{cases} \\
% &\quad +\frac{1}{Pe}\int_\Omega\text{d}\Omega~\nabla\left(\frac{v}{\mathcal{D}_\phi(\phi)}\right)\cdot\nabla(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\mathcal{D}_c(c)) \\
% &\quad - \frac{1}{Pe}\int_{\mathcal{F}/\partial\Omega}\text{d}\Gamma~\left\{\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\left(\frac{v}{\mathcal{D}_\phi(\phi)}\right)\right\}\cdot\llbracket\mathcal{D}_c(c)\rrbracket\textbf{n}+\left\{\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\mathcal{D}_c(c)\right\}\cdot\llbracket \frac{v}{\mathcal{D}_\phi(\phi)}\rrbracket\textbf{n} \\
% &\quad +\frac{1}{Pe}\int_{\mathcal{F}/\partial\Omega}\text{d}\Gamma~ \frac{\alpha_{\text{DG}}}{h}\llbracket \frac{v}{\mathcal{D}_\phi(\phi)}\rrbracket\cdot\llbracket\mathcal{D}_c(c)\rrbracket\\
% &\quad -\frac{1}{Pe}\quad \int_{\partial\Omega_{\text{D},c}}\text{d}\Gamma~\left(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\left(\frac{v}{\mathcal{D}_\phi(\phi)}\right)\right)\cdot(\mathcal{D}_c(c)-c_{\text{D}})\textbf{n}\\

% &\quad -\frac{1}{Pe}\quad \int_{\partial\Omega_{\text{D},c}}\text{d}\Gamma~(\mathcal{D}_{\mathsf{D}}(\mathsf{D})\cdot\nabla\mathcal{D}_c(c))\cdot\left(\frac{v}{\mathcal{D}_\phi(\phi)}\textbf{n}\right)\\

% &\quad +\frac{1}{Pe}\quad \int_{\partial\Omega_{\text{D},c}}\text{d}\Gamma~\frac{\gamma_{\text{DG}}}{h} v(c-c_{\text{D}})\\
% &\quad -\frac{1}{Pe}\quad \int_{\partial\Omega_{\text{N},c}}\text{d}\Gamma~\frac{vc_{\text{N}}}{\mathcal{D}_\phi(\phi)}\\
% &\quad -Ki\int_\Omega\text{d}\Omega~v\,\frac{\mathcal{D}_{R,c}(Rc)}{\mathcal{D}_\phi(\phi)} \\

\end{align*}
$$

## Linear algebra 

$$
\textbf{u}^n=\sum_jU_j^n\boldsymbol{\xi}_j^{\textbf{u}}~,~
p^n=\sum_jP_j\xi_j^{p} ~,~
\psi^n=\sum_j\Psi_j^n\xi_j^{\psi} ~,~
\theta^n =\sum_j\Theta_j^n\xi_j^{\theta} ~,~
c^n=\sum_jC_j^n\xi_j^{c} ~,~
s^n =\sum_jS_j^n\xi_j^{s}
$$

$$\text{span}\{\xi_j^w\}_j\subset V_j\quad\forall w\in\{\textbf{u}, p,\psi,\theta,c,s\}$$

$$
\begin{align*}
&\text{Find}~\textbf{U}^n,~\textbf{P}^n,~\boldsymbol{\Theta}^{n+1},~\textbf{C}^{n+1},~\textbf{S}^{n+1}~\text{such that}~\forall n\geq0 \\
&\mathbb{A} 
\begin{cases}
\mathsf{A}^n_{\textbf{u},p}\cdot(\textbf{U}^n, \textbf{P}^n)^{\mathsf{T}}=\textbf{b}^n_{\textbf{u},p} \\
\mathsf{A}^{n}_{\theta}\cdot\boldsymbol{\Theta}^{n+1}\cdot\textbf{C}^n=\textbf{b}^n_{\theta} \\
\mathsf{A}^{n}_{c}\cdot\textbf{C}^{n+1}=\textbf{b}^n_{c} \\
\mathsf{A}^{n}_{s}\cdot\textbf{S}^{n+1}=\textbf{b}^n_{s} \\
\end{cases}
\end{align*}
$$