# Early-time modelling

## Time-dependent base state 

For any horizontally-invariant inital conditions $c_0(y)$ and $s_0(y)$, the system will initially in a horizontally-invariant and approximately motionless base state described by

$$
\begin{align*}
\textbf{u}&=\textbf{0} \\
p&=p_{\text{b}}(y,t) \\
c&=c_{\text{b}}(y,t) \\
s&=s_{\text{b}}(y,t)
\end{align*}
$$

until velocity is no longer negligible at the onset of convection. In the following we consider the $\varepsilon=0$ case, in effect neglecting the evolution of $s$ such that $s=s_0(y)$ for all $t$. To derive the concentration $c_{\text{b}}(y,t)$ base state, we solve

$$
\begin{align*}
\frac{\partial c_{\text{b}}}{\partial t}&=\frac{1}{Ra}\frac{\partial^2 c_{\text{b}}}{\partial y^2}-\frac{1}{Ra}\frac{\frac{\text{d}s_0}{\text{d}y}}{1-s_0(y)}\frac{\partial c_{\text{b}}}{\partial y}+Da\frac{s_0(y)}{1-s_0(y)}(1-c_{\text{b}})
\end{align*}
$$

subject to initial conditions $c_{\text{b}}(y,t=0)=c_0(y)$ and boundary conditions $\frac{\partial c_{\text{b}}}{\partial y}|_{y=0,1}=0$ by the method of separation of variables with $c_{\text{b}}(y,t)=1-Y(y)T(t)$. This leads to an eigenvalue problem

$$
\begin{align*}
\left(-\frac{\text{d}^2}{\text{d}y^2} + \frac{\frac{\text{d}s_0}{\text{d}y}}{1-s_0(y)}\frac{\text{d}}{\text{d}y} + RaDa\frac{s_0(y)}{1-s_0(y)}\right)Y_n(y)=\lambda_nY_n(y)
\end{align*}
$$

or equivalently in Sturm-Liouville form with weight function $1-s_0(y)$

$$
\begin{align*}
\bigg(-\frac{\text{d}}{\text{d}y}\bigg((1-s_0(y))\frac{\text{d}}{\text{d}y}\bigg) + RaDa\,s_0(y)\bigg)Y_n(y)=\lambda_n(1-s_0(y))Y_n(y)
\end{align*}
$$

to solve for the eigenvalues $\{\lambda_n\}$ and eigenfunctions $\{Y_n(y)\}$ satisfying boundary conditions $\frac{\text{d}Y}{\text{d}y}|_{y=0,1}=0$. It follows from the theory of self-adjoint operators [citation] that $\lambda_n\in\mathbb{R}$ and eigenfunctions corresponding to distinct eigenvalues are orthogonal with respect to the weighted inner product, so $\int_0^1 (1-s_0(y))Y_n(y)Y_{n'}(y)~\text{d}y=0$ for $n\neq n'$. The time-dependence follows the exponential decay $T_n(t)\propto e^{-\lambda_n t/Ra}$. By the principle of linear superposition, 

$$
\begin{align*}
c_{\text{b}}(y,t)=1-\sum_n c_n e^{-\lambda_n t/Ra}Y_n(y)
\end{align*}
$$

with constants $\{c_n\}$ determined from the initial condition $c_0(y)$ as

$$
\begin{align*}
c_n = \frac{\int_0^1(1-s_0(y))(1 - c_0(y))Y_n(y)~\text{d}y}{\int_0^1(1-s_0(y))|Y_n(y)|^2~\text{d}y}
\end{align*}
$$

Though we have in this section assumed $\varepsilon=0$ and hence a time-independent $s=s_0(y)$ to derive $c_{\text{b}}(y,t)$, a first-order approximation to the time-dependent $s_{\text{b}}(y,t)$ may still be found by substituting $c_{\text{b}}(y,t)$ in the saturation evolution equation

$$
\begin{align*}
s_{\text{b}}(y,t)\approx\exp\bigg(-\varepsilon DaRa\sum_n\frac{c_n}{\lambda_n}(s_n-e^{-\lambda_n t/Ra})Y_n(y)\bigg)
\end{align*}
$$

where $\{s_n\}$ are constants determined from the initial condition $s_0(y)$.

## Separation of variables solution

In general the eigenvalue problem is analytically intractable, however if $c_0(y)$ and $s_0(y)$ are Heaviside step functions, an exact solution may be pursued. The eigenvalue problem simplifies to 

$$
\begin{align*}
\lambda_nY_n(y) = \begin{cases}
(-\frac{\text{d}^2}{\text{d}y^2} + \Lambda)Y_n(y) & y>\zeta_0	 \\
-\frac{\text{d}^2Y_n}{\text{d}y^2} & y<\zeta_0
\end{cases}	
\end{align*}
$$

where for convenience we define the constant $\Lambda=RaDa\,s_r/(1-s_r)$, with eigenfunctions

$$
\begin{align*}
Y_n(y) &= \begin{cases}
A_n^+\cos((y-\zeta_0)\sqrt{\lambda_n-\Lambda}) + B_n^+\sin((y-\zeta_0)\sqrt{\lambda_n-\Lambda}) & y>\zeta_0 \\
A_n^-\cos((y-\zeta_0)\sqrt{\lambda_n}) + B_n^-\sin((y-\zeta_0)\sqrt{\lambda_n})  & y<\zeta_0
\end{cases}
\qquad\lambda_n>\Lambda\\
Y_n(y) &= \begin{cases}
A_n^+\cosh((y-\zeta_0)\sqrt{\Lambda-\lambda_n}) + B_n^+\sinh((y-\zeta_0)\sqrt{\Lambda-\lambda_n}) & y>\zeta_0 \\
A_n^-\cos((y-\zeta_0)\sqrt{\lambda_n}) + B_n^-\sin((y-\zeta_0)\sqrt{\lambda_n})  & y<\zeta_0
\end{cases}
\qquad\Lambda>\lambda_n>0\\
Y_n(y) &= \begin{cases}
A_n^+\cosh((y-\zeta_0)\sqrt{\Lambda+|\lambda_n|}) + B_n^+\sinh((y-\zeta_0)\sqrt{\Lambda+|\lambda_n|}) & y>\zeta_0 \\
A_n^-\cosh((y-\zeta_0)\sqrt{|\lambda_n|}) + B_n^-\sinh((y-\zeta_0)\sqrt{|\lambda_n|})  & y<\zeta_0
\end{cases}
\qquad\lambda_n<0
\end{align*}
$$

The upper and lower eigenfunctions must satisfy the boundary conditions $\frac{\text{d}Y}{\text{d}y}|_{y=1}=0$ and $\frac{\text{d}Y}{\text{d}y}|_{y=0}=0$ respectively, and additionally at $y=\zeta_0$ there must be continuity of $Y_n$ and $\frac{\text{d}Y}{\text{d}y}$. Without loss of generality we set $A_n^\pm=1$. Applying continiuity of $Y_n$, we find

$$
\begin{align*}
B_n^+&=
\begin{cases}
\tan((1-\zeta_0)\sqrt{\lambda_n - \Lambda}) & \lambda_n>\Lambda \\
\tanh((1-\zeta_0)\sqrt{\Lambda-\lambda_n }) & \Lambda>\lambda_n>0 \\
-\tanh((1-\zeta_0)\sqrt{C+|\lambda_n|}) & \lambda_n<0 \\
\end{cases} \\
B_n^-&=
\begin{cases}
-\tan(h\sqrt{\lambda_n}) & \lambda_n>\Lambda \\
-\tan(h\sqrt{\lambda_n}) & \Lambda>\lambda_n>0 \\
\tanh(h\sqrt{|\lambda_n|}) & \lambda_n<0 \\
\end{cases}
\end{align*}
$$

Applying continiuity of $\frac{\text{d}Y_n}{\text{d}y}$, we obtain the algebraic equations required to solve for the eigenvalue spectum $\{\lambda_n\}$

$$
\begin{align*}
\sqrt{|\lambda|} = \begin{cases}
-\sqrt{\lambda-\Lambda}\frac{\tan((1-\zeta_0)(\lambda-\Lambda))}{\tan(h\sqrt{\lambda})} & \lambda>\Lambda \\
\sqrt{\Lambda-\lambda}\frac{\tan((1-\zeta_0)(\Lambda-\lambda))}{\tan(h\sqrt{\lambda})} & \Lambda>\lambda>0 \\
-\sqrt{\Lambda+|\lambda|}\frac{\tan((1-\zeta_0)(\Lambda+|\lambda|))}{\tanh(h\sqrt{|\lambda|})} & \lambda<0
\end{cases}
\end{align*}
$$

which can be obtained by numerical root-finding methods. 

In the limit $Da\to\infty$, concentration in upper subdomain reaches its equilibrium value very quickly and the boundary layer above $y=\zeta_0$ is vanishing small, so an alternative solution can be obtained by the method of separation of variables. If we approximate $c(y>\zeta_0,t)\approx 1$ for all $t>\mathcal{O}((1-s_r)/Da\,s_r)$, we then need only solve for the concentration in the lower subdomain subject to boundary conditions $\frac{\partial c_{\text{b}}}{\partial t}(y=0, t)=0$ and $c_{\text{b}}(y=\zeta_0, t)= 1$. The solution obtained by the method of separation of variables is 

$$
c^-(y,t)=1-\sum_{n\in\mathbb{N}}\frac{4(-1)^n}{\pi}\cos\left(\frac{(2n+1)\pi y}{2\zeta_0}\right)e^{-\frac{(n+1)^2\pi^2t}{4\zeta_0^2Ra}}$$

which at times $t<\frac{1}{Ra}$ is well-approximated by the error function

$$
c^-(y,t)=1-\mathrm{erf}\left(\frac{y-\zeta_0}{2Ra\sqrt{t}}\right)
$$

## Similarity solution

Working with the same Heaviside step function initial conditions, we may also pursue a similarity solution by noting that the concentration base state is approximately constant in the upper subdomain, except within the vicinity of a small internal boundary layer around $y=\zeta_0$. If the dissolution reaction dominates in the upper subdomain, then we solve

$$\frac{\partial c_{\text{b}}}{\partial t}\approx\frac{Da\,s_r}{1-s_r}(1-c_{\text{b}})$$

to obtain 

$$c_{\text{b}}(y>\zeta_0, t) \approx c^+(t)\equiv 1-\exp\left(-\frac{Da\,s_r}{1-s_r}t\right)$$

valid for $t<\mathcal{O}((1-s_r)/Da\,s_r)$. The saturation at these early times is  

$$s(y>\zeta_0,t)\approx s^+(t)\equiv s_r\exp\left(\varepsilon\frac{1-s_r}{s_r}\left(e^{-\frac{Da\,s_r}{1-s_r}t}-1\right)\right) $$

and of course $s(y<\zeta_0,t)=0$ at all times. As the concentration in the upper subdomain tends towards equilibrium, diffusion will become more relevant, especially around $y=\zeta_0$, so we instead consider 

$$0\approx\frac{1}{Ra}\frac{\partial^2 c_{\text{b}}}{\partial y^2}+\frac{Da\,s_r}{1-s_r}(1-c_{\text{b}})$$

Since we expect a thin boundary layer near $y=\zeta_0$, we operate under the assumption that $RaDa\gg1$ and rescale the vertical coordinate as $Y=\sqrt{RaDa\,s_r/(1-s_r)}(y-\zeta_0)$ to solve the quasi-steady approximation

$$0=\frac{\partial^2 c_{\text{b}}}{\partial Y^2}+(1-c_{\text{b}})$$

with boundary conditions $\frac{\partial c_{\text{b}}}{\partial Y}|_{Y\to\infty}=0$ to obtain $c_{\text{b}}(Y,t)=1-(1-c_{\zeta_0})e^{-Y}$, where the interfacial concentration $c_{\zeta_0}$ is to be determined from continuity conditions between solutions obtained for the upper and lower subdomains.

In the lower subdomain, we pursue a similarity solution by transforming to the similarity variable $\eta=\sqrt{Ra}(y-\zeta_0)/t^{1/2}$ and seeking solutions of the form $c_{\text{b}}(y,t)=f(\eta)$ to

$$
\begin{align*}
\frac{\partial c_{\text{b}}}{\partial t}=\frac{1}{Ra}\frac{\partial^2 c_{\text{b}}}{\partial y^2} 
\implies \frac{\text{d}^2f}{\text{d}\eta^2} +\frac{\eta}{2}f(\eta)&=0
\end{align*}
$$

subject to the boundary condition $\frac{\text{d}f}{\text{d}\eta}(\eta\to-\infty)=0$ and continuity condition $f(\eta=0)=c_{\zeta_0}$ to obtain 

$$
\begin{align*}
f(\eta)=c_{\zeta_0}\text{erf}(\eta)+1) 
\implies c_{\text{b}}(y,t) = c_{\zeta_0}\bigg(\text{erf}\bigg(\frac{\sqrt{Ra}(y-\zeta_0)}{t^{1/2}}\bigg)+1\bigg)
\end{align*}
$$

Finally we match the interfacial concentration gradients $\frac{\partial c_{\text{b}}}{\partial y}|_{y=\zeta_0}$ to obain 

$$c_{\zeta_0}(t)=\frac{t^{1/2}}{t^{1/2}+\sqrt{\frac{4(1-s_r)}{\pi Da\,s_r}}}$$

and hence our full similarity solution for the concentration base state is 

$$
\begin{align*}
c_{\text{b}}(y,t)=\begin{cases}
1-(1-c_{\zeta_0}(t))e^{-\sqrt{RaDa\,s_r/(1-s_r)}} & y>\zeta_0 \\
c_{\zeta_0}(t)\bigg(\text{erf}\bigg(\frac{\sqrt{Ra}(y-\zeta_0)}{t^{1/2}}\bigg)+1\bigg) & y<\zeta_0
\end{cases}
\end{align*}
$$