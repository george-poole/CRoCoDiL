
## Algorithm

1. Initialize namespace $\mathcal{Q}=\{t, \Delta t, s\}$
1. Initialize $n\gets0$
1. Initialize $t^n\gets0$
1. **if** $\exists\psi$ **then**
1. &emsp; Update namespace $\mathcal{Q} \gets \mathcal{Q}\cup\{\psi, \textbf{u}\}$
1. **else**
1. &emsp; Update namespace $\mathcal{Q} \gets \mathcal{Q}\cup\{\textbf{u}, p\}$
1. **if** $\exists\theta$ **then**
1. &emsp; Update namespace $\mathcal{Q} \gets \mathcal{Q}\cup\{\theta\}$
1. &emsp; Initialize $\theta^n \gets \theta_0$
1. **if** $\exists c$ **then**
1. &emsp; Update namespace $\mathcal{Q} \gets \mathcal{Q}\cup\{c\}$
1. &emsp; Initialize $c^n \gets c_0$
1. **while** $t^n<t_{\text{stop}}$ **and** $n<n_{\text{stop}}$ **do**
1. &emsp; **if** $\exists\psi$ **then**
1. &emsp;&emsp; Update $\psi^n$ by solving ...
1. &emsp;&emsp; Update $\textbf{u}^n$ by projecting or interpolating ...
1. &emsp; **else**
1. &emsp;&emsp; Update $\textbf{u}^n, p^n$ by solving ...
1. &emsp; **if** $n<n_{\text{init}}$ **then**
1. &emsp;&emsp; Update $\Delta t^n\gets\Delta t_{\text{init}}$
1. &emsp; **else**
1. &emsp;&emsp; Update $\Delta t^n\gets\Delta t^n_{\text{CFLR}}$
1. &emsp; Update $t^{n+1}\gets t^n + \Delta t^n$
1. &emsp; **if** $\exists\theta$ **then**
1. &emsp;&emsp; Update $\theta^{n+1}$ by solving ... 
1. &emsp; **if** $\exists c$ **then**
1. &emsp;&emsp; Update $c^{n+1}$ by solving ...
1. &emsp; **if** $\varepsilon\neq 0$ **then**
1. &emsp;&emsp; Update $s^{n+1}$ by solving ...
1. &emsp; Update $q^{n-j}\gets q^{n-j+1}~\forall q\in\mathcal{Q}$ and $\forall j\geq 0$ 
1. &emsp; Update $n\gets n+1$
1. **end while**