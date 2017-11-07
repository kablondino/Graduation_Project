# Full Model Notes from Tim Stap's Thesis

The full model acts on the domain of
$$\Omega \,=\, \left\{x, t \,\in\, \mathbb{R}^2 \,|\, (0 \leq x \leq L) ~\text{and}~ (t \geq 0)\right\}$$

The model can be reduced down to the following form:
$$\dfrac{\partial\mathbf{v}}{\partial t} \,=\, \dfrac{\partial}{\partial x} F\left(x, t, \mathbf{v}, \dfrac{\partial\mathbf{v}}{\partial x}\right) + S\left(x, t, \mathbf{v}, \dfrac{\partial\mathbf{v}}{\partial x}\right)$$

The vectors $\mathbf{v}$, $F$, and $S$ represent the following:
$$\mathbf{v}(x, t) \,=\,\begin{bmatrix} n(x, t) \\ T(x, t) \\ Z(x, t) \end{bmatrix}~,~~~~~~~~~~
\mathbf{F} \,=\, \begin{bmatrix}
			D(\partial_x Z)\,\dfrac{\partial n}{\partial x} \\
			\dfrac{D(\partial_x Z)}{\zeta}\,\dfrac{\partial T}{\partial x} \\
			\dfrac{\mu D(\partial_x Z)}{\epsilon}\,\dfrac{\partial Z}{\partial x}
			\end{bmatrix}~,~~~~~~~~~~
\mathbf{S} \,=\, \begin{bmatrix}
			0 \\
			\left(\dfrac{\zeta + 1}{\zeta}\right) \dfrac{D(\mathcal{E})}{n} \dfrac{\partial n}{\partial x} \dfrac{\partial T}{\partial x} \\
			\dfrac{c_n T}{\epsilon n^2} \dfrac{\partial n}{\partial x} + \dfrac{c_T}{\epsilon n} \dfrac{\partial T}{\partial x} + \dfrac{G(Z)}{\epsilon}
			\end{bmatrix}.
$$

The diffusivity function $D(\mathcal{E})$ is given as
$$D(\partial_x Z) \,=\, D_\text{min} + \dfrac{D_\text{max} - D_\text{min}}{1 + \alpha_\text{sup}\cdot(\partial_x Z)^2},$$
with $\alpha_\text{sup}$ as the suppression rate coefficient. This arises from the nonlinear suppression model (2.3).

### Domain Boundary and Boundary Conditions

$$\delta \Omega \,=\, \left\{x, t \,\in\, \Omega ~|~ x = 0 ~~\text{and}~~ x = L\right\}$$

The Neumann and Robin boundary conditions can be expressed in the form of
$$p\left(x, t, \mathbf{v}\right) + F\left(x, t, \mathbf{v}, \dfrac{\partial\mathbf{v}}{\partial x}\right) \,=\, 0 ~~~~~ \text{for} ~~~~~ (x, t) \in \delta\Omega$$

The $p$ term is expressed, with $\Gamma_c$ and $q_c(t)$ as actuation paramters:
$$p(0, t) \,=\, -\begin{bmatrix}
				D(\partial_x Z) \dfrac{n}{\lambda_n}\\
				\dfrac{D(\partial_x Z)}{\zeta} \dfrac{T}{\lambda_T} \\
				\dfrac{\mu D(\partial_x Z)}{\epsilon} \dfrac{Z}{\lambda_Z}
				\end{bmatrix}_{x = 0}
~~~~~ \text{and} ~~~~~
p(L, t) \,=\, \begin{bmatrix}
				\Gamma_c(t) \\
				\dfrac{(\gamma - 1) q_c - T\Gamma_c(t)}{n} \\
				0
				\end{bmatrix}_{x = L}$$

A good initial condition for $Z$:
$$Z(x, 0) \,=\, Z_S\left[1 - \tanh\left(\dfrac{L\,x - L)}{2}\right)\right] \,=\, Z_S\left[1 - \frac{\exp(L\,x - L) - 1}{\exp(L\,x - L) + 1}\right]$$

<!---
Parameters used by Staps:
 $\Gamma_c$  $\gamma$  $\lambda_n$  $\lambda_T$  $\lambda_Z$  $D_\text{min}$  $D_\text{max}$  $\zeta$  $c_n$  $c_T$  $q_c$  $a$  $b$  $c$  $Z_S$  $\alpha_\text{sup}$  $\mu$  $\epsilon$ 
 -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  ----- 
 $-\dfrac{4}{5}$  $\dfrac{5}{3}$  $\dfrac{5}{4}$  $\dfrac{3}{2}$  $\dfrac{5}{4}$  $\dfrac{2}{5}$  $2$  $\dfrac{1}{2}$  $-1.1$  $-0.9$  $-4$  $\dfrac{3}{2}$  $2$  $-1$  $-\dfrac{3}{2}$  $\dfrac{1}{2}$  $\dfrac{1}{20}$  $\dfrac{1}{25}$ 
--->
-----------------------------------------------------------

## Expanded Model for All Plasma Variables
\begin{align}
	\dfrac{\partial n}{\partial t} \,&=\, \dfrac{\partial}{\partial x}\left(D\left(\partial_x Z\right)\dfrac{\partial n}{\partial x}\right) + 0 \\
	\dfrac{\partial T}{\partial t} \,&=\, \dfrac{\partial}{\partial x}\left(\dfrac{D\left(\partial_x Z\right)}{\zeta} \dfrac{\partial T}{\partial x}\right) + \left(\dfrac{\zeta + 1}{\zeta}\right) \dfrac{D(\mathcal{E})}{n} \dfrac{\partial n}{\partial x} \dfrac{\partial T}{\partial x} \\
	\epsilon \dfrac{\partial Z}{\partial t} \,&=\, \dfrac{\partial}{\partial x}\left(\mu\,D\left(\partial_x Z\right) \dfrac{\partial Z}{\partial x}\right) + \dfrac{c_n T}{n^2} \dfrac{\partial n}{\partial x} + \dfrac{c_T}{n} \dfrac{\partial T}{\partial x} + G(Z)
\end{align}

## $Z$ Model, with Substitutions
\begin{align}
	&\dfrac{m_i}{e \rho_{pi}} \,n T\, \left(\dfrac{B_\theta}{B}\right)^2 \dfrac{\partial Z}{\partial t} \,=\, \dfrac{m_i \mu_i}{e \rho_{pi}} \,n T\, \dfrac{\partial^2 Z}{\partial x^2} \\
	&+\, B_\theta^2 \left[\left(g_n^\text{an} - g_n^\text{cx} - g_n^{\pi\parallel}\right) \dfrac{1}{n} \dfrac{\partial n}{\partial x} + \left(g_T^\text{an} - g_T^\text{cx} - g_T^{\pi\parallel}\right) \dfrac{1}{T} \dfrac{\partial T}{\partial x} + \left(g_Z^\text{an} - g_Z^\text{cx} - g_Z^{\pi\parallel}\right) Z - f^\text{OL}\right]
\end{align}

### Terms in $Z$ Model

+ Electron Anomalous Diffusion

$$D^\text{an} \,=\, \dfrac{\epsilon^2 \sqrt{\pi}}{2 a_m} \dfrac{\rho_{pe} T}{B}$$

$$g_n^\text{an} \,=\, -e \,n\, D^\text{an}~,~~~~~~~~~~ g_T^\text{an} \,=\, -e \,n\, \alpha^\text{an}\, D^\text{an}~,~~~~~~~~~~ g_Z^\text{an} \,=\, \dfrac{-e \,n\, D^\text{an}}{\rho_{pi}}$$

+ Charge Exchange Friction

$$g_n^\text{cx} \,=\, -\dfrac{m_i \,n_0 \langle\sigma v\rangle_\text{cx} \,n T}{B_\theta^2}~,~~~~~~~~~~ g_T^\text{cx} \,=\, \alpha^\text{cx}\,g_n^\text{cx}~,~~~~~~~~~~ g_Z^\text{cx} \,=\, -\dfrac{g_n^\text{cx}}{\rho_{pi}}$$

+ Ion Orbit Loss

$$g^\text{OL} \,=\, e \,n\, \nu_\text{eff} \sqrt{\epsilon} \,\rho_{pi}~,~~~~~~~~~~ f^\text{OL} \,=\, \dfrac{g^\text{OL}\,\exp\left[-\sqrt{\nu_{*i} + Z^4}\right]}{\sqrt{\nu_{*i} + Z^4}}$$

+ Ion Bulk Viscosity

$$N \,=\, \dfrac{\nu_{*i}\,\epsilon^{3/2}\,\nu_{ei}}{\nu_{ii}}$$

\begin{align}
	\begin{pmatrix}\xi_\theta \\ \xi_\phi \end{pmatrix} \,&=\, \dfrac{1}{\pi} \int_0^{\sqrt{\nu_{*i}}} \begin{pmatrix} 1 \\ \frac{5}{2} - x \end{pmatrix} x^2 \exp(-x) \left[\int_{-1}^{+1} \dfrac{N / \sqrt{x} ~~ \text{d}y}{\left(y + Z / \sqrt{x}\right)^2 + \left(N / \sqrt{x}\right)^2}\right] \text{d}x \\
	&=\, \dfrac{1}{\pi} \int_0^{\sqrt{\nu_{*i}}} \begin{pmatrix} 1 \\ \frac{5}{2} - x \end{pmatrix} x^2 \exp(-x) \, \tan^{-1}\left(\dfrac{2 N \sqrt{x}}{N^2 + Z^2 - x}\right) \text{d}x
\end{align}

$$\text{To consolidate, use the following:}~~~\eta \,=\, \dfrac{\epsilon^2 \sqrt{\pi}}{8 a_m} m_i \,n\, (v_{T_i})^2$$

$$g_n^{\pi\parallel} \,=\, \eta \, \rho_{\pi} B_\theta \, \xi_\theta~,~~~~~~~~~~ g_n^{\pi\parallel} \,=\, \eta \, \rho_{pi} \left(B_\theta\,\xi_\theta - B\,\xi_\phi\right)~,~~~~~~~~~~ g_Z^{\pi\parallel} \,=\, 2\eta \, B_\theta \, \xi_\theta $$

