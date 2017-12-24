# Full Model Notes

The full model acts on the domain of
\begin{align}
	\Omega \,=\, \left\{x, t \,\in\, \mathbb{R}^2 \,|\, (0 \leq x \leq L) ~\text{and}~ (t \geq 0)\right\}.
\end{align}

The model can be reduced down to the following form:
\begin{align}
	\dfrac{\partial}{\partial t} \mathbf{v}&(x,t) \,=\, \dfrac{\partial}{\partial x} F\left(x, t, \mathbf{v}, \dfrac{\partial\mathbf{v}}{\partial x}\right) \,+\, S\left(x, t, \mathbf{v}, \dfrac{\partial\mathbf{v}}{\partial x}\right) \\
\mathbf{v}& \,=\,\begin{bmatrix} n(x, t) \\ T(x, t) \\ Z(x, t) \end{bmatrix}~,~~~~
\mathbf{F} \,=\, \begin{bmatrix}
			D\,\cdot n^\prime \\[1ex]
			\dfrac{D}{\zeta}\,\cdot T^\prime \\[2ex]
			\dfrac{\mu D}{\epsilon}\,\cdot Z^\prime
			\end{bmatrix}~,~~~~
\mathbf{S} \,=\, \begin{bmatrix}
			0 \\[1ex]
			\left(\dfrac{\zeta + 1}{\zeta}\right) \dfrac{D}{n} \cdot n^\prime \, T^\prime \\[2ex]
			\dfrac{c_n T}{\epsilon n^2} \cdot n^\prime \,+\, \dfrac{c_T}{\epsilon n} \cdot T^\prime \,+\, \dfrac{G(Z)}{\epsilon}
			\end{bmatrix}~.
\end{align}

The diffusivity function $D(\mathcal{E})$ is given in a few forms:
\begin{align}
	D(Z) \,&=\, \dfrac{D_\text{max} + D_\text{min}}{2} + \dfrac{(D_\text{max} - D_\text{min})\tanh(Z)}{2} ~~~~~~ &\text{Zohm} \\
	D(Z^\prime) \,&=\, D_\text{min} \,+\, \dfrac{D_\text{max} - D_\text{min}}{1 + \alpha_\text{sup}\cdot(Z^\prime)^2} ~~~~~~ &\text{Staps} \\
	D(Z, Z^\prime) \,&=\, D_\text{min} + \dfrac{D_\text{max} - D_\text{min}}{1 + a_1\,Z^2 + a_2\,Z \cdot Z^\prime + a_3\left(Z^\prime\right)^2} ~~~~~~ &\text{Flow-Shear}
\end{align}

### Domain Boundary, Initial Conditions, and Boundary Conditions

Paquay's initial conditions for density and temperature, and Stap's initial condition for $Z$:
\begin{align}
	n(x,0) \,&=\, -\dfrac{\Gamma_\infty \lambda_n}{D} \, \left(1 + \frac{x}{\lambda_n}\right)~, \\
	T(x,0) \,&=\, q_\infty \, \dfrac{\gamma - 1}{\Gamma_\infty} \, \left[1 - \frac{\lambda_n}{\zeta \lambda_T + \lambda_n} \, \left(1 + \frac{x}{\lambda_n}\right)^{-\zeta}\right]~, \\
	Z(x,0) \,&=\, Z_S\left[1 - \tanh\left(\dfrac{L\,x - L}{2}\right)\right] \,=\, Z_S\left[1 - \frac{\exp(L\,x - L) - 1}{\exp(L\,x - L) + 1}\right]~.
\end{align}

\begin{align}
	\text{Boundary:}~~~~\delta \Omega \,=\, \left\{x, t \,\in\, \Omega ~|~ x = 0 ~~\text{and}~~ x = L\right\}
\end{align}

For Matlab, the Neumann and Robin boundary conditions can be expressed in the form of
\begin{align}
	&p\left(x, t, \mathbf{v}\right) \,+\, F\left(x, t, \mathbf{v}, \dfrac{\partial\mathbf{v}}{\partial x}\right) \,=\, 0 ~~~ \text{for} ~~~ (x, t) \in \delta\Omega \\
&p(0, t) \,=\, -\begin{bmatrix}
				D \cdot \dfrac{n}{\lambda_n}\\[2ex]
				\dfrac{D}{\zeta} \cdot \dfrac{T}{\lambda_T} \\[2ex]
				\dfrac{\mu D}{\epsilon} \cdot \dfrac{Z}{\lambda_Z}
				\end{bmatrix}_{x = 0}
~~ \text{and} ~~
p(L, t) \,=\, \begin{bmatrix}
				\Gamma_c(t) \\[1ex]
				\dfrac{(\gamma - 1) q_c - T\Gamma_c(t)}{n} \\[2ex]
				0
				\end{bmatrix}_{x=L}~.
\end{align}

Generalized versions for boundary conditions at the plasma edge ($x=0$):
\begin{align}
	\frac{\partial n}{\partial x} \,=\, \frac{n}{\lambda_n}~, ~~~~\frac{\partial T}{\partial x} \,=\, \frac{n}{\lambda_T}~, ~~~~\frac{\partial Z}{\partial x} \,=\, \frac{Z}{\lambda_Z}~,~~~~ \frac{\partial^2 Z}{\partial x^2} \,=\, 0
\end{align}

...towards the core ($x=L$):
\begin{align}
	\frac{\partial n}{\partial x} \,=\, -\frac{\Gamma_c}{D}~, ~~~~ \frac{\partial T}{\partial x} \,=\ \frac{\zeta\left(T \Gamma_c - (\gamma - 1) q_c\right)}{n\,D}~, ~~~~\frac{\partial Z}{\partial x} \,=\, 0
\end{align}

### Steady-State Solutions

\small
\begin{align}
	&\bar{n}(x) \,=\, \bar{n}(0) - \int_0^x \frac{\Gamma_c}{D(\bar{Z}(x))}~\text{d}x, ~~ \bar{T}(x) \,=\, \frac{(\gamma - 1) q_c}{\Gamma_c} \left(1 - \lambda_g\left(\frac{\bar{n}(x)}{\bar{n}(0)}\right)^{-\zeta}\right), ~~ G(\bar{Z}(x)) \,=\, \theta\,D(\bar{Z}(x)) \\
	&\bar{n}(0) \,=\, -\frac{\Gamma_c \lambda_n}{D(\bar{Z}(0))}, ~~ \lambda_g \,=\, \frac{\frac{\lambda_n}{\zeta \lambda_T}}{1 + \frac{\lambda_n}{\zeta \lambda_T}}, ~~ \theta \,=\, \frac{(\gamma - 1) q_c}{\Gamma_c^2 \lambda_n^2} (c_n + c_g), ~~ c_g \,=\, \frac{\zeta c_T - c_n}{1 + \zeta \frac{\lambda_T}{\lambda_n}}
\end{align}
\normalsize

<!---
Parameters used by Staps:
 $\Gamma_c$  $\gamma$  $\lambda_n$  $\lambda_T$  $\lambda_Z$  $D_\text{min}$  $D_\text{max}$  $\zeta$  $c_n$  $c_T$  $q_c$  $a$  $b$  $c$  $Z_S$  $\alpha_\text{sup}$  $\mu$  $\epsilon$ 
 -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  ----- 
 $-\dfrac{4}{5}$  $\dfrac{5}{3}$  $\dfrac{5}{4}$  $\dfrac{3}{2}$  $\dfrac{5}{4}$  $\dfrac{2}{5}$  $2$  $\dfrac{1}{2}$  $-1.1$  $-0.9$  $-4$  $\dfrac{3}{2}$  $2$  $-1$  $-\dfrac{3}{2}$  $\dfrac{1}{2}$  $\dfrac{1}{20}$  $\dfrac{1}{25}$ 
--->

## Expanded Model for All Plasma Variables
\begin{align}
	\dfrac{\partial n}{\partial t} \,&=\, \dfrac{\partial}{\partial x}\left[D \cdot n^\prime\right] \,+\, 0 \\
	\dfrac{\partial T}{\partial t} \,&=\, \dfrac{\partial}{\partial x}\left[\dfrac{D}{\zeta} \cdot T^\prime\right] \,+\, \left(\dfrac{\zeta + 1}{\zeta}\right) \dfrac{D}{n} \cdot n^\prime \, T^\prime \\
	\epsilon \dfrac{\partial Z}{\partial t} \,&=\, \dfrac{\partial}{\partial x}\left[\mu\,D \cdot Z^\prime\right] \,+\, \dfrac{c_n T}{n^2} \cdot n^\prime \,+\, \dfrac{c_T}{n} \cdot T^\prime \,+\, G(Z)
\end{align}

## Alternate $Z$ Model, with Substitutions
\begin{align}
	&\dfrac{m_i}{e \rho_{pi}} \,n T\, \left(\dfrac{B_\theta}{B}\right)^2 \dfrac{\partial Z}{\partial t} \,=\, \dfrac{m_i \mu_i}{e \rho_{pi}} \,n T\, \dfrac{\partial^2 Z}{\partial x^2} \\
	&+\, B_\theta^2 \left[\left(g_n^\text{an} - g_n^\text{cx} - g_n^{\pi\parallel}\right) \dfrac{n^\prime}{n} + \left(g_T^\text{an} - g_T^\text{cx} - g_T^{\pi\parallel}\right) \dfrac{T^\prime}{T} + \left(g_Z^\text{an} - g_Z^\text{cx} - g_Z^{\pi\parallel}\right) Z - f^\text{OL}\right]
\end{align}

+ Electron Anomalous Diffusion

\begin{align}
	D^\text{an} \,&=\, \dfrac{\epsilon^2 \sqrt{\pi}}{2 a_m} \dfrac{\rho_{pe} T}{B}~,~~~~ g_n^\text{an} \,=\, -e \,n\, D^\text{an}~,~~~~ g_T^\text{an} \,&=\, -e \,n\, \alpha^\text{an}\, D^\text{an}~,~~~~ g_Z^\text{an} \,=\, \dfrac{-e \,n\, D^\text{an}}{\rho_{pi}}
\end{align}

+ Charge Exchange Friction

\begin{align}
	g_n^\text{cx} \,=\, -\dfrac{m_i \,n_0 \langle\sigma v\rangle_\text{cx} \,n T}{B_\theta^2}~,~~~~ g_T^\text{cx} \,=\, \alpha^\text{cx}\,g_n^\text{cx}~,~~~~ g_Z^\text{cx} \,=\, -\dfrac{g_n^\text{cx}}{\rho_{pi}}
\end{align}

+ Ion Bulk Viscosity

\begin{align}
	\text{To consolidate:}&~~~ N \,=\, \dfrac{\nu_{*i}\,\epsilon^{3/2}\,\nu_{ei}}{\nu_{ii}} ~~~\text{and}~~~ \eta \,=\, \dfrac{\epsilon^2 \sqrt{\pi}}{8 a_m} m_i \,n\, (v_{T_i})^2 \\
	\begin{pmatrix}\xi_\theta \\[1ex] \xi_\phi \end{pmatrix} \,&=\, \dfrac{1}{\pi} \int_0^{\sqrt{\nu_{*i}}} \begin{pmatrix} 1 \\ \frac{5}{2} - x \end{pmatrix} x^2 \exp(-x) \, \tan^{-1}\left(\dfrac{2 N \sqrt{x}}{N^2 + Z^2 - x}\right) \text{d}x \\
	g_n^{\pi\parallel} \,=\, \eta \, \rho_{\pi}& B_\theta \, \xi_\theta~,~~~~ g_n^{\pi\parallel} \,=\, \eta \, \rho_{pi} \left(B_\theta\,\xi_\theta - B\,\xi_\phi\right)~,~~~~ g_Z^{\pi\parallel} \,=\, 2\eta \, B_\theta \, \xi_\theta
\end{align}
<!--- Original line for Ion Bulk Viscosity integrals
	\begin{pmatrix}\xi_\theta \\[1ex] \xi_\phi \end{pmatrix} \,&=\, \dfrac{1}{\pi} \int_0^{\sqrt{\nu_{*i}}} \begin{pmatrix} 1 \\ \frac{5}{2} - x \end{pmatrix} x^2 \exp(-x) \left[\int_{-1}^{+1} \dfrac{N / \sqrt{x} ~~ \text{d}y}{\left(y + Z / \sqrt{x}\right)^2 + \left(N / \sqrt{x}\right)^2}\right] \text{d}x \\
--->

+ Ion Orbit Loss:

\begin{align}
	g^\text{OL} \,=\, e \,n\, \nu_\text{eff} \sqrt{\epsilon} \,\rho_{pi}~,~~~~ f^\text{OL} \,=\, \dfrac{g^\text{OL}\,\exp\left[-\sqrt{\nu_{*i} + Z^4}\right]}{\sqrt{\nu_{*i} + Z^4}}
\end{align}

-----------------------------------------------------------

## Extra Information
\begin{align}
	n\,T \,=\, \left(\gamma - 1\right)U, ~~~~ \chi \,=\, \frac{D}{\zeta(\gamma - 1)}
\end{align}

