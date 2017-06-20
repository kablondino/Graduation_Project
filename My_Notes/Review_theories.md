# Notes on Electric Field in _A Review of Theories of the L-H Transition_ by J.W. Connor and H.R. Wilson
## 1. Introduction

+ ASDEX found an empirical law for the heating required for H-mode (with some nuance):
	$$P_{h} \,=\, 0.04 \, \bar{n}_e \, B \, S$$

	+ Above, in MW, with $\bar{n_e}$ as the line-averaged electron density in units of $10^{20}$ m$^{-3}$, $B$ as the magnetic field strength, and $S$ as the plasma surface area.

+ L-H transition requires less heating when the magnetic field config. is such that the ion-$\nabla B$ drift direction is towards a single null formed by the divertor geometry, rather than away from it.

+ Wagner [11] assumes, for the above empirical law, that the transition occurs with the ion diamagnetic velocity achieves some critical value:
	$$V_{*i} \,=\, -\frac{1}{n e B} \frac{\text{d}p_i}{\text{d}r} \,=\, V_{crit}$$

+ Since the transition is an edge phenomenon, it's possible an addition atomic physics parameter is involved. This means one cannot scale the plasma temperature.

#### (i) _Local edge parameters_

+ Also, since the transition occurs at the edge, one would expect that the condition for it to occur can be expressed in terms of local parameter one would expect that the condition for it to occur can be expressed in terms of local parameters and length scales of the edge.

	+ Models are usually expressed in terms of such quantities, or equivalent dimensionless parameters such as the normalized Larmor radius $\rho_{*} \,=\, \rho_i \,/\, a$, collisionality, and $\beta$.

	+ In terms of local plasma parameters, hysteresis is less obvious, but the critical temperature responds to the ion-$\nabla B$ drift in the same manner as $P_h$.

#### (ii) _Spatial and temporal characteristics_

+ How the edge region scales is a potential test of theories.

	+ Transition usually occurs at sub-millisecond or microsecond scales, with a slower evolution in the millisecond range.

	+ The temporal sequence of events for testing the theory: the creation of steep electric field gradients, suppression of turbulence, and steepening of density and temperature profiles.

	+ It is tackled in 3 different ways:

		+ Dimensional analysis, with confinement properties and energy balance involving cross-field transport into the SOL.

		+ Stabilization of instabilities to cause transport in L-mode

		+ **!!** The mechanisms producing bifurcations, such as electric fields and flows. SEE REFERENCES!

+ While the stabilization of some linear instability responsible for transport may be part of it, it is unlikely to account for the rapid changes usually observed in the L-H transition since stability criteria usually depend on parameters that change on a transport timescale.

+ The theories are divided into 3 basic groups:

	+ Those that rely on achieving conditions of plasma parameters for the stabilization of particular instabilities or suppression of the associated turbulence. They are then subdivided into the edge of the core and those in the SOL.

	+ Those that rely on the suppression of turbulent transport by sheared radial electric fields or plasma flow, further classified by the physical processes that drive and limit these quantities.

	+ Those that are considered miscellaneous mechanisms, such as effects arising from divertor action, neoclassical theory, and particular transport models that don't fit above; for example, the role of neoclassical drifts.

+ All of this is summarized in two tables in the appendix.

---------------------------------------

## 2. Theoretical Considerations
### 2.1 _A dimensional analysis_

+ In terms of plasma physics, expect to be able to represent the L-H transition condition by $\rho \propto \beta^x \nu^y$

### 2.2 _Edge power balance_

+ A more consistent approach to obtaining the transition condition is to used a local edge power balance to determine radial length scales, relating to the SOL width $\Delta$.

	+ You can relate the outflow of power into the SOL to a local model for the heat flux across the magnetic field $q_\perp$ at the edge:
		$$P \propto a R n q_\perp \sim a R n \chi_\perp T / \Delta$$

	+ The power into the SOL must be transported along the field lines to the divertor plates ($q_\parallel$ is heat flux parallel to the magnetic field):
		$$P \sim \frac{q_{\parallel} \Delta a}{q}$$

	+ Some assumptions are needed for the perp and parallel heat transport in the SOL: give estimates when perp transport is either Bohm or gyro-Bohm and when parallel transport is dominated by either collisional diffusion or by flow at sound speed.

	+ Balance perp and parallel heat flows in the SOL give:
		$$\Delta(SOL) \sim \sqrt{\frac{n T \chi_\perp R q}{q_\parallel}}$$

+ The results of taking each assumption is given by Eqs. (2.24) through (2.27).

### 2.3 _The shear flow paradigm_

+ Sheared radial electric fields reduce turbulent transport by either stabilization of linear modes or through a reduction of turbulence amplitudes, correlation lengths, or a change in phases between the fluctuations responsible for the turbulent transport fluxes.

+ $V^{\prime}_E$ denotes the radial derivative of the equilibrium $\mathbf{E}\times \mathbf{B}$ velocity, $V_E$.

+ A limit on sheared rotation was found: $\omega_s > \omega_t$:

	+ $\omega_s = k_\theta \Delta r V_E^\prime$, with $\Delta r$ as the radial correlation length and $k_\theta$ as a characteristic poloidal wavenumber of the turbulence.

	+ $\omega_t = 4D / (\Delta r)^2$ as the diffusive scattering rate, with $D$ being the turbulent diffusivity.

	+ Turbulent eddies are torn apart by the differing $\mathbf{E}\times\mathbf{B}$ flows at different radial locations within the eddy, shown in Fig. 7.

	+ Further analysis of a model equation showed that turbulence suppression occurs when $|\omega_t| < |\omega_s|$, which is known as the Biglari-Diamond-Terry (BDT) criterion, which can also be written as:
		$$V_E^\prime > (\Delta r k_\theta \tau_c)^{-1}$$

		+ $\tau_c$ is the turbulent decorrelation time, which can be considered on the order of the inverse of the diamagnetic frequency.

+ A more complete expression for $\omega_s$ is shown in Eq. (2.31), where $\Delta \phi$ and $\Delta \eta$ are correlation angles for the turbulence around the torus and along the field line, respectively:
	$$\omega_s^2 \,=\, (\Delta r)^2 \left[\frac{1}{(\Delta \phi)^2} \left(\frac{\partial}{\partial r} \left(\frac{q V_e}{r}\right)\right)^2 + \frac{1}{(\Delta \eta)^2}\left(\frac{\partial}{\partial r} \left(\frac{V_\parallel}{q R}\right)\right)^2\right]$$

	+ After some generalizations, and assuming the correlation 'lengths' and electrostatic potential $\Phi$ are constant on a flux surface, it can be reduced (the proportionality hold when $\Delta r$ and $\Delta \theta$ are fixed):
		$$\omega_s \,=\, \frac{R^2 B_\theta^2}{r B_\phi} \left(\frac{\Delta r}{\Delta \theta}\right) \frac{\partial^2 \Phi}{\partial \psi^2} \,\propto\, R^3$$

	+ When turbulence is isotropic ($\Delta r \sim r\Delta \theta$), it can be written that (Eq. (2.34)):
		$$\omega_s \,=\, \frac{RB_\theta}{B} \frac{\partial}{\partial r}\left(\frac{E_r}{RB_\theta}\right)$$

+ The curvature of $E_r$ does not have a generic effect on turbulence, but can affect specific instabilities.

+ The transport coefficients are taken to be suppressed by the presence of $\omega_s$:
	$$D \,=\, \frac{D_0}{1 + c\left(\frac{\omega_s}{\omega_t}\right)^\gamma}$$

+ Another mechanism for reducing transport is the effect of $V_E$ on the cross-phase of the fluctuations contributing to the transport flux (source [67]).

### 2.4 _Radial electric fields_

+ $V_E^\prime$ is determined by the radial electric field $E_r$, which is then deduced from the radial force balance for any plasma species $j$ (Eq. (2.38)):
	$$E_r \,=\, -\frac{1}{n_j e_j} \frac{\text{d} p_j}{\text{d} r} + V_{\theta j} B_\phi - V_{\phi j} B_\theta$$

	+ $e_j$ is charge, $p_j$ is pressure, $n_j$ is density, and $V_{\theta j}$ and $V_{\phi j}$ are the poloidal and toroidal velocities.

	+ Therefore, changes in $E_r$ can be associated with changes in $V_{\theta j}$, $V_{\phi j}$, or $p_j^\prime$ (the prime is the radial derivative).

	+ The ones that vary in the poloidal velocity correspond to bifurcations in the solutions of the poloidal momentum balance equation, requiring a momentum source/torque and sink.

		+ Sources: ion-orbit loss, non-ambipolar electron loss, Stringer spin-up due to poloidal asymmetries in turbulent transport, Reynolds stress from turbulence, etc.

		+ Sinks: ion parallel neoclassical viscosity, charge exchange on neutrals, etc.

	+ Models involving $V_{\phi j}$ or $p_j^\prime$ depend on sources of particles (neutrals surrounding the plasma, neutral beams), energy (heating), or toroidal momentum (neutral beams) driving $E_r^\prime$ until a transport bifurcation occurs.

+ If one balances $\omega_s$ from Eq. (2.34) using $E_r$ derived from the $p^\prime$ contribution in Eq. (2.38), then the criterion $\omega_s = \gamma_{max}$ reduces to $\rho_{*s} > \rho_{*s,\text{crit}}$, in which $\rho_{*s,\text{crit}}$ depends on the scaling of $\gamma_{max}$.

### 2.5 _Bifurcations_

+ Bifurcations can take place when the equation has multi-valued or non-monotonic solutions as a function of an 'order' variable, such as density or temperature.

	+ Multi-valued solutions cause 'hard' bifurcations. When the solution is merely non-monotonic, it is a 'soft bif., also known as a 'first-order phase transition.' When the solution remains monotonic but undergoes a change in slope, it is called a 'second-order phase transition.'

+ The poloidal torque balance has the potential to lead to a bif. in poloidal flow, and therefore $E_r$.

	+ The following section of the paper covers this in significant detail. This includes that the particle flux has a nonlinear dependence on the radial electric field.

+ The simplest models for the bif. physics are local in radius, but one can introduce effects to determine the radial structure of the bif. layer.

	+ A small discussion follows on determining the bif. layer width $\Delta$.

		+ If the bif. in $E_r$ arises from $V_\theta$, it could be that the poloidal gyroradius and viscosity determine it.

		+ However, it could be set by transport equations for density and temperature, whose gradients drive $E_r$. (READ)

---------------------------------------

## 3. Instabilities and Turbulence at the Plasma Edge

+ In this section, some cases where the reduction in transport fluxes resulting from a self-consistent generation of sheared flows will be noted.

### 3.1 _Edge Region of the Core_
#### 3.1.1 _Ideal MHD Modes_

+ Bishop [77] proposes ideal MHD ballooning mode stability has a role in L-H transitions.

	+ He concludes that stability properties become progressively better in a separatrix configuration in which the X-point position moves into the region of favorable curvature and that a finite edge current allows the possibility of complete stability to both ballooning and interchange modes for surfaces close to the separatrix.

	+ This is due partly because regions of zero local magnetic shear are move to regions of favorable curvature.
	+ Suggests that the L-H transition occurs well below the simple ballooning limit, Eq. (3.1)

+ Edge plasma current can also drive instabilities (DUH). The peeling mode occurs when a resonant surface for an ideal external kink mode lies just outside the plasma surface. In cylindrical geometry, this mode is destabilized by a finite plasma current at the edge.

	+ These (mentioned in the paper) properties suggest that L-mode is considered to be unstable to the peeling mode, resulting in a high level of magnetic fluctuations at the tokamak edge and therefore an enhanced thermal diffusivity.

	+ As heating is increase, $\alpha$, a parameter, will rise and eventually pass beyond the critical value required for stability; the peeling mode will then switch off and the edge confinement will improve, which can be interpreted as necessary for H-mode.

	+ Whether the peeling mode is stabilized as power increases depends on the competition between the rise in $\alpha$ and the rise in temperature, providing a possible explanation in L-H transitions between large and small tokamaks.

	+ (Look into context) (iv) H-mode can be accessed by reducing the edge current by current ramp-down and by pellet injection (decreases edge temperature, which increases resistivity).

#### 3.1.2 _Resistive Ballooning Modes_

+ _Skipped, for now_

#### 3.1.3 _Tearing Modes_

+ _Skipped, for now_

#### 3.1.4 _Drift Waves_

+ _Skipped, for now_

### 3.2 _Scrape-off Layer Region_
#### 3.2.1 _Resistive Interchagees_

+ _Skipped, for now_

#### 3.2.2 _Electron Temperature Gradient Modes_

+ _Skipped, for now_

#### 3.2.3 _Drift Waves_

+ _Skipped, for now_

---------------------------------------

## 4. Models Involving Sheared Radial Electric Field and Flows

+ It is possible to classify the models by sheared radial electric fields according to the processes generating it: poloidal flow drives and damping, toroidal flows, and transport processes.

+ In addition, one can separate the models according to whether or not they treat the evolution of gradients and fluctuations on the same footing.

### 4.1 _Poloidal Torque Balance_

+ Balance of torques driven by radial currents against poloidal viscosity or damping.

	+ The balance of these torques can lead to multiple solutions for the poloidal flow (or radial electric field); can be associated with low and high transport, _i.e._ L- and H-modes.

#### 4.1.1 _Ion-Orbit Loss Theories_

+ Earliest model that predicts a bif. in $E_r$, due to Itoh and Itoh [69].

	+ They used ion-orbit loss at the edge and balanced it against a non-ambipolar electron diffusion to derive an expression for the equilibrium radial electric field.

	+ The non-ambipolarity comes from the emission of wave momentum which is not absorbed in the plasma but lost in the SOL.
	$$\lambda_e \,=\, -\frac{T_e}{T_i} \rho_{pi} \left(\frac{n_e^\prime}{n_e} + c_e \frac{T_e^\prime}{T_e}\right)$$

	+ $c_e$ is a constant of O(1) that depends on the nature of the microturbulence

		+ For $\lambda_e < \lambda_c \approx \rho_{pi}^2 \nu_{ii} / (\sqrt{\epsilon} D_e)$, only one solution exists, which corresponds to a low (negative) electric field.

		+ As $\lambda_e$ increases to $\approx \lambda_c$, two new solutions exist with higher $E_r$.

		+ For $\lambda_e$ increasing beyond $\lambda_c$, the plasma is forced to jump to the high $E_r$ solution. This is considered equivalent to $\nu_{*i}$ being below a critical value.
			$$\lambda_e = \lambda_c ~~~~\longrightarrow~~~~ \frac{D_e \tau \sqrt{\epsilon}}{\nu_{ii}\rho_{pi} L_n} \,\gtrapprox\, 1$$

		+ When the previous threshold is exceeded, the resulting field is positive, at odds with experimental measurements.

	+ Usually, $\lambda_c$ is of order unity, giving $\rho_{pi} / L_n \gtrapprox 1$, where $L_n$ is the density gradient length.

+ A later model [120] found that for a negative radial field $E_r$, its radial derivative $E_r^\prime$ flips from being small and positive in L-mode to large and negative in H-mode. This jump is hypothesized to be the confinement improvement and associated reduction in turbulence level at the time of transition. This model also has hysteresis.

+ <span style="color:red">Another model has the generation of toroidal flow by ion-orbit loss with similar results, with an extension to a situation where additional heating generates energetic electrons, which are then lost due to magnetic field ripple trapping [121,122].</span>

	+ A bifurcation condition arises, in which $h$ signifies 'hot' electrons, \delta$ is the ripple amplitude, and $V$ is plasma volume:
		$$\frac{n_h}{n_i}\left(\frac{T_h}{T_i}\right)^{7/2} \,>\, 0.17 \frac{\epsilon^3}{\lambda_h \delta^{3/2}} \sqrt{\frac{m_i}{m_e}} \nu_{*i}^2, ~~~~~~ \lambda_h \,=\, \frac{\rho_{pi} n_h^\prime}{n_h}$$

	+ This give a power threshold:
		$$P_{th} \,=\, 10\epsilon^2 \frac{T_h \rho_{pi}}{T_i a} \nu_i n_i T_i V$$

+ The current diffusive ballooning mode (CDBM) model has been examined [123,124], which is based on a resistive ballooning mode turbulence, but where dissipation arises from renormalized viscosity and current diffusion coefficients associated with electron inertia (supported by simulations of turbulence) [125,126].

	+ This includes a reduction in the transport for larger $E_r^\prime$, independent of sign. Using the results from this, the threshold for the L-H transition becomes [127]:
		$$\sqrt{\frac{a}{m_i}} \, \frac{m_e}{e^2 B^2} \, \frac{a_0 q^2}{f(s)} \, \frac{|\text{d}T / \text{d}r|^{5/2}}{\rho_{pi} \nu_{ii} T} \,=\, 1$$

		+ $a_0$ is a numerical coefficient and $f(s)$ is a function of magnetic shear $s$, depending on the mode structure.

		+ When this model is incorporated into a transport code which models the self-consistent evolution of $p^\prime$ and $E_r^\prime$, it can produce a stable edge transport barriere [128].

+ Model of Shaing and Crume also appealed to ion-orbit loss; they balanced the change in poloidal momentum, resulting from loss of the faster ions in the tail of a Maxwellian velocity dist., against the neoclassical poloidal viscous damping ('magnetic pumping'), which is dominated by slower thermal ions.

	+ The result give solutions for the poloidal flow velocity which bifurcates as the temperature is increase: a low velocity (L-mode) and high velocity (H-mode) solution.

	+ The improvement in confinement is assumed to result from the effect of the sheared flow on the correlation length of the turbulence.

	+ The threshold conditions again can be expressed as falling below a critical ion collisionality $\nu_{*i}$, where $c$ is a constant:
		$$\frac{n R q}{\epsilon^{3/2} T_i^2} \,<\, c$$

	+ Predictions of this theory include a poloidal spin-up of the plasma, hysteresis, and an associated negative radial field within a region at the plasma edge of width $\sqrt{\epsilon} \rho_{pi} / \sqrt{s}$._

		+ In addition, the original theory was based on the assumption that the ions take up a Maxwellian velocity dist., that predicts high collisionality plasmas should not exhibit H-mode. *However*, the presence of a fast ion source invalidates that assumption.

+ Shaing and Zhang [130] modified the previous to include stress arising from the interaction of a magnetic field perturbation $b_r$ with the vacuum vessel in the equation for the ion poloidal velocity. Relative to the orbit-loss term, they find an effective stress which can remove the possibility of a bif. in poloidal flow leading to sheared radial electric fields:
	$$\pi_{eff} \,=\, \frac{C_A^2 a \omega \tau_W (r_s / r_W)^{2m} (b_r / B)^2}{V_{Th i}^2 w (1 + (\omega\tau_W)^2 [1 - (r_s / r_W)^{2m}]^2)}$$

	+ $w$ is the magnetic island width arising from $b_r$, $\omega$ is the frequency in the lab frame, $\tau_W$ is the wall time constant, $r_s$ is the resonant surface for the islands, and $r_W$ is the wall radius.

	+ This has the consequence that increasing $b_r$ increases $P_{Th}$ for the L-H transition; conversely, for a given $b_r$, increasing $n$ or $T_i$ assists the transition. Also, the bifurcated flow is smaller.

	+ The aforementioned effects are significant when $\pi_{eff} \geq 1$, which can occur typically for $b_r / B \sim 4\times 10^{-4}$.

+ Hinton and Kim [131] have demonstrated how to link $E_r$ from the SOL to the core across the separatrix, where a radial current $j_r$ must flow, which allows divertor biasing to affect the core.

+ Miyamoto [132] has considered ion-orbit loss in a divertor/separatrix geometry and demonstrated that bifs. in electrostatic potential are possible.

+ Complementarily, Krasheninnikov and Yushmanov [133] show that a potential step can increase ion-orbit loss, allowing an unstable growth of potential and a L-H transition.

+ Kim et al [134] and Hinton et al [135] noted that orbit squeezing, due to electric field gradients comparable with the ion poloidal Larmor radius, can increase ion poloidal flows, thus modifying the effect of poloidal flow damping.

	+ In addition, the same effect reduces the bootstrap current so that the radial field shear can inhibit the bootstrap current drive for MHD instabilities at the plasma edge.

+ Change [136] has examined ion-orbit loss at the X-point, where it is most effective since the poloidal field vanishes there. A critical collisionality depending on the details of the X-point geometry and condition is found, and a layer width $\Delta \sim \bar{\rho_{pi}}$, where $\bar{\rho_{pi}}$ is an average over the flux surface.

#### 4.1.2 _The Effects of Neutral Particles_



#### 4.1.3 _Stringer Spin-up_



#### 4.1.4 _Neoclassical Flows_



#### 4.1.5 _Turbulent Reynolds Stress and Viscosity_



### 4.2 _Transport Bifurcation Theories_



### 4.3 _Phase-Transition Models_



---------------------------------------

## 5. Other Theories for the L-H Transition



---------------------------------------

## 6. Discussion



---------------------------------------

## 7. Conclusions



---------------------------------------

## Tables!
The following is the local L-H transition criteria in terms of dimensionless variables, and also in the form $T_{crit} \,=\, C n^{\alpha_n} B^{\alpha_B}$.

| Theory | Transition Criterion | Critical Temperature |
|-------------------- |-------------------- |--------------------- |
| Peeling mode stability by Connor et al. | $$\frac{\beta q^2}{\epsilon_p} > \alpha_c(s))$$ | $$C \propto \frac{\tau}{1 + \tau} \frac{\epsilon_p}{q^2}$$ |
| Resistive ballooning by Drake et al. | $$\hat{\nu}_{*e} < \sqrt{\frac{m_i}{m_e} \frac{\epsilon_p}{\epsilon_n^2}} \frac{1}{2\pi^2 q_a} \frac{\tau}{1 + \tau}$$ | $$C = 0.019 q \sqrt{\frac{1 + \tau}{\tau} R \epsilon_n} \left[\frac{1}{\epsilon_p A_i}\right]^{1/4}, ~~ \alpha_n = \frac{1}{2}, ~~ \alpha_B = 0$$ |
| Drift resistive ballooning, velocity shear, Stringer spin-up, by Rogers et al. | $$\frac{\epsilon}{q s \sqrt{\epsilon_p}} > 1$$ | Geometric criterion |
| Drift resistive ballooning, velocity shear, $\mathbf{E}\times\mathbf{B}$ flow by Rogers et al. | $$\rho_{*s} > \sqrt{\tau} \frac{s \epsilon_p^{3/2}}{\epsilon}$$ | $$C = 9.58\times 10^4 \frac{s^2 \epsilon_p^3 R^2 \tau}{A_i}, \\ \alpha_n = 0, ~~ \alpha_B = 2$$ |
| Ideal and drift resistive ballooning, ion temperature gradient effects by Rogers et al. | $$\beta > 0.5 \frac{\epsilon_p}{q^2}, ~~~~ \hat{\nu}_{*} < \frac{2\sqrt{2}}{9} \sqrt{\frac{m_i}{\epsilon_n m_e}} \frac{1}{q}$$ | $$n_{cr} = 629 \left[\frac{\tau}{1 + \tau}\right]^{2/3} \frac{\epsilon_p^{2/3} B^{4/3}}{R^{1/3} q^2} \left(\frac{A_i}{\epsilon_n}\right)^{1/6} \\ n > n_{cr}: C = 0.008 q \sqrt{R} \left(\right)^{1/4}, ~~ \alpha_n = 1/2, ~~ \alpha_B = 0 \\ n < n_{cr}: C = 124 \frac{\tau}{1 + \tau} \frac{\epsilon_p}{q^2}, ~~ \alpha_n = -1, ~~ \alpha_B = 2$$ |
| Resistive ballooning, Reynolds stress, Stringer spin-up by Guzdar and Hassam | $$\frac{\hat{\nu}_{*i}}{\rho_{s}(s)\beta} < \sqrt{\frac{\tau}{2}} \frac{D}{D_B} \frac{q^3}{\epsilon_n^2 \epsilon}$$ | $$ \text{Bohm diffusion:} ~~~C = 3.28 \frac{\tau^{5/7}}{(1 + \tau)^{2/7}} \left(\frac{\epsilon \epsilon_n R}{q}\right)^{4/7} \frac{1}{A_i^{1/7}}, ~~ \alpha_n = 0, ~~ \alpha_B = 6/7 \\ \text{Gyro-Bohm diffusion:} ~~~ C = 11.9 \frac{\tau^{5/8}}{(1 + \tau)^{1/4}} R^{3/4} \sqrt{\frac{\epsilon_n}{q}} \frac{\epsilon^{3/4}}{A_i^{1/4}}$$ |
| Drift resistive ballooning, 3D nonlinear simulations by Zeiler et al. | See reference [86], Drake et al. |  |
| Microtearing modes by Ohyabu | $$\hat{\nu}_{*e} < \sqrt{\frac{b}{2} \frac{m_e}{m_i}} \frac{q}{\epsilon_T}$$ | $$\text{For } b = \text{const.:} ~~~ C = 0.225 \sqrt{\epsilon_T R} \left(\frac{A_i}{b}\right)^{1/4}, ~~ \alpha_n = 1/2, ~~ \alpha_B = 0 \\ \text{For } m = \text{const.:} ~~~ C = 3.01 \left[\frac{\epsilon_T \epsilon R^2}{m}\right]^{2/5}, ~~ \alpha_n = 2/5, ~~ \alpha_B = 2/5$$ |
| Low $m$ tearing modes stabilized by drift effects by Strauss | $$\hat{\nu}_{*e} < 0.49 \sqrt{\frac{m_i}{m_e}} \frac{\beta \tau}{1 + \tau} \frac{s\epsilon}{\epsilon_T(r\Delta^\prime)}$$ | $$C = 0.21 \left(\frac{r\Delta^\prime}{s\epsilon}\right)^{1/3} \frac{(R\epsilon_n q)^{1/3}}{A_i^{1/6}}, ~~ \alpha_n = 0, \alpha_B = 6/7$$ |
| Low $m$ tearing modes stabilized by sound wave coupling by Strauss | $$\hat{\nu}_{e*} < 0.39 \frac{q}{s\epsilon_n} \frac{k_y a}{\rho_s (a\Delta^\prime)^2} \sqrt{\frac{m_i}{m_e}} \frac{Z\beta^2 \tau^{7/2}}{(1 + \tau)^2}$$ | $$C = 0.27 \left[\frac{s\epsilon_n}{\epsilon Z}\right]^{2/7} \left[\frac{(a\Delta^{\prime})^2}{k_y a}\right]^{2/7} \tau^{-3/7}, ~~ \alpha = -2/7, ~~ \alpha = 6/7$$ |
| Finite $\beta$ effects on collisional drift wave turbulence by Scott et al. | $$\beta > \frac{m_e}{m_i}$$ | $$C 0.135 \frac{\tau}{(1 + \tau)} \frac{1}{A_i}, ~~ \alpha_n = -1, ~~ \alpha_B = 2$$ |
| Drift-Alfvén turbulence, $k_\parallel \sim 1 / (R q)$, $L_p \sim \Delta_{SOL}$ | $$\sqrt{\frac{m_i}{2 m_e}} \frac{\beta q}{\epsilon_p} > 1 + \left[\hat{\nu}_{*e}^2 \sqrt{\frac{2 m_i}{m_e} \frac{\epsilon_p}{q}}\right]^{1/3}$$ | $$\text{Collisionless:} ~~~ C \sim \frac{1}{(R q)^{4/5}} \frac{1}{A_i^{3/5}}, ~~ \alpha_n = -7/5, ~~ \alpha_B = 2 \\ \text{Collisional:} ~~~ C \sim \frac{(R q)^{2/15}}{A_i^{1/15}}, ~~ \alpha_n = 1 /15, ~~ \alpha_B = 2/5$$ |
| Drift wave instabilities driven by neutrals by Rogister | $$\rho_{*s} \propto \sqrt{\frac{\ell}{a \hat{\nu}_{e}}}$$ | $$C \propto \frac{A_i q}{\epsilon}, ~~ \alpha_n = 2, ~~ \alpha_B = -2$$ |
| Resistive interchange SOL instability, ion-orbit-loss viscocity by Pogutse et al. | $$\rho_{*s} > 10 \sqrt{\frac{m_i}{m_e}} \frac{\Delta_{SOL}^2}{R^2 q \epsilon^{3/2}} \tau$$ | $$C = 1.76\times 10^{10} \frac{\tau^2 \Delta_{SOL}^4}{R^2 \epsilon q^2}, ~~ \alpha_n 0, ~~ \alpha_B = 2$$ |
| Resistive interchange SOL instability, sheath resistivity dominates by Pogutse et al. | $$\rho_{*s} > \frac{31.6}{\epsilon} \left(\frac{\Delta_{SOL}}{R}\right)^{3/2}$$ | $$C = 9.58\times 10^7 \frac{\Delta_{SOL}^3}{R A_i}, ~~ \alpha_n =0, ~~ \alpha_B = 2$$ |
| Resistive MHD, drift/interchange SOL turbulence by Cordey et al. | $$\rho_{*s} > C_1 \frac{\Delta_{SOL}}{\epsilon R}$$ | $$C \propto \frac{\Delta_{SOL}^2}{a_i}, ~~ \alpha_n = 0, \alpha_B = 2$$ |
| Electron temperature gradient modes in SOL, incorporating self-consistent calculation of $\epsilon T_e$, by Cohen an Xu| $$\beta > \beta_c \propto \rho_{*s}^{2/5} \left(\frac{a}{L_\parallel}\right)^{2/5} \frac{\tau^{2/15}}{A_i^{1/10}}$$ | $$C \propto \frac{\tau^{17/12}}{(1 + \tau)^{5/4}} \sqrt{\frac{1}{L_\parallel}} A_i^{1/8}, ~~ \alpha_n = -4/5, ~~ \alpha_B = 2$$ |
| Resistive skin effect reduces turbulence length scales in SOL, by Chankin | $$\beta > \beta_c \propto \left(\frac{\epsilon}{q}\right)^{2/3} \frac{\hat{\nu_{*e}}}{\sqrt{A_i}} \rho_{*s}^{2/3} \sqrt{\frac{1 + \tau}{\tau}}$$ | $$C \propto \left(\frac{\tau}{1 + \tau}\right)^{3/16} \frac{(R q)^{1/8}}{A_i}^{1/16}, ~~ \alpha_n = 0, ~~ \alpha_B = 1/2$$ |
| Ion-orbit-loss torque balancing neoclassical poloidal viscous damping, by Shaing and Crume | $$\nu_{*i} < \nu_c$$ | $$C \propto \frac{\tau}{\epsilon^{3/4}} \sqrt{R q}, ~~ \alpha_n = 1/2, ~~ \alpha_B = 0$$ |
| Ion-orbit-loss, non-ambipolar electron diffusion, by Itoh and Itoh | $$$$ | $$$$ |
| More fundamental condition on ion-orbit loss by Ohkawa and Hinton | $$$$ | $$$$ |
|  | $$$$ | $$$$ |
|  | $$$$ | $$$$ |
|  | $$$$ | $$$$ |


