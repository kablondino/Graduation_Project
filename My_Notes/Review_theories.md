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

		+ $\tau_c$ is the turbulent decorrelation time, which can be considered on the order of the inverse of the diamagnetic freq.

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

	+ $w$ is the magnetic island width arising from $b_r$, $\omega$ is the freq in the lab frame, $\tau_W$ is the wall time constant, $r_s$ is the resonant surface for the islands, and $r_W$ is the wall radius.

	+ This has the consequence that increasing $b_r$ increases $P_{Th}$ for the L-H transition; conversely, for a given $b_r$, increasing $n$ or $T_i$ assists the transition. Also, the bifurcated flow is smaller.

	+ The aforementioned effects are significant when $\pi_{eff} \geq 1$, which can occur typically for $b_r / B \sim 4\times 10^{-4}$.

+ Hinton and Kim [131] have demonstrated how to link $E_r$ from the SOL to the core across the separatrix, where a radial current $j_r$ must flow, which allows divertor biasing to affect the core.

+ Miyamoto [132] has considered ion-orbit loss in a divertor/separatrix geometry and demonstrated that bifs. in electrostatic potential are possible.

+ Complementarily, Krasheninnikov and Yushmanov [133] show that a potential step can increase ion-orbit loss, allowing an unstable growth of potential and a L-H transition.

+ Kim et al [134] and Hinton et al [135] noted that orbit squeezing, due to electric field gradients comparable with the ion poloidal Larmor radius, can increase ion poloidal flows, thus modifying the effect of poloidal flow damping.

	+ In addition, the same effect reduces the bootstrap current so that the radial field shear can inhibit the bootstrap current drive for MHD instabilities at the plasma edge.

+ Change [136] has examined ion-orbit loss at the X-point, where it is most effective since the poloidal field vanishes there. A critical collisionality depending on the details of the X-point geometry and condition is found, and a layer width $\Delta \sim \bar{\rho}_{pi}$, where $\bar{\rho}_{pi}$ is an average over the flux surface.

#### 4.1.2 _The Effects of Neutral Particles_

+ Neutrals can influence ion-orbit-loss theories in 2 ways:

	+ Increasing both poloidal flow damping and ion-orbit losses through charge exchange with ions.??

+ Itoh and Itoh [121] introduced poloidal flow damping due to momentum losses from charge exchange with neutrals. This modifies the bifurcation condition:
	$$\frac{\lambda_e}{\lambda_c} \,=\, 1 + d_n \lambda_i, ~~~~~~~~ \lambda_i \,=\, -\rho_{pi}\left(\frac{n_i^\prime}{n_i} + \frac{c_i T_i^\prime}{T_i}\right) \\
	c_i \sim 1, ~~~~~~~ d_n \,=\, \sqrt{\epsilon} \, \frac{n_o \langle\sigma_{cx} v \rangle}{\nu_{ii}F} \frac{\ell_n}{\rho_{pi}}$$

	+ $n_o$ is the neutral density, $\sigma_{cx}$ is the charge exchange cross-section, $\ell_n$ is the neutral penetration length, and $F$ is a form factor $\sim 1$. This effect increases the value of $\lambda_c$ needed for a bifurcation:
		$$n_o \,>\, \frac{c\nu_{ii}}{\sqrt{\epsilon} \langle \sigma_{cx} v \rangle} \frac{\rho_{pi}}{\ell_n}$$

	+ These two conditions provide limits on both electron and neutral particle densities (Eqs. 4.14 and 4.15)

	+ Neutrals were also found to increase the time for a transition to take place. The paper also considered the effects of an impurity ion.

+ The same authors have also related the neutral population to the nature of the wall materials, showing that the condition for H-mode can be written as:
	$$P(v_f)r_f \,<\, \left[\frac{\langle \sigma_{ion} v \rangle n_e + \langle \sigma_{cx} v \rangle n_i}{\langle \sigma_{recomb} v \rangle n_e + (\nu_i \rho_{pi} / \ell_{recomb})}\right] \frac{\nu_i}{\langle \sigma_{cx} v \rangle n_i} \frac{\rho_{pi}}{\ell_{recomb}}$$

	+ $r_f$ is the reflection coefficient of fast neutrals at the wall; $P(v_f)$ is the probability that the neutral has a velocity $v_f$; $\sigma_{ion}$, $\sigma_{cx}$, and $\sigma_{recomb}$ are the cross-sections for ionization, charge exchange, and recombination of neutrals, respectively; and $\ell_{recomb}$ is the mean-free path for recombination.

	+ Since $r_f$ is higher for wall materials with higher $Z$, this could explain why these mitigate against H-mode.

+ In [138], neutrals were investigated deeper at the X-point. It was found to lead to an increase in the effective collision freq for ion-orbit loss.

	+ For a set of parameters, there is a critical value for the neutral density near the X-point for a transition: thus 'condensed' neutrals can trigger H-mode, which is consistent with experiment. It can be written:
		$$\xi \,<\, \frac{d_o \left(\frac{1}{2} + \sqrt{\nu_{*i}}\right)^3 \exp\left(\frac{1}{2} + \sqrt{\nu_{*i}}\right)}{2 \left(\frac{1}{4} + \sqrt{\nu_{*i}}\right)^{3/4} \left(\frac{3}{2} + \sqrt{\nu_{*i}}\right)} \,\equiv\, \xi_c, ~~~~~~~~~~ d_o \,=\, n_o^{main}\langle \sigma_{cx} v \rangle \frac{1 + 2q^2}{\sqrt{\epsilon} \omega_{bi} q^2}$$

		+ $\xi = \nu_{in_o}(X) / \omega_{bi}$, with $\nu_{in_o}(X)$ as the ion-neutral collision freq at the X-point and $\omega_{bi}$ as the ion bounce freq.

+ Shaing and Hsu [139] considered charge exchange losses of poloidal momentum in the ion-orbit-loss model. They found the critical parameter is:
	$$\hat{\nu} \,\equiv\, \frac{n_o \langle \sigma_{cx} v \rangle R q}{V_{Th,i}}$$

+ Xiao et al [140] investigated the effect of charge exchange on ion orbits. They found that a negative value of $E_r$ produces an inward ion mobility so that a positive bias helps H-mode.

	+ The opposite bias is less helpful, relying on detrapping and compressing of passing particle orbits.

#### 4.1.3 _Stringer Spin-up_

+ Pfirsch-Schlüter transport can drive a poloidal spin-up. This mechanism is damped by the magnetic pumping which exceeds the drive, so that no spin-up is predicted by neoclassical theory.

	+ However, Hassam et al [70] has proposed, in reaction to higher observed transport, a 'spontaneous poloidal spin-up' associated with anomalous transport which is not constant on a flux surface. Balancing this drive against the damping from magnetic pumping gives a threshold:
	$$\frac{1}{n r} \frac{\partial}{\partial r}(n r \widetilde{v}_r) + \epsilon \frac{1 + 2q^2}{q^2} \gamma_{MP} \,<\, 0$$

	+ $\gamma_{MP}$ represents the magnetic pumping and $\widetilde{v}_r = v_r - \langle v_r \rangle$, with $v_r$ as the radial particle flow velocity and the angled bracketed one as a flux surface average.

	+ A simplified condition can be deduced for the onset of the poloidal flow instability, by $v_r$ represented by a diffusive term and a pinch term that is constant on a flux surface. Also, if we assume $\gamma_{MP} \propto \nu_{ii}$, then:
		$$\frac{\epsilon}{q^2} (1 + 2q^2) \nu_{ii} \,\lessapprox\, \delta D \left(\frac{n^{\prime\prime}}{n} + \frac{n^\prime}{n r}\right)$$

	+ Hassam and Antonsen [142] have extended this to include toroidally directed, poloidally asymmetric momentum sources, but did not include diamagnetic flows or gyro-viscosity.

+ McCarthy et al [143] did an analytic and numerical investigation; they noted that if turbulent velocities are comparable to diamagnetic ones, then the Stringer spin-up effect is greater than the Reynolds stress by a factor of $L_n^3 / R \rho_s^2$.

	+ Their model was toroidally axisymmetric but included a primitive SOL with a particle sink; a flux consisting of a poloidally asymmetric diffusive part and a pinch, and flow damping from magnetic pumping were also included.

	+ Strong magnetic pumping is taken to correspond with L-mode, when it's found that poloidal flow balances the diamagnetic velocity, as in neoclassical theory.

	+ For the case of low magnetic pumping, large flows in the opposite direction to the diamagnetic velocity develop, with a shear flow scale length ($\lambda_{e,mfp}$ is the electron mean-free path and $\rho_e$ is the electron Larmor radius):
		$$L_v \,\sim\, \left(\frac{a \rho_e}{L_n \lambda_{e, mfp}}\right)^{1/4} L_n$$

		+ The L-H transition would therefore be expected to relate to a change of collisionality.

+ The effects of a divertor X-point geometry were considered by Strauss [144] using a dissipative reduced MHD model with viscous stresses in which poloidal variations are opposed by sound wave propagation. He found that the condition for spin up was:
	$$\frac{\chi_\perp}{L_p^2} \, \widetilde{G} \,>\, \frac{B_\theta}{B} \frac{C_s}{r}$$

	+ $\chi$ is the thermal diffusivity and $\widetilde{G}$ is a form factor representing the poloidal variations of the geometry. For a Bohm-like diffusivity, this corresponds to a critical value of $\rho$.

#### 4.1.4 _Neoclassical Flows_

+ Poloidal viscous forces in neoclassical theory produce a poloidal ion flow proportional to the ion temperature gradient:
	$$V_{\theta,i} \,=\, \frac{k}{e B} \frac{\text{d}T_i}{\text{d}r}$$

	+ $k$ depends on the ion collisionality $\nu_{*i}$. This results in a negative $E_r$, observed experimentally.

	+ Plasma microturbulence is expected to produce an 'anomalous' viscosity and inertia and this might affect the poloidal momentum balance equations.

+ Rozhansky and Tendler [145] have taken these into account to derive an equation for the poloidal flow which balances the neoclassical visc. with the anomalous visc. and inertia, excluding toroidal flow.

	+ The neoclassical terms dominates for low edge temperature gradients, which means the equation has a low poloidal flow solution (L-mode).

	+ For higher edge temperature, the poloidal flow is dominated by the anomalous terms which predicts a much higher flow velocity (H-mode). The transition condition:
		$$\frac{\text{d}^2 T_i}{\text{d}r^2} \,\gtrapprox\, \sqrt{\pi} \,\epsilon q\, \frac{n(a) C_s \tau_p}{\bar{n} a^2 (1 + 2q^2)} \frac{\text{d}T_i}{\text{d}r}$$

		+ $\tau_p$ is the particle confinement time, $n(a)$ is the edge density, and $\bar{n}$ is the average density.

+ The role of sheared radial electric field associated with the flow (4.23) in suppressing particle diffusivity was investigated by Rozhansky [146,147]. In theirs, there is a steep drop in the diffusion coefficient when the shear suppression parameter reaches a critical value. The shear suppression parameter:
	$$\kappa \,=\, \rho_{pi}^2 \frac{\text{d}^2 (\ln n)}{\text{d}r^2}$$

	+ The steep drop in the diffusion coefficient is most easily satisfied at the edge, where $n$ is smallest. The resulting continuity equation is completed by a condition on the flux from the core and a boundary condition at the separatrix, which can represent the SOL.

	+ They found a L-H transition can be triggered by a change in the core density gradient; the flux must exceed a critical value which increase $T_i$ and decreases with $n$.

		+ The transition front propagates faster than any diffusion time, and shows hysteresis.

		+ Introducing time scales $\tau_1$ for core changes due to the formation of the barrier and $\tau_2$ for the response of the SOL, this model can cause dithering.

	+ The sensitivity of $\kappa$ to the density profile means that the transition can be triggered by pellet injection, adiabatic compression, and changes in the SOL.

+ Hinton [149] invoked the velocity (4.23) to develop a bif. model for the transition by retaining the collisionality dependence of $\mu$, the parallel visc., in the Pfirsh-Schlüter regime and keeping only contributions to $V_{\theta,i}^\prime$ from $T_i^\prime$. He used a transport model:
	$$\chi \,=\, \chi_0 + \frac{\chi_1}{1 + \alpha V_{\theta,i}^{\prime \gamma}}$$

	+ $\chi_0$ is residual H-mode thermal diffusivity, $\chi_1$ is L-mode diffusivity. The value of $\alpha$ can be chosen to reflect the BDT criterion, equation (2.30), for, say, drift wave turbulence, and $\gamma$ is some parameter.

	+ This form leads to a transport flux which exhibits a maximum as a function of $T_i^\prime$ when the heat flux reaches a critical value; this allows a bif. to an 'H-mode' state with higher $T_i^\prime$.

	+ Global solutions of the transport equations give an edge region of steep gradients whose radial extent is determined by where the heat flux from the core just exceeds the critical value. The corresponding solutions for the energy confinement time $\tau_E$ display hysteresis.

	+ This model is considered generic, needing more details for $\chi$.

#### 4.1.5 _Turbulent Reynolds Stress and Viscosity_

+ The turbulence itself can drive flows through Reynolds stress [150].

	+ The non-ambipolar electron transport from [69] can be considered (for the alternate view of turbulent torque coming from the MHD dynamo [151]) since it corresponds to a radial current which can produce a torque.

	+ [72] calculated the effects of the turbulent Reynolds stress; they found only turbulence which supports radially propagating waves is able to drive rotation (_e.g._ drift wave turbulence).

		+ However, since flow shear can itself lead to this propagation, this can be amplified for other types of turbulence such as resistive interchanges [152].

	+ The turbulence should be radially inhomogeneous, which tends to limit rotation to the plasma edge where either

		a. density and temp. gradients are sufficiently steep that they have a significant variation across a radial correlation length of turbulence or...

		b. the radial correlation length is cut off by the plasma boundary.

	+ No detailed calculations (as of this review) of the resulting steady-state flow have been done and no bif. condition is evident.

+ Because the turbulence to generate a Reynolds stress is directional, it occurs most readily at the plasma edge (supported by experiments [153]) and exists over a width $\Delta$ related to the radial mode structure of the turbulence.

+ This width is likely to be microscopic (such as $\rho_i$) in the cylindrical approx., the theory of edge ballooning modes show that it could be larger in toroidal geometry.

	+ Particularly, for drift wave or ion temp. gradient turbulence, one expects $\Delta \sim \rho_i^{2/3} a^{1/3}$ [154,155].

		+ Indirect measurements of $\Delta$ on JET, based on the pedestal height and the assumption that the steep edge gradient is at the ideal MHD ballooning limit, are broadly consistent with this. It also offers an explanation of the isotope scaling of the pedestal temperature [20].

+ The $\eta_i$ (ion temp gradient) driven mode is a strong candidate for the instability responsible for the L-mode anomalous transport [159].

	+ Considering a slab model, [160] have demonstrated that a radial electric field gradient has a strong stabilizing influence on the $\eta_i$ mode and might explain the improvement in transport observed in H-mode, although it could not explain the origin of the radial field.

	+ Later, the same authors [161] calculated the turbulent viscosity resulting from the non-linearly saturated $\eta_i$ instability. They obtained the result that this anomalous viscosity increases with radial electric field shear to some critical value beyond which it fails. This allows for the picture of H-mode as follows:

		+ Some drive (_i.e._ momentum source) spins up the plasma in a continuous way to give an equilibrium flow determined by the balance between the drive and the anomalous viscous damping.

		+ At some critical value of the drive, the radial electric field derivative reaches the critical value $E_{r,crit}^\prime$ where the viscous damping is max.:
			$$E_{r,crit}^\prime \,=\, 2.02\times 10^4 B \sqrt{T_e / A_i} ~ \tau^{-0.77} \, L_n^{0.7} \, L_{T_i}^{-0.26} \, L_s^{-1.44} ~ \text{V m}^{-2}$$

			+ $B$ is in T, $T_e$ in eV, and lengths are in m. This is just a prediction.

			+ The authors postulate a smaller 'background' viscous damping which increases linearly with the electric field gradient and does not possess a max.

	+ At the critical point, the plasma jumps to the part of the curve representing background viscosity, which then determines the radial electric field in H-mode. The result corresponds to a flow shear scaling like Eq. (2.28).

	+ The improvement in confinement results from the stabilization of the $\eta_i$ mode with the sudden increase in the electric field gradient. The power threshold corresponds to the point at which the critical electric field is reached, estimated as follows:

		+ Two mechanisms may be considered for the momentum input to balance the turbulent viscosity: ion-orbit loss as [68] and a NBI source.

		+ The balance of the momentum source $M$ and the turbulent viscosity gives the following, with $V$ as the plasma volume, and $H$ as a form factor assumed to be a function of only $E_r^{\prime}$ and represents a variation of the turbulent viscosity with $E_r^\prime$:
			$$M \,\sim\, \frac{n T_e V \rho_s}{L_n^2} \, H$$

		+ Neglecting the poloidal rotation (so that $\mathbf{E}\times\mathbf{B}$ poloidal flow is balanced by the diamagnetic poloidal flow) and assuming $M$ to be given by the ion-orbit-loss mechanism, gives a condition that corresponds to a minimum collisionality:
			$$f \, v_{*i}^{1/2} \, \gtrapprox \, \frac{\rho_s a}{L_n^2}$$

		+ This condition results in a power threshold assuming forms for the scaling of temps and the fraction of fast particles $f$ with power.

		+ For NBI sources, we can assume $M \propto P$ and obtain:
			$$P \, \gtrapprox \, P_{Th} \,\sim\, \frac{n T V \rho_s C_s}{L_n^2}$$

	+ This theory also allows for a smooth transition. If the $\eta_i$ mode is completely stabilized by $E_r^\prime$ before the max in the anomalous viscosity is reached, the no spontaneous spin-up of the plasma occurs and the threshold power corresponds to that power which provides sufficient radial electric field gradient to stabilize the $\eta_i$ mode.

		+ The opposite case has a bifurcation-type H-mode with the spontaneous creation of the radial electric field.

### 4.2 _Transport Bifurcation Theories_

+ In transport bif. theories, the transport suppression from $V_E^\prime$ is determined by $E_r^\prime$ arising from $p_i^\prime$ through radial force balance.

	+ Practically, these are closely related to [149], where $E_r^\prime$ is driven by the neoclassical $V_{\theta,i}$ equation (4.23)._

+ Bifurcation becomes possible when the sources of energy and particles can drive sufficient gradients.

+ When the poloidal and toroidal flows can be neglected, the $\mathbf{E}\times\mathbf{B}$ velocity is given by:
	$$V_E \,=\, = -\frac{1}{n_i e B} \frac{\text{d}p_i}{\text{d}r}$$

	+ This expression is correct if standard neoclassical theory is valid at the edge and when the ion temp gradient is neglected. However, the inclusion of orbit-squeezing effects which result from the existence of a radial gradient in the radial electric field can lead to significant poloidal rotation (see [135]).

	+ In the case of orbit-squeezing, $V_E$ would be suppressed by a factor of $S = 1 - m_i E_r^\prime / (e B_{\theta}^2)$, describing the orbit-squeezing, and leads to a mod of the parameters $\alpha$ and $\gamma$ in Eq. (4.26)._

+ Using the above eq, the particle and energy fluxes can be expressed solely in terms of density and pressure gradients (and the two parameters). Particle and energy conservation then yield expressions for the fluxes in terms of the heat and particle sources [162]. The resulting system can be solved to determine pressure and temp profiles.

	+ For a given particle flux $\Gamma$ and heat flux $q$ less than a critical value, only one solution exists for the ion density and pressure gradient pair, both of which are low (L-mode).

	+ When $q$ is above the critical value, the solution jumps to the high gradient solution, interpreted as H-mode. It is easiest to access H-mode at the plasma edge where $\Gamma$ is largest.

	+ A threshold condition is given that corresponds to the product of fluxes, where $D_0$ and $D_1$ are H- and L-mode particle diffusivities, respectively, and $\chi_0$ and $\chi_1$ are the corresponding thermal diffusivities:
		$$(q\Gamma)_{crit} \,=\, \frac{9}{16} \frac{D_1 \chi_1}{\sqrt{3\alpha}} \left(1 + \frac{4 D_0}{3 D_1}\right) \left(1 + \frac{4 \chi_0}{3 \chi_1}\right), ~~~~~ \text{(normalized)}$$

		+ Applying the condition at the edge means that $q$ becomes the total input power per area and $\Gamma$ can be calculated from the source resulting from the neutrals diffusing in from the edge and becoming ionized inside the plasma.

	+ If one interprets $D_0$ and $\chi_0$ as the neoclassical diffusivities and $D_1$ and $\chi_1$ as anomalous diffusivities, the above eq complicated scaling for a power threshold, even in the limit that the neoclassical diffusivities are neglected.

		+ The model also depends on the model which is chosen for recycling.

		+ As a note, bifurcations are only possible when $\gamma > 1/2$ and $D_1/ D_0 > 2\gamma / (\gamma - 1/2)^2$.

	+ Therefore, a bif. transition is expected where the edge transport is dominated by turbulent processes, and a smooth transition would occur where there is more neoclassical transport.

+ This work has been extended to include poloidal and toroidal flows, as well as density and temp gradients [62,163]. The reduction in drift wave transport is:
	$$D \,=\, D_0 + \frac{D_1}{1 + \alpha |S_\perp|^\gamma}, ~~~~~~~~~ S_\perp \,=\, \frac{B_\phi}{B} V_\theta^\prime + \frac{B_\theta}{B} \frac{\Pi_\theta}{\mu_\perp} + \frac{q\Gamma}{e B \chi_\perp D_\perp}$$

	+ $q$, $\Gamma$, and $\Pi_\theta$ are heat, particle, and toroidal momentum fluxes, with corresponding diffusivities $\chi_\perp$, $D_\perp$, and $\mu_\perp$.

	+ For transitions, if $\Gamma$ is expressed in terms of the edge recycling particle flux, then the recycling coefficient $r < 0.98$.

	+ The width of the edge transport barrier is $\Delta \propto \ell_i \ln(P)$, meaning it is only weakly varying with additional power and proportional to the neutral ionization depth. The barrier is formed in a few milliseconds: $\Delta \tau \sim (\Delta / a)^2 \tau_p$, for $\tau_p$ as particle conf time.

	+ It also addresses the generation of poloidal flows by biased probes or ion-orbit-loss mechanisms. The threshold power reduces when $B_\phi (\text{d}V_\theta / \text{d}r) > 0, which is expected for ion-orbit loss [163].

	+ The condition for orbit-squeezing to produce a large poloidal flow $\rho_{pi} / L_p > 1$ typically exceeds that needed for a transition arising from the normal pressure gradient contributions to $S_\perp$ in (4.34)._ [164]

+ [165] have given a simplified generic model for the transition based on ion temperature gradient turbulence and allowing for rotational shear stabilization.

	+ It is the increase in rotational shear stabilization with $\nabla T$, or $q$, which allows for a bif. at a normalized power <span style="color:red">(4.35)</span>.

+ [166] considered situations with RF electron heating where the preferential loss of hot electrons can generate sheared radial electric fields.

	+ Two cases, one with thermalized and the other with suprathermal electrons, were investigated using transport modelling and both led to transitions, but it had no specific predictions to compare with experiment other than a positive sign of $E_r$.

### 4.3 _Phase-Transition Models_

+ 'Phase-transition' models evolve in time a set of spatially-averaged, coupled equations for the turbulence level, the driving gradients, and the consequent radial electric field, with transport coefficients responding to the fluctuations and the radial electric field (or flow) shear [120].

	+ If turbulence dynamics time scale is removed, they collapse to the transport bif type.

	+ A turbulent dynamo effect amplifies the flow shear until the fluctuations are self-regulated as the flow shear suppression exerts itself.

	+ The threshold power follows from the condition that the energy input rate in the absence of flow shear exceeds neoclassical or charge exchange poloidal flow damping.

+ [167,168] investigated the interaction between flows and turbulence in a generic system of equations. The fluctuation energy $E = |\widetilde{n}_k / n |^2$ satisfies:_
	$$\frac{1}{2} \frac{\text{d}E}{\text{d}t} \,=\, \gamma_0 E - \alpha_1 E^2 - \alpha_2 U E$$

	+ $\gamma_0$ is the linear growth rate, $\alpha_1$ is the nonlinear damping due to coupling to other helicities, and $\alpha_2 U \equiv \alpha_2 V_E^{\prime 2}$ represents the effect of a sheared $\mathbf{E}\times\mathbf{B}$ flow in stabilizing the instability.

	+ The sheared $\mathbf{E}\times\mathbf{B}$ flow itself is governed by the competition between the Reynold stress drive due to the fluctuations and damping due to the parallel viscosity $\mu$ since:
		$$V_E \,=\, -\frac{e}{B} \,=\, V_\theta - \frac{\rho_s C_s}{p_i} \frac{\text{d}p_i}{\text{d}r}$$

		+ This leads to $\dfrac{1}{2} \dfrac{\text{d}U}{\text{d}t} \,=\, -\mu U + \alpha_3 U E$, where the ion pressure gradient contribution is dropped for simplicity.

	+ The $\alpha_i$'s can be estimated for different types of turbulence. This system forms a 'predator-prey' model, in which the fluctuation intensity $E$ is the prey to the shear flow U (predator), and depends only on the parameters $a = \alpha_3 / \alpha_1$ and $b = \mu / \gamma_0$.

	+ There are two fixed points: $E = \gamma_0 / \alpha_1$, $U = 0$ and $E = \mu / \alpha_3$, $U = (\gamma_0 - \alpha_1 \mu / \alpha_3) / \alpha_2$.

		+ The first is stable if $\gamma_0 > \alpha_1 \mu / \alpha_3$ and can be identified as H-mode.

	+ For a generic drift wave model, $\gamma_0 \sim k_\theta \rho_s C_s / L_T$. Power balance at the edge implies $L_T^{-1} \sim P / (n T a R \chi)$, which give the bif condition <span style="color:red">(4.39)</span>.

		+ A scaling for $\alpha_1 can be obtained by assuming a decorrelation freq of order $\omega_{*}$, a radial mode width of order $\rho_s$, and taking $\chi \sim C_s \rho_s^2 / L_n$, giving <span style="color:red">(4.40)</span>. Furthermore [167,168] states that $\alpha_3 \sim C_s / L_s$ for generic drift wave turbulence, giving:
			$$P_{th} \,=\, \mu n T L_s a R$$

	+ A more direct estimate of the bif condition [167,169] leads to a critical value of the local quantity, where $\Delta r \sim \rho_i L_s / L_n \approx 2$ cm:
		$$\lambda \,=\, \frac{V_{Th,i}}{L_s \mu} \left(\frac{\Delta r}{L_n}\right)^2$$

		+ This gives a critical value of $\nu_{*i}$ in the banana region (where $\nu_{*} < 1$ in (4.42))and $\rho_{*}$ in the plateau regime (where $\nu_{*} > 1$) and is more easily satisfied in a divertor config where $L_s$ is smaller.

		+ This model has also been extended to include an equation for the evolution of the ion pressure gradient and its contribution to $V_E$ [170].

		+ This leads to a L-H transition (first-order phase transition) in $V_E^\prime$ and $E_r$ as a function of heat flux $q$. However, there is a 2nd, later transition to zero fluctuations.

	+ The inclusion of $p_i^\prime$ (ion pressure gradient) leads to a preferred direction for a seed radial electric field gradient $E_r^\prime$. Near the transition, the $V_\theta$ contribution to $E_r$ dominates due to $p_i^\prime$ and conversely for $P \gg P_{Th}$; thus $P_{Th}$ is unaffected.

	+ After the transition, the fluctuations are quenched, and therefore, the drive for $V_{\theta}^\prime$ vanishes so that only the $p_i^\prime$ contribution remains.

		+ A rapid power ramp compresses the time duration of the sheared flow generation phase, which may render it unobservable.

		+ On the other hand, the transition time becomes logarithmically singular at the transition point $P \sim P_{th}$, which is like a 2nd-order phase transition, _i.e._ there are jumps in temporal gradients of quantities, not in the quantities themselves.

+ An elaboration of the model [171] is to introduce an external torque $T_{ext}$ to drive $V_\theta$, which can be done by a biased limiter, probe, or RF waves. It results in a reduction of the critical power:
	$$P_{Th} \,=\, P_{Th}^0 - K T_{ext}^{2/3}$$

	+ Another extension [169,172] was made to describe the radial propagation of the L-H transition front by including radial diffusion effects in the equations for fluctuation energy and poloidal velocity shear. The velocity of propagation is give, with $D_1$ as the diffusion coeff and $L_p$ as the pressure gradient length in L-mode:
		$$V_f \,=\, \frac{2D_1}{L_p} \sqrt{\frac{P}{P_{Th}} - 1}$$

	+ With the 3-eq model described, one obtains a propagating form for $P \gtrapprox P_{Th}$, but at $P \gg P_{Th}, there is no propagation and the fluctuations collapse uniformly everywhere.

+ [173] have introduced neutrals, which increase the effective damping of $V_\theta$ by charge exchange, decrease the effective energy flux through charge exchanged and ionization, and provide a particle source. The local particle and heat fluxes act as order parameters._

	+ In a model that neglects $q$ and $V_\theta$, the bif condition is symmetric in particle flux per area from the core $\Gamma_c$ and the neutral influx $\Gamma_o$:
		$$\Gamma_c + \Gamma_o \,>\, \frac{n k_\theta \rho_s C_s}{4}, ~~~~~~~ \text{for drift wave turbulence: } k_\theta \rho_s \leq 1$$

		+ This more-complete model has 4 parameters, one being $g_1 = q / T\Gamma$.

	+ It was found that the threshold particle flux (the above eq) decreases with $q$ in the absence of neutrals but increases in their presence.

	+ If the energy losses in the edge due to ionization and charge exchange compete with the core heat flux ($g_1 < 1$), then increasing the neutral density lowers the threshold because $n^\prime$ dominates $E_r$ and the neutral source increases $n^\prime$.

		+ This case is like the simple model (4.46); it has the implication that H-modes with $L_T > L_n$ are more stable against changes in $\Gamma_o$ inducing a back transition.

		+ For cases $g_1 \sim 1$, the $n^\prime$ and $T^\prime$ contributions to $E_r$ cancel and there is a max in the threshold.

		+ The threshold is raised if $T_i / T_e$ is decreased, so that edge electron heating is unfavorable.

+ The phase-transition model can be elaborated to obtain a position and width for the barrier. The location of the transition can be obtained from the Maxwell construction [75] and in a time-dependent calculation, $x$, the radial position of the barrier at time $t$, propagates as:
	$$x \,=\, \sqrt{\frac{D\,t\,(\Gamma - \Gamma_M)}{\Gamma_M}}$$

	+ $D$ is the diffusion coeff, $\Gamma$ the flux, and $\Gamma_M$ is the Maxwell flux.

	+ In the case of an edge neutral source, [174] deduced that the barrier (pedestal) width is:
		$$\Delta_{ped} \,=\, \sqrt{\frac{D}{\nu_{ion}}} \, F\left(\frac{\Gamma_o}{\Gamma_M}\right)$$

		+ $\nu_{ion}$ is the ionization rate for neutral atoms at the plasma surface and $F$ is a weak function of the incident neutral flux $\Gamma_o$ that provides an ionization source at the plasma edge.

+ The basis for these phase-transition models have been investigated using non-linear simulations of resistive pressure gradient driven turbulence in [175].

	+ It is found that a strong cancellation between turbulent Reynolds stress drive for $V(\theta)$ and damping due to turbulent viscosity means that the transition can only occur at low fluctuation amplitudes near marginal stability. The reduction in transport in H-mode is partly due to the changed phase relationship between fluctuations.

+ The possibility of active control of the transition is done in [176], where it is shown that low-freq poloidal flows are more weakly damped by neoclassical viscosity than are static ones, which could allow the driving stabilization oscillating poloidal flows.

+ Other ideas [177], they develop a dynamical model to describe electrostatic resistive pressure gradient driven turbulence in slab geometry. This consists of 3 ODE's describing the evolution of electrostatic potential, plasma flow, and fluctuation energy.

	+ It is derived from energy balance by modelling the transfer and dissipative terms, with the energy input as an order parameter (rather than having $p_i^\prime$).

	+ The transition is triggerd by an increase in $T_i$, leading to a decrease in the neoclassical poloidal viscosity that opposes Reynolds stress.

		+ For the transition to occur, the edge diffusivity must exceed the viscosity ($T_i \geq 50 - 100$ eV), giving Eq <span style="color:red">(4.49)</span> as a power threshold. The time for transition is usually $\sim 50 \mu$s.

	+ At low densities (_i.e._ banana regime), one finds a 1st-order phase transition showing hysteresis.

	+ For higher densities (plateau), a 2nd-order phase transition occurs. This model was completed and verified by [178].

+ [179] is a phase-transition model called the edge turbulent layer (ETL) model.

	+ A transition from zero plasma flow to a finite flow with reduced turbulence occurs a critical $\nabla T$, by solving a Lorenz-like system:

		+ It is done by solving a Lorenz-like system of equations for velocity shear, turbulent kinetic energy, thermal energy, and the temperature profile, with the incident heat flux and temperature $T_b$ (temperature gradient across the turbulent layer) at the core-edge interface.

	+ The transition can only occur if the poloidal viscosity varies inversely with temperature. This requires $T > \sqrt{0.7 n q R} \, Z(\text{eff})^2 \epsilon^{-3/4} \times 10^-2$ keV, typically 0.26 keV.

	+ This bif corresponds to the transition and can be expressed as a critical power $P_{th} = 1.23 \pi^2 n T_b R a L_t \nu_i$. It can be fitted to:
		$$P_{th} \,\propto\, n B^2 \left(1 + \frac{c}{B n}\right)$$

	+ This ETL model can be used to provide a boundary condition for a core transport code [180]. The bif occurs at a critical value of $T_b$ and on a much faster time scale that $\tau_E$.

---------------------------------------

## 5. Other Theories for the L-H Transition



---------------------------------------

## 6. Discussion



---------------------------------------

## 7. Conclusions



---------------------------------------

## Tables!
<!-- The following is the local L-H transition criteria in terms of dimensionless variables, and also in the form $T_{crit} \,=\, C n^{\alpha_n} B^{\alpha_B}$.

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

-->
