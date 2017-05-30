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

#### (ii) Spatial and temporal characteristics

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

+ The results of taking each assumption is given by Eqs. 2.24 through 2.27.

### 2.3 _The shear flow paradigm_

+ $V_E^\prime$ denotes the radial derivative of the equilibrium $\mathbf{E}\times \mathbf{B}$ velocity, $V_E$.

+ When turbulence is isotropic ($\Delta r \sim r\Delta \theta$), it can be written that:
	$$\omega_s \,=\, \frac{RB_\theta}{B} \frac{\partial}{\partial r}\left(\frac{E_r}{RB_\theta}\right)$$

### 2.4 _Radial electric fields_

+ $V_E^\prime$ is determined by the radial electric field $E_r$, which is then deduced from the radial force balance for any plasma species $j$:
	$$E_r \,=\, -\frac{1}{n_j e_j} \frac{\text{d} p_j}{\text{d} r} + V_{\theta j} B_\phi - V_{\phi j} B_\theta$$

	+ $e_j$ is charge, $p_j$ is pressure, $n_j$ is density, and $V_{\theta j}$ and $V_{\phi j}$ are the poloidal and toroidal velocities.

	+ Therefore, changes in $E_r$ can be associated with changes in $V_{\theta j}$, $V_{\phi j}$, or $p_j^\prime$ (the prime is the radial derivative).

	+ The ones that vary in the poloidal velocity correspond to bifurcations in the solutions of the poloidal momentum balance equation, requiring a momentum source/torque and sink.

		+ Sources: ion-orbit loss, non-ambipolar electron loss, Stringer spin-up due to poloidal asymmetries in turbulent transport, Reynolds stress from turbulence, etc.

		+ Sinks: ion parallel neoclassical viscosity, charge exchange on neutrals, etc.

	+ Models involving $V_{\phi j}$ or $p_j^\prime$ depend on sources of particles (neutrals surrounding the plasma, neutral beams), energy (heating), or toroidal momentum (neutral beams) driving $E_r^\prime$ until a transport bifurcation occurs.

---------------------------------------

## 3. Instabilities and Turbulence at the Plasma Edge



---------------------------------------

## 4. Models Involving Sheared Radial Electric Field and Flows



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


