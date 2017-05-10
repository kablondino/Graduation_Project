# Notes on _Bifurcation Theory of the L-H Transition in Fusion Plasmas_
## Chapter 1: Introduction
### 1.2 L-H Transition

+ L-mode: turbulence enhances radial transport

+ H-mode: local reduction in transport (at edge)

	+ The edge therefore has an "edge transport barrier"

+ There is NO consensus on what plasma physics mechanism causes spontaneous transition from L-mode to H-mode

+ 3 different types of dynamics of transition are observed: sharp (sudden), oscillatory (dithering, or I-phase), and smooth

+ Bifurcation theory is the mathematical study of qualitative changes in the solutions of dynamical systems.

### 1.3 Research Questions

+ The MAIN question: How can we employ bifurcation theory to unravel the L-H transition mechanism?

+ Since physical systems can only have a limited number of bifurcations, what bifurcation structure can be recognized in L-H transition dynamics?

+ How can we identify the co-dimension 3 bifurcation in 1-D models?

+ How do these 1-D bifurcating models combine transitions in time and space?

+ What is the bifurcation structure of the model proosed by Zohm?

+ How does the bifurcation structure change when changing the transport reduction mechanism?

+ How does the bifurcation structure change when adding an extra dynamical equation for turbulence? **AND** What is the best dynamical description of the turbulence reduction by sheared $\mathbf{E}\times\mathbf{B}$-flows?

---------------------------------------

## Chapter 2: L-H Transitions in Magnetically-Confined Plasmas
### 2.1 Experimental Observations

+ "The complex, nonlinear behavior of the plasma and very fast timescales of the transition makes it very hard to discriminate the cause and effect relations between the many evolving physical quantities."

	+ Also, usually diagnostics are not quick enough, so they only see the initial L-mode and the final H-mode.

+ Transport in most of the tokamak is the same regardless of mode; only in the edge is turbulence quenched in H-mode.

	+ The $n$, $T$, and $p$ profiles locally steepen inside the transport barrier. This makes H-mode look like L-mode on a pedestal.

+ The most promising observation on why the turbulence is reduced is the acceleration of flows during the transition.

+ The 3 different types of dynamics:

	1. Sharp: most common. When heating is slowly increased above a threshold, it suddenly goes from L- to H-mode; additionally, it will stay in H-mode if the heating is reduced below the initial threshold. This means there is hysteresis characteristics.

	2. Smooth: Low-density discharge in JT-60U and DIII-D with slow heating power increase caused a very smooth transition, with no clear bifurcation.

	3. Oscillatory: self-explanatory, but can take many different forms, with many different names. Each different form may or may not have different mechanisms, or combinations of mechanisms

### 2.2 Related Physical Mechanisms

+ Many proposed mechanisms are one of two things: the trigger of the transition, or the sustaining of the H-mode.

+ One mechanism for sustainment is the reduction of turbulence by sheared flows. The shearing of turbulent eddies until they break into smaller eddies is well known.

+ However, there are different flows in plasmas because of mass difference.

	+ Each species' flows can be decomposed into orthogonal directions: poloidal and toroidal directions or parallel and perpendicular to the magnetic field.

	+ Further, these mass flows are driven by different effects, so flows can be decomposed into driving terms:

		+ Diamagnetic flow, driven by $\nabla p$

		+ $\mathbf{E}\times\mathbf{B}$-flow, driven by the $\mathbf{E}$-field. This is the **SAME** for all particles, independent of mass and charge: $$v_{\mathbf{E}\times\mathbf{B}} \,=\, \frac{\mathbf{E}\times\mathbf{B}}{B^2}$$

		+ According to Burrell's paper of 1994, the $\mathbf{E}\times\mathbf{B}$-flow reduces (ion temperature-driven) turbulence, tearing apart eddies. Also, because the flow works the same on all particles, it can stabilize all possible modes.

		+ This makes it a very universal mechanism for stabilization. There **is** a consensus on this and is accepted. A radial electric field well has be observed near the edge of an H-mode plasma.

	+ The main experimental parameter is heating; heat flux at the outer bits is more relevant. Ions are shown to be the ones to contribute more, since at low $n$, the energy exchange between $e^-$ and ions is limited and lots of ECRH is needed for H-mode.

		+ However, since heating is more financially and practically dealt with, the threshold values need to be found in terms of the ion heat flux.

+ Another decomposition of flows is encountered: the division between zonal and mean flows.

	+ Zonal flows are defined as driven by turbulence itself, and result in small fluctuations of the radial $\mathbf{E}$-field on top of the amount that causes the mean $\mathbf{E}\times\mathbf{B}$-flow.

	+ There is no clear separation in length scales between the two... so it is unclear how they are distinct.

---------------------------------------

## Chapter 3: Bifurcation Theory

+ A bifurcation boundary divides the entire parameter space of a certain model into regions with qualitatively the same type of solutions.

	+ Since L-mode and H-mode are qualitatively distinct, they are most likely separated by a bifurcation.

### 3.1 Introduction

+ A bifurcation is a topological change in the dynamical solution when a small and smooth change of a parameter is made.

+ Dynamical systems are described in terms of DE's as $\dot{x} \,=\, f(x)$, where $x$ is a dynamic variable, and $f(x)$ is an arbitrary function of that variable and probably depends on some parameters.

	+ The classic (and simple) example of a bifurcation is when $f$ is quadratic with $x$, as the following: $\dot{x} \,=\, a + x^2$, where $a$ is the control parameter.

	+ As long as $a$ is negative, there are two solutions: $x_0 \,=\, \pm\sqrt{-a}$ (negative is stable, positive is unstable). However, when $a = 0$, there is only 1 equilibrium point, and is called a saddle-node fixed point. If $a > 0$, there are no equilibrium points. This is known as the _fold bifurcation_.

+ The simplest form of any bifurcating system can be reduced is called the _topological norm form_.

+ The least number of parameters needed to construct this topological norm form is called the _co-dimension_.

+ Bifurcation theory describes BOTH the existence of stationary states but also analyzes the stability of them.

+ Let's introduce a perturbation $x \,\mapsto\, x_0 + x_1 e^{\lambda t}$, and $x_0 \gg x_1$ and $\lambda$ as the eigenvalue, which sign determines if the perturbation grows or not.

	+ For a system of $N$ equations, the eigenvalues become $N$-dimensional complex vectors, with the signs of the real parts determining what is attracting or repelling.

+ The _Hopf bifurcation_ is another with co-dimension 1, with the real part of a pair of complex conjugated eigenvalues vanishing.

	+ This means that the referred steady-state changes from unstable to stable, or vice versa.

	+ _Supercritical Hopf bifucation_ is when there is a stable limit cycle surrounding the steady-state.

	+ _Subcritical_ is when there is an unstable (repelling) limit cycle surrounding the steady-state.

	+ When an unstable limit cycle grows and touches the stable limit cycle, they both vanish, causing a _global fold bifurcation_.

### 3.2 Bifurcations vs L-H Transition Dynamics

+ Different types of observed transition dynamics can be recognized as certain types of bifurcations.

	+ Fold bifurcations are the natural way to describe sharp transitions, such as L-H AND H-L transitions; the natural way to combine two folds is shown in Fig 3.3.

	+ This gives TWO critical values for $a$; this hysteresis behavior is characteristic of two co-dimension 1 fold bifurcations coming from a co-dimension 2 cusp bifurcation.

		+ The norm form of this: $\dot{x} \,=\, a - bx - x^3$; $b$ determines the size of the hysteresis.

	+ The cusp bifurcation separates the parameter space into two regions: one with sharp transitions with hysteresis, and one with smooth transitions without fold bifurcations (with $a$ being what's scanned).

+ Therefore, the cusp bif. organizes two types of L-H transition dynamics: smooth and sharp. Oscillatory behavior can arise due to a Hopf bif.

	+ There is a way that the Hopf bif. can be combined with the cusp bif, which can be view as a _degenerate Bogdanov-Takens bif_. This is **the** co-dimension 3 bif. referenced here.

	+ There are several equivalent unfoldings of this bif., seen as Eq. 3-5

+ The original cusp bif. equation is coupled to a damped variable; as long as the coupled constant $c$ is below a critical value, the parameter space is the same as indicated in Fig. 3.3c.

	+ Any additional coupling will cause the model to be very sensitive to perturbations, making it less robust of a candidate for the L-H transition.

	+ At the critical value of $c$, the co-dimension 3 bif. is encountered at the position of the cusp bif. A regime of limit cycle solutions opens up, covering the original cusp. See Fig. 3.4.

		+ The limit of cycle solutions are produced by Hopf bifs. If both steady-states are unstable, the system will oscillate according to this limit cycle. However, if one steady-state turns stable, the system will transit towards it.

		+ Oscillatory solutions only occur in the region surrounding the cusp bif. point where there are no steady-states.

	+ The specific arrangement of smooth and sharp transitions separated by the oscillatory transitions is characteristic of an underlying co-dimension 3 bif. The quality of an H-mode model can be judged by the existence of the co-dimension 3 bif.

		+ If the model does not contain this bif., it can not describe all different types of L-H transitions.

		+ If a model does describe all of the types of  transitions without having this co-dimension 3 bif., there is a parameter in the model that could pull behavior into inappropriate regions of parameter space, _e.g._ oscillations far away from the L-H transition point. Therefore, co-dimension 3 bif is required to be robust.

---------------------------------------

## Chapter 4: Bifurcation Theory for the L-H Transition in Magnetically-Confined Fusion Plasmas

<!-- ------------ THE ABSTRACT ------------

	The mathematical field of bifurcation theory is extended to be applicable to 1-dimensionally resolved systems of nonlinear partial differential equations, aimed at the determination of a certain specific bifurcation.
	This extension is needed to be able to properly analyze the bifurcations of the radial transport in magnetically confined fusion plasmas.
	This is of special interest when describing the transition from the low-energy-confinement state to the high-energy confinement state of the radial transport in fusion plasmas (i.e. the L-H transition), because the nonlinear dynamical behavior during the transition corresponds to the dynamical behavior of a system containing such a specific bifurcation.
	This bifurcation determines how the three types (sharp, smooth, and oscillating) of observed L-H transitions are organized as function of all the parameters contained in the model.

-->

### 4.1 Introduction

+ Two separate fold bifurcations are necessary to describe the hysteresis, since the heating to trigger the L-H transistion is different than for H-L transition.

+ Two types of parameters affect the existence and magnitude of the hysteresis:

	+ One controls the existence by causing the two fold bifs. to meet in the cusp bif. This causes smooth transitions.

	+ The second causes the hysteresis to be replaced by limit cycle oscillations from a Hopf bif.

	+ These two parameters branch off in parameter space out of the underlying co-dimension 3 bif. The analysis of this bif. it is possible to find how the parameters affect the evolution.

+ The lowest-order system containing this co-dimension 3 bif. is the FitzHugh-Nagumo: $$\dot{x} \,=\, -a - bx - x^3 + cy \\ \dot{y} \,=\, -x - y$$

	+ For $c = 0$, the steady state solutions can have one or multiple possibilities depending on $a$ and $b$.

	+ For $c \neq 0$, the bif. structure stays the same until $c$ is above some critical value, to which the cusp turns into the region of limit cycle solutions (oscillatory solutions).

	+ If a detailed model for the edge transport barrier dynamics contains this co-dimension 3 bif., it is proven that the regions of parameter space with L-mode, H-mode, hysteresis, and dithering are organized in the same way as the FitzHugh-Nagumo model.

### 4.2 Generalized Bifurcation Theory

+ A system of PDEs can be viewed as an infinite system of ODEs, with each ODE describing the evolution of a single point coupled to its neighboring points.

+ Summarizing, to have a general dynamical system with the co-dimension 3 bif., two vectors $\mathbf{v}_1$ and $\mathbf{u}_1$ can be found that satisfy:
	$$M_1 \mathbf{v}_1 \,=\, 0 ~~~~~~~~~~~~~~~~ \mathbf{u}_1^T M_1 \,=\, 0 \\ \mathbf{u}_1^T (M_2 \mathbf{v}_1) \mathbf{v}_1 \,=\, 0 ~~~~~~~~ \mathbf{u}_1^T \cdot \mathbf{v}_1 \,=\, 0$$

+ To obtain the parameter threshold values for the different transitions, the fold bif. is found by analyzing steady-state conditions. The Hopf bif. need further investigation.

	+ The Hopf bif. can be found by unfolding the Bogdanov-Takens bif. in all its parameters and by identifying the specific direction which keeps the eigenvalues purely imaginary, _i.e._ the Hopf condition: $M(\mathbf{v} - \mathbf{v}_{new}) \,=\, i\omega(\mathbf{v} - \mathbf{v}_{new}), ~~~ M = M_1 + \delta M$

	+ After combining, we can come up with a rewritten condition that is invariant under transformation:
	$$\mathbf{u}_1^T \cdot \mathbf{v}_2 (\mathbf{u}_1^T M_3 \mathbf{v}_2 + \mathbf{u}_2^T M_3 \mathbf{v}_1) \,=\, \mathbf{u}_2^T \cdot \mathbf{v}_2 (\mathbf{u}_1^T M_3 \mathbf{v}_1)$$

### 4.3 Finite Dimensional Case

+ Cusp: $a^2 \,=\, -\dfrac{4}{27}(b + c)^3$

+ Hopf: $a^2 \,=\, -\dfrac{4}{27}(b + 1)\left(b + \dfrac{3}{2}c - \dfrac{1}{2}\right)^2$

### 4.4 Transport Model for the L-H Transition

+ Continuity equation for mass (density) and energy, with a single temperature and all the particle and energy deposition into the plasma is somewhere inthe the core outside of the model:
	$$\frac{\partial n}{\partial t} \,=\, \frac{\partial \Gamma}{\partial r} \\ \frac{\partial}{\partial t} \left(\frac{nT}{\gamma - 1}\right) \,=\, -\frac{\partial q}{\partial r}$$

 + Particle and heat flux ($\gamma$ is the adiabatic index):
	$$\Gamma \,=\, -D \frac{\partial n}{\partial r} \\ q \,=\, -\chi n \frac{\partial T}{\partial r} + \frac{\Gamma T}{\gamma - 1}$$

+ Going to higher confinement can be described as a reduction of transport coefficients: particle $D$ and heat $\chi$

+ In the turbulent transport model used here, only a mean flow due to $\mathbf{E}\times\mathbf{B}$-drift is used, such that the transport coefficients become a direct function of the normalized radial electric field:
	$$Z \,=\, \frac{\rho_p e E_r}{T_i}$$

+ To properly describe dynamics, the evolution of the radial electric field must be taken into account via Ampère's law:
	$$\epsilon \frac{\partial Z}{\partial t} \,=\, \mu \frac{\partial^2 Z}{\partial r^2} + c_n \frac{T}{n^2} \frac{\partial n}{\partial r} + \frac{c_T}{n} \frac{\partial T}{\partial r} - G(Z)$$

	+ $\epsilon \,=\, B_p^2 / (B^2 \nu_i)$ is the dielectric constant

	+ The 1st term on the RHS is the radial currents, and $\mu \sim \rho_p^2$, the ratio of viscosity to collision frequency

	+ The 2nd and 3rd terms of the RHS are due to the bipolar part of the anomalous cross field flux, _i.e._ the excess flux of electrons compared to ions.

	+ The last term $G$ is a catch-all for every other effect, and is Taylor expanded because we need an inflection point to obtain the cusp bif.: $G(Z) \,\approx\, a + b(Z - Z_s) + (Z - Z_s)^3$

	+ The model must be the correct size: the outer edge of the plasma (SOL) is fixed at $r = 0$; the inner boundary is at $r = -\infty$; and at $r = \infty$, the $T$ and $n$ are forced to drop towards zero, giving conditions in Eq. 4-38 (or Eq. 5-9).

+ An assumption can be made about the transport coefficients, which allows us to solve the steady-state $n$ and $T$ profiles as a function of particle diffusivity alone:
	$$\chi(Z) \,=\, \frac{D(Z)}{\zeta(\gamma - 1)}$$

	+ The value of the radial electric field is determined by the roots of Eq. 4-43, shown in Eq. 4.44

		+ The LHS is a function purely of the radial electric field, and the RHS is a function of the densition, and is monotonic.

		+ In addition, the edge of the radial electric field is given by Eq. 4-45

		+ An angle $\theta$ can be looked at in the parameter space, shown in Fig 4.5 and Eq. 4-46. This can be considered a single bifurcation parameter.

+ The fold condition is $\mathbf{u}^T_1 M_1 \,=\, 0$, with $M_1$ given in Eq. 4-47. This leads to the following condition:
	$$\frac{\text{d}}{\text{d}Z} \left(\frac{G}{D}\right)\biggr\rvert_e \,=\, 0$$

+ The cusp condition is defined as $\mathbf{u}_1^T (M_2 \mathbf{v}_1) \mathbf{v}_1 \,=\, 0$, which is taken care of by differentiating the 3-tesor $M_2 \,=\, \dfrac{\partial M_1}{\partial \mathbf{v_0}}$. This leads to a similar condition:
	$$\frac{\text{d}^2}{\text{d}Z^2} \left(\frac{G}{D}\right)\biggr\rvert_e \,=\, 0$$

+ The condition for the Bogdanov-Takens bif.:
	$$\int_{-\infty}^0 (u_n v_n + u_T v_T + u_Z v_Z) \, \text{d}x \,=\, 0$$

+ For the Hopf bif., it is quite complicated, with 4 different cases shown in Fig. 4.7; the condition boils down to:
	$$\frac{\text{d}}{\text{d}Z} \left(GD\right)\biggr\rvert_e \,=\, 0$$

+ A complete, but 1-D control parameter space of the bifs. is in Fig. 4.8. The solid curve surrounds the oscillating regime and the large dashed curve corresponds to the sharp hysteresis-like transitions, both part of the Hopf bif. The short-dashed curve corresponds to the two-fold bif. that merges at the cusp bif. point.

### 4.5 Conclusion and Discussion

+ This model, as expected, shows that increasing the heating power increases heat flux from the core, $q(-\infty)$, which decreases $\theta$ towards H-mode.

	+ However, the model also predicts that increasing particle flux has the opposite effect. This is not observed in experiment, because increasing particle flux usually also means increased heating.

	+ Also, this model does not take into account extra momentum and flow due to extra particles, and some SOL physics are excluded.

+ The point of this model is to compare different L-H transition mechanisms, by either doing a different bif. analysis and comparing, or by incorporating new mechanism.

+ The size of the transport barrier and its parametric dependencies can be determined from this model.

+ A list of possible physics topics suggestions to be included in the model is given.

---------------------------------------

## Chapter 5: Bifurcation Theory of a One-Dimensional Transport Model fo the L-H Transition

<!-- ------------ THE ABSTRACT ------------

	Transitions between low and high-confinement (L-H transitions) in magnetically confined plasmas can appear as three qualitatively different types: sharp, smooth, and oscillatory.
	Bifurcation analysis unravels these possible transition types and how they are situated in parameter space.
	In this paper the bifurcation analysis is applied to a 1-dimensional model for the radial transport of energy and density near the edge of magnetically confined plasmas.
	This phenomenological L-H transition model describes the reduction of the turbulent transport by E × B-flow shear self-consistently with the evolution of the radial electric field.
	Therewith, the exact parameter space, including the threshold values of the control parameters, of the possible L-H transitions in the model is determined.
	Furthermore, a generalised equal area rule is derived to describe the evolution of the transport barrier in space and time self-consistently.
	Applying this newly developed rule to the model analysed in this paper reveals a naturally occurring transition to an extra wide transport barrier that may correspond to the improved confinement known as the very- high-confinement mode.

-->

### 5.1 Introduction

+ Some models based on sets of 0-dimensional dynamical equations **can** describe global temporal evolution around L-H transitions, but they lack a description of radial structure of the transport barrier.

	+ There is more need of models that can predict such spatial and temporal observations with their threshold parameters. Bifurcation analysis to the rescue!

+ The considered model assumes L-mode radial transport is dominated by turbulence. It is also assumed that shear in the $\mathbf{E}\times\mathbf{B}$-flow is capable of tearing apart turbulent eddies.

	+ This means it is necessary to include the evolution of the radial electrice field and corresponding flow profile. In this model, there is not small-scale tearing and possible back-reaction of turbulence-generating zonal flows.

+ Generalized Equal-Area Rule, which can be applied to other areas of science. In this context, it applies to the spatial and temporal evolution of the transport barrier.

### 5.2 Transport Model for the L-H Transition

_NOTE:_ This section is very similiar to 4.4, and is only more specific with $\mathbf{E}\times\mathbf{B}$-flow.

+ A well known effect is the reduction of turbulence by the generation of sheared flows, generated by something external or by the turbulence itself.

	+ This kind of self-organizing mechanism could be responsible for the self-sustained transport barrier; the sheared flows are identified as $\mathbf{E}\times\mathbf{B}$-flows.

+ The quenching mechanism is frequently modeled as an effective diffusivity depending on the $\mathbf{E}\times\mathbf{B}$-flow shear:
	$$D \,=\, D_{min} + \frac{D_{max} - D_{min}}{1 + \widetilde{\alpha}\left(V^\prime_{\mathbf{E}\times\mathbf{B}}\right)^2}$$

	+ The prime indicates the radial derivative, and the square of the flow shear means that both signs of flow shear can suppress turbulence.

	+ A similar expression can be used to express the thermal conductivity.

	+ Approximation: $V_{\mathbf{E}\times\mathbf{B}} \approx E_r / B$

+ We cannot expect that the L-H transition will be initiated simply by a difference in the two transport coefficients; therefore, we can make a simplification: $\chi \,=\, D \,/\, \zeta(\gamma - 1)$, with $\zeta$ as a proportionality factor.

	+ This gives the transport equations (5-6a,b) and evolution of the field (5-6c).

### 5.3 Bifurcation Analysis

+ The above is similar to Zohm's model, but differs in the description of effective diffusivity:

	+ Zohm's takes only the value of the radial electric field for diffusivity (shown in Eq. 5-10), with the resulting steady-state density profile as Eq. 5-11.

+ The model for this paper gives a new steady-state in Eq. 5-12

	+ The angle $\theta$ is also found, in Eq. 5-13; Where the two intersect is the state of the system.

+ A description of the transitions is then shown on page 68-70 of the pdf.

+ In Fig. 5.3, the parameter space of both models are plotted with the same values (_a_ is shear, _b_ is Zohm):

	+ In Zohm's model, the oscillations during an oscillatory L-H transition last a lot longer.

	+ The onset of the oscillatory behavior for Zohm's is at values with a wider range, where the flow-shear model would have already gone into H-mode.

	+ The values for $\theta(L-H)$ are generally higher (_i.e._ lower heating threshold) for the flow-shear model, supporting that sheared flow is more efficient in reducing turbulence.

### 5.4 The Transport Barrier: Space and Time Consistently

+ Temporal transitions are a first-order derivative, and spatial transitions are a second-order derivative (Eq. 5-14):
	$$-\epsilon\frac{\partial X}{\partial t} \,+\, \mu \frac{\partial^2 X}{\partial r^2} \,=\, F(X) - c(r,t)$$

+ Temporal transitions between roots correspond to sudden jumps from L- to H-mode and back; spatial transition is when the core of the plasma exhibits L-mode-like transport and the edge exhibits H-mode-like transport.

	+ The limit of purely temporal or purely spatial transitions are well known.

	+ The limit of $\mu \rightarrow 0$ gives maximum hysteresis for the temporal transition.

	+ The limit of $\epsilon \rightarrow 0$ describes how a high-transport core can be connected to a low-transport edge.

		+ In the time-independent case, the equation can be integrated over space, which must vanish at $X_+$, which leads to Maxwell's equal area rule (Eqs. 5-15 and 5-16).

+ For the whole system, assume the jumps in time ans pace are rapid ($\epsilon, \mu \ll 1$), so that transitions happen in an almost 1-D zone in $(r,t)$-space.

	+ This lets us use $\dfrac{\text{d}}{\text{d}t} \,\rightarrow\, v \dfrac{\text{d}}{\text{d}r}$

	+ This allows a new function $K(X)$, which gives the generalized equal area rule (see Eqs. 5-17 through 5-19)
	+ The GEA rule determines the position in space and time of the transition between L- and H-mode transport corresponding to the temporal growth of the barrier region.

_NOTE:_ This section was definitely confusing.

### 5.5 Two Different Regimes in Transport Barrier Sizes

+ This section applies the GEA rule to Zohm's model; Eq. 5-11 is rewritten in the form of Eq. 5-14, shown in 5-20.

+ Figure 5.6 shows the typical evolution of the solution, as $\theta_1 \rightarrow \theta_2 \rightarrow \theta_3 \rightarrow \theta_4$, with the subfigures showing the state. Read the caption of the figure for more info.

	+ The transition moves into the plasma to build a _thick-barrier H-mode_ only once $\theta$ crosses the GEA condition. Else, the H-mode is only in the _thin-barrier_ regime.

		+ The thick barrier has width increasing with input power, while the thin barrier has a constant width as a function of input power.

	+ The thin-barrier width is of the order of the viscocity, _e.g._ several gyro-radii.

### 5.6 Conclusion and Discussion

+ The thick-barrier is hypothesized to be the mechanism responsible for VH-mode; however, tokamaks might already reach other limits before the required heating power is reached, so this is not a claim.

---------------------------------------

## Chapter 6: Comparison of Bifurcation Dynamics of Turbulent Transport Models for the L-H Transition

<!-- ------------ THE ABSTRACT ------------

	In more than three decades a large amount of models and mechanisms have been proposed to describe a very beneficial feature of magnetically confined fusion plasmas: the L-H transition.
	Bifurcation theory can be used to compare these different models based on their dynamical transition structure.
	In this paper we employ bifurcation theory to distinguish two fundamentally different descriptions of the interaction between turbulence levels and sheared flows.
	The analytic bifurcation analysis characterises the parameter space tructure of the transition dynamics.
	Herewith, in these models three dynamically different types of transitions are characterised, sharp transitions, oscillatory transitions and smooth transitions.
	One of the two models has a very robust transition structure and is therefore more likely to be more accurate for such a robust phenomenon as the L-H transition.
	The other model needs more fine-tuning to get non-oscillatory transitions.
	These conclusions from the analytic bifurcation analysis are confirmed by dedicated numerical simulations, with the newly developed code Bifurcator.

-->

### 6.1 Introduction



### 6.2 Turbulent Transport Models for the L-H Transition

+ The minimum amount of transport is determined by neoclassical effects; anomalous transport depends on the turbulence level on top of that.

	+ The anomalous transport increases linearly with the turbulence level $\mathcal{E}$:
		$$D \,=\, D_{min} + \left(D_{max} - D_{min}\right)\frac{\mathcal{E}}{\mathcal{E}_{max}}, ~~~~~~~~ \mathcal{E}_{max} \,=\, \frac{\gamma_L}{\alpha_{sat}}$$

		+ $\mathcal{E}_{max}$ is the steady-state turbulence level without any flow shear suppresion; $\gamma_L$ is the linear growth rate of the turbulence; $\alpha_{sat}$ depends on the saturation mechanism corresponding to the turbulence.

+ Eqs. 6-4 through 6-6 show the evolution of the turbulence, either linearly (Eq. 6-5), or nonlinearly (Eq. 6-6); both are cited.

+ The evolution of the radial electric field is determined by the sum of all possible radial currents at the edge: $\epsilon_0 \dfrac{\partial E_r}{\partial t} \,=\, -\sum J_r$. All possible mechanisms for generating $J_r$ is summed in Appendix 6.A.

_NOTE:_ The rest of the section has been covered. It talks about the evolution of the field, along with $G(Z)$ and the conditions.

### 6.3 Bifurcation Analysis

+ The difference in the two models (linear vs nonlinear) arises due to the turbulence level ODE. The steady-state differences are given. See Eqs. 6-15 and 6-16. The results for $D$ of which are given in 6-17 and 6-18.

	+ Linear: whether there is turbulence only affects the range of $Z$ in which it is stable.

	+ Nonlinear: no turbulence results in an always unstable system, and nonzero turbulence is stable.

	+ The bifurcation analysis has been done on these, and are expected to be equivalent. However, with turbulence, the bifurcation structure does change(?).

+ In the nonlinear model, $b$ (from the Taylor expansion of $G$) affects what type of transition will occur. In the linear model, $b$ has no effect, and the transition always has an oscillatory phase. See Figs. 6.3 and 6.4.

+ The nonlinear turbulence model is in steady-state exactly the same as the flow-shear model.

+ Fig. 6.5 summarizes much about the possible transition dynamics. The following figures show smaller values of $\alpha$ for the linear model.

+ Summary: similar transition dynamics can be found in both models. However, the linear model is very sensitive to $\alpha$; in contrast, the nonlinear model is very robust.

### 6.4 Bifurcator

+ Bifurcator is a numerical solver for nonlinear ODEs, optimized for bifurcating systems. PDEs need to be discretized first to be used.

+ It uses various implicit Runge-Kutta methods; implicit meaning the time integration involves solving nonlinear systems.

	+ Using Newton iteration and requires the user to define the Jacobian matrix of the discretized system.

+ A bif. detection scheme is implemented: it obtains the steady-state solution for a given set of parameter and varies one of them over a interval with user-defined increment.

### 6.5 Numerical Bifurcation Analysis

+ For a parameter scan, it is important to start every simulation from the steady-state profile of the previous step.

+ Simulations of the nonlinear model shows L-modes where the radial electric field profile is close to zero everywhere and the turbulence $\mathcal{E}$ is close to $\mathcal{E}(max)$, and H-modes where the electric field well is formed near the edge and locally the turbulence is reduced. (AS EXPECTED)

+ Read this section again.

### 6.6 Conclusion and Discussion



### 6.A Appendix: Radial Currents



---------------------------------------

## Chapter 7: Evaluation and Future Prospects

### 7.1 Conclusions and Discussion



### 7.2 Outlook



---------------------------------------

## Summary



## My Own Questions

+ Could the oscillatory transition explain ELMs?
