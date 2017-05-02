# Notes on _Bifurcation Theory of the L-H Transition in Fusion Plasmas_
## Chapter 1: Introduction
### 1.2 L-H Transition

+ L-mode: turbulence enhances radial transport

+ H-mode: local reduction in transport (at edge)

	+ The edge therefore has an "edge transport barrier"

+ There is NO consensus on what plasma physics mechanism causes spontaneous transition from L-mode to H-mode

+ 3 different types of dynamics of transition are observed: sharp (sudden), oscillatory, and smooth

+ Bifurcation theory is the mathematical study of qualitative changes in the solutions of dynamical systems.

### 1.3 Research Questions

+ The MAIN question: How can we employ bifurcation theory to unravel the L-H transition mechanism?

+ Since physical systems can only have a limited number of bifurcations, what bifurcation structure can be recognized in L-H transition dynamics?

+ How can we identify the co-dimension 3 bifurcation in 1-D models?

+ How do these 1-D bifurcating models combine transitions in time and space?

+ What is the bifurcation structure of the model proosed by Zohm?

+ How does the bifurcation structure change when changing the transport reduction mechanism?

+ How does the bifurcation structure change when adding an extra dynamical equation for turbulence? **AND** What is the best dynamical description of the turbulence reduction by sheared $\mathbf{E}\times\mathbf{B}$-flows?

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

## Chapter 4: Bifurcation Theory for the L-H Transition in Magnetically-Confined Fusion Plasmas
### 4.1 Introduction

+ Two separate fold bifurcations are necessary to describe the hysteresis, since the heating to trigger the L-H transistion is different than for H-L transition.

+ Two types of parameters affect the existence and magnitude of the hysteresis:

	+ One controls the existence by causing the two fold bifs. to meet in the cusp bif. This causes smooth transitions.

	+ The second causes the hysteresis to be replaced by limit cycle oscillations from a Hopf bif.

	+ These two parameters branch off in parameter space out of the underlying co-dimension 3 bif. The analysis of this bif. it is possible to find how the parameters affect the evolution.

+ The lowest-order system containing this co-dimension 3 bif. is the Fitz-Hugh-Nagumo: $$\dot{x} \,=\, -a - bx - x^3 + cy \\ \dot{y} \,=\, -x - y$$

	+ For $c = 0$, the steady state solutions can have one or multiple possibilities depending on $a$ and $b$.

	+ For $c \neq 0$, the bif. structure stays the same until $c$ is above some critical value, to which the cusp turns into the region of limit cycle solutions (oscillatory solutions).

	+ If a detailed model for the edge transport barrier dynamics contains this co-dimension 3 bif., it is proven that the regions of parameter space with L-mode, H-mode, hysteresis, and dithering are organized in the same way as the Fitz-Hugh-Nagumo model.

### 4.2 Generalized Bifurcation Theory

+ A system of PDEs can be viewed as an infinite system of ODEs, with each ODE describing the evolution of a single point coupled to its neighboring points.
