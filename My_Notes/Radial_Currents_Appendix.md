# Notes on Chapter 6 Appendix of Weymiens's Thesis

$$\epsilon_0 \frac{\partial E_r}{\partial t} \,=\, -\sum J_r \,=\, e\sum (\Gamma_e - \Gamma_i)$$

## 1. Neoclassical Polarization Current
### Notes on _Neoclassical Dielectric Property of a Tokamak Plasma_ by F.L. Hinton and J.A. Robertson

<!-- ------------ THE ABSTRACT ------------

	The response of a tokamak plasma to a time-dependent axisymmetric radial electric field is considered.
	The bounce-average motion of magnetically trapped ions is shown to consist of a radial drift, in addition to the well-known toroidal precession.
	This is a neoclassical polarization drift, which is larger than the standard one by a factor of B^2 / B_\theta^2, the square of the ratio of total to poloidal magnetic fields.
	The resulting low-frequency dielectric constant is larger than the standard one by approximately the same factor, when the ions are in the banana regime of neoclassical theory.

-->

+ I. Introduction

There are several mechanisms which would set up [a quasi-electrostatic field], such as the loss of energetic ions or runaway electrons.

Time-dependent radial electric field in the context of momentum input from NBI. The plasma toroidal velocity is related to the radial electric field by:
	$$u_{\phi} \,=\, \frac{c E_r}{B_\theta}$$

Angular momentum input in the same direction as the plasma current (co-injection) results in an increasing radial electric field, while counter-injection results in a decreasing one.

Macroscopically, a radial current must flow in the plasma to provide the torque to change the plasma angular momentum, with the required value of the current density as ($\rho$ is mass density):
	$$\frac{j_r B_\theta}{c} \,=\, \rho \frac{\partial u_\phi}{\partial t}$$

This give a relationship between $j_r$ and $\dot{E_r}$, which **is** the neoclassical polarization current:
	$$j_r \,=\, \frac{\rho c^2}{B_\theta^2} \frac{\partial E_r}{\partial t}$$

It is argued that both trapped and free (untrapped) particles participate with the same radial drift velocity.

The expression for current can be combined with Maxwell's equation, with the radial component of $\nabla\times\mathbf{B}$ averaged to zero, and $j_r^\text{ext}$ as any external current (NBI):
	$$\frac{\partial E_r}{\partial t} \,=\, -4\pi(j_r - j_r^\text{ext})$$

---------------------------------------

## 2. Shear Viscosity
### Notes on _Effect of Electrode Biasing on the Radial Electric Field Structure Bifurcation in Tokamak Plasmas_ by N. Kasuya, K. Itoh, and Y. Takase

<!-- ------------ THE ABSTRACT ------------

	The mechanism for formation of a steep structure in the radial electric field is a key issue in plasma confinement.
	Properties of the radial electric field bifurcation are studied taking into account the effect of electrode biasing.
	The radial electric field structure is determined by the charge conservation equation.
	From the nonlinear mechanism associated with local current due to ion bulk viscosity, a transition can take place.
	Various types of radial electric field structures with multiple peaks are allowed for the same boundary condition.
	The ion orbit loss term breaks the symmetry of the radial current similarly to the ambipolar radial electric field.
	A radial current driven by the electrode plays the role of a control parameter in a transition similarly to the pressure gradient.
	A phase diagram is given in the spontaneous drive vs external drive space.
	Differences in the radial shape of solitary electric field structures are demonstrated in the presence of spatial varying components.
	This study clarifies the mechanisms of nonlinear structure formation in transport barriers.

-->

+ "In this paper, properties of the radial electric field transition between L-mode and H-mode are studied, taking into accout the effect of electrode biasing."

Electrode biasing is a method in controlling the radial electric field externally.

This model also includes neoclassical bulk viscosity and orbit losses as contributors. These two components have different dependencies on the radial electric field, and their competition gives structural bifurcation of the field.

+ The radial electric field structure is determined by the charge conservation equation:
	$$\frac{\partial E_r}{\partial t} \,=\, -\frac{1}{\epsilon_0 \epsilon_\perp} \left(J_\text{visc} + J_r + J_\text{other} - J_\text{ext}\right)$$ with $\epsilon_\perp$ as the dielectric constant of magnetized plasma, $J_\text{visc}$ is current driven by shear viscocity, $J_r$ is local current due to ion bulk viscosity, and $J_\text{other}$ encompassing other, such as orbit losses and charge exchange losses.

The shear viscocity current can be written (with $\mu_i$ as the shear viscocity of ions) as:
	$$J_\text{visc} \,=\, -\epsilon_0\epsilon_\perp \nabla\cdot \mu_i \nabla E_r$$

+ The charge conservation equation above can be normalized as the following:
	$$\frac{\partial^2 X}{\partial x^2} - (X - X_a) f(X,y) - g(X,y) + I = 0 \\
	X \,=\, \frac{E_r}{v_{i,therm} B_\theta}, ~~~~~ y \,=\, \frac{r \nu_{ii} B}{v_{i,therm} B_\theta}, ~~~~~ I \,=\, \frac{J_\text{ext}}{v_{i,therm} B_\theta \sigma(0)}, \\
	x \,=\, \frac{r - r_0}{l}, ~~~\text{and}~~~ l \,=\, \sqrt{\frac{\mu_i \epsilon_0 \epsilon_\perp}{\sigma(0)}}$$
$X$ is (another) normalized radial electric field, $X_a$ is the ambipolar radial electric field, $\sigma(0)$ is the conductivity for zero radial electric field, $y$ is the normalized collision frequency, $I$ is the normalized external current, and $x$ is the minor radius normalized by $l$. $r_0$ is chosen to be the mid-point between the electrode and limiter.

	+ The first term in the normalized DE corresponds to the shear viscocity and acts as a diffusion term.

	+ The second term arises from $J_r$ and is given by the neoclassical ion particle flux.

	+ The function $f(X,y)$ is the electrical conductivity, which is nonlinearly-dependent on the radial electric field.

	+ The function $g(X,y)$ comes from $J(\text{other})$, specifically ion orbit losses.

	+ $X_a$ depends on the gradient of each plasma parameter. See Eq. (6). This term is used as the control parameter that describes the plasma.

+ The number of particles in the loss region gives the magnitude of the orbit loss term:
	$$g(X,y) \,=\, g_0(X,y) \,=\, A \frac{y}{\Pi} \, \exp\left(-\Pi\right) ~~~~~\text{with}~~~~\Pi \,\equiv\, \sqrt{\epsilon^{-3/2} y + (\alpha X)^4}$$

	+ $\alpha = 0.5$ is the numerical constant that accounts for orbit shape, $\epsilon = 0.33$ is the inverse aspect ratio, and $A = 3$ is a constant that represents the ratio of the magnitude of the orbit loss term to the neoclassical term.

Two nonlinear terms affect the radial electric field in this case. The ion orbit loss term breaks the symmetry of the radial current similarly to the ambipolar radial electric field that depends on $\nabla p$.

As $\nabla p$ increases, a spontaneous sharp transition takes place when the collision frequency is small and squeezing of the banana orbit is strong.

The radial electric field is negative when $I = 0$. When there is current in the electrode, it acts like a $\nabla p$. The threshold for transition by electrode biasing changes in accordance with $\nabla p$.

+ Previously, we assumed the ambipolar radial electric field term to be constant in space, so the DE had translational invariance in space. In reality, there is spatial variation in the plasma parameter that affect the structure of the radial electric field **through** the ambipolar radial electric field.

+ Also, the magnitude of the orbit loss terms depends on the spatial position (assume a simple form of exponential decay):
	$$g(X,y,x) \,=\, g_0(X,y) \, \exp\left[\left(-\frac{x - d}{\rho_p}\right)^2\right]$$

---------------------------------------
