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

---------------------------------------
