# Notes on _Edge radial electric field structure and its connections to H-mode confinement in Alcator C-Mod plasmas_
## I. Introduction

+ H-mode plasmas have a negative radial electric field well forms just inside the last-closed flux surface (LCFS)

+ On DIII-D, the evolution of the field has been observed prior to the transition, by measuring an increase in the edge impurity poloidal velocity in the electron diamagnetic direction.

+ The details of $E_r$ well structure vary from machine to machine, although many report similar observations.

+ _EDA_ stands for Enhanced D-Alpha H-mode, and are quiescent, without ELMs but have increased particle transport.

## II. The Edge CXRS Diagnositic on Alcator C-Mod

$$E_r \,=\, \frac{1}{n_i Z_i e} \frac{\text{d}p_i}{\text{d}r} - V_{\theta i} B_\phi + V_{\phi i} B_\theta$$

+ The above ignores the contributions from the Reynolds stress and off-diagonal terms in the pressure tensor.

	+ Measurements of the impurity ion dist. function indicate that it is Maxwellian, justifying the simplification of the stress tensor to the gradient of the pressure.

	+ The three terms are:

		+ Diamagnetic $\dfrac{1}{n_i Z_i e} \dfrac{\text{d}p_i}{\text{d}r}$

		+ Poloidal velocity $V_{\theta i} B_\phi$

		+ Toroidal velocity $V_{\phi i} B_\theta$

	+ To calculate $E_r$ accurately, one needs the toroidal and poloidal magnetic fields and velocities, as well as the temperature and density for a particular species.

## III. $E_r$ Behaviour in H-mode Plasmas
### A. EDA H-Modes

+ The EDA H-modes have the radial field forming within 10 mm inside the LCFS, with the full-width half-max (FWHM) between 4 and 6 mm.

+ At smaller minor radii, the radial field is always positive and exhibits a smooth, slowly-varying profile.

+ Shearing rate:
	$$\omega_{\mathbf{E}\times\mathbf{B}} \,=\, \frac{r}{q(r)} \frac{\text{d}}{\text{d}r} \left[\frac{q(r)}{r} \frac{E}{B}\right]$$

+ Generally, the contributions to $E_r$ in EDA H-modes are similar to each other.

+ The poloidal velocity well depth is usually greater and narrower than the diamagnetic well depth. This means the poloidal velocity contribution tends to set the $E_r$ well width and the $\mathbf{E}\times\mathbf{B}$ shearing rate.

+ Inside the pedestal, the positive radial field is contributed to by a combo of cocurrent toroidal rotation and slightly ion diamagnetic poloidal flow. These far outweigh the slightly negative contribution from the diamagnetic component.

	+ In the core, the field remains positive as a result, which is definitely dominant.

### B. ELM-free H-modes

+ These are inherently transient because impurities accumulate, which causes radiative collapse.

+ The wells can be significantly deeper than the EDA H-mode wells. This also means that the ELM-free modes can have energy confinement at 30% higher and substantially higher particle confinement than EDA modes.

	+ The width of the well, however, is the same, causing a much larger shearing rate, up to 5 times as high.

	+ The depth is not constant, however. It degrades as expected.

+ Due to particle confinement that is too high, the temperature goes down. This causes the shearing rate to decline, which indicates that the $\mathbf{E}\times\mathbf{B}$ shear is linked more strongly to the edge energy barrier than to the particle barrier.

	+ Also, an additional mechanism is needed to explain this; $\mathbf{E}\times\mathbf{B}$ shear is not enough.

## IV. Comparison to Main Ion $E_r$ Estimates

+ Since the radial electric field is the same for all species, info about the main ion population can be inferred by comparing the total radial field to an estimate of the main ion diamagnetic contribution.

	+ The difference between $E_r$ and the main ion diamagnetic contribution provides a measure of the main ion flow perpendicular to the magnetic field.

	+ Also this comparison can give info about the source or formation physics of the field, suggesting guidelines for creating self-consistent field profiles for modeling codes.
