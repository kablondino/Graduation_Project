# Notes on _Backstepping control of H-L confinement transition in fusion plasma model_
## Chapter 2: Physics of the confinement transition and H-mode
### 2.1 Order of scales

+ Transport of magnetic flux is the slowest process within a tokamak plasma; therefore, the magnetic field can be treated as a constant background field during transition.

+ Magnetic reconnection (_e.g._ sawtooth instability) occurs on a fast time scale, but is not considered for this thesis.

+ Pressure gradients drive the turbulent modes, giving anomalous transport of particles, energy, and momentum.

	+ Changes in turbulence occur in the smallest of time scales; they are often inferred from fluctuations in density and electric field.

### 2.2 Turbulence
#### 2.2.1 Dynamics

+ The dominant modes are the ion temperature gradient mode and trapped electron mode in a collisionless plasma.

+ Micro-instabilities give drift waves because electrons neutralize fluctuations in finite time.

	+ The coupling of drift waves to other plasma modes yields different types of drift-wave turbulences such as the drift resistive ballooning mode and drift Alfvén turbulence.

+ Drift waves are the only type of turbulence able to drive a radial flux and therefore, drift wave turbulence determines the level of edge turbulence and associated transport. Basic model for this:
	$$\frac{\text{d}\mathcal{E}}{\text{d}t} \,=\, \gamma_1 \mathcal{E} - \alpha_{sat}\mathcal{E}^2$$

	+ $\mathcal{E}$ is the turbulence level, $\gamma_1$ is the growth rate, and $\alpha_{sat}$ is the saturation rate._

#### 2.2.2 Flow shear suppression of turbulence

+ A radial shear in plasma flow decorrelates turbulent eddies, which are small compared to the flow length scale. As eddies are extended by the sheared plasma flow, their wave energy is transferred to larger flow structures ad the expense of small turbulent structures:
	$$\frac{\partial V_{\mathbf{E}\times\mathbf{B}}}{\partial x} \,=\, \frac{\partial}{\partial x} \left(\frac{E_r}{B}\right) \,\approx\, \frac{1}{B_\phi} \frac{\partial E_r}{\partial x}$$

+ Both zonal and mean $\mathbf{E}\times\mathbf{B}$ flow are involved in the transition.

	+ Zonal flows are produced by the shear of Reynolds stress arising from correlated radial and poloidal fluctuations.

		+ Its shear suppresses the drift-wave turbulence, but the mechanism is damped by ion-ion collisions or flow instabilities. Therefore, zonal flows provide a temporary sink of energy, which show their signature before the transition.

	+ The mean flow remains in H-mode, while the zonal flow vanishes. Shear flow suppression 
