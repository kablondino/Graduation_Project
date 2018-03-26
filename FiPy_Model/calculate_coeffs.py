"""
	This file contains the $g$ coefficients and plasma
	parameters for use in the alternate $Z$ equation
	developed by Staps.
	It has been written as a function.
"""

from boundary_init_cond import *

import scipy.special		# For the Faddeeva (plasma dispersion) function
import numpy

## Plasma Parameters
#e_p = 1.0 + (m_i*density_adjusted + m_e*density_adjusted) / (epsilon_0 * B**2)

def update_g_coeffs():
	# Neutrals density in use for CX friction
	n_0.setValue(4.0e17 * a_in0 * (temperature / 100.0)\
			**(3.0/4.0))										# [m^-3]

	# Thermal velocities (most probable)
	v_Ti.setValue((2.0 * temperature / m_i)**(1.0/2.0))			# [m/s]
	v_Te.setValue((2.0 * temperature / m_e)**(1.0/2.0))			# [m/s]

	# Poloidal gyro-(Larmor) radii
	rho_pi.setValue(m_i * v_Ti / (charge * B_theta))		# [m]
	rho_pe.setValue(m_e * v_Te / (charge * B_theta))		# [m]

	# Banana orbit bounce frequencies
	omega_bi.setValue(aspect**(3.0/2.0) * v_Ti / (q * R))		# [s^-1]
	omega_be.setValue(aspect**(3.0/2.0) * v_Te / (q * R))		# [s^-1]

	# Collision frequencies within electrons and ions
	nu_ei.setValue(1.33e5*(density*1.0e-20)\
			/ (temperature)**(3.0/2.0))			# [s^-1]
	nu_ii.setValue(1.2 * (m_e / m_i)**(1.0/2.0) * nu_ei)		# [s^-1]

	# Collision frequency of trapped ions and neutrals
	nu_in0.setValue(a_in0 * omega_bi)							# [s^-1]

	# Effective detrapping frequency
	nu_eff.setValue(nu_ii + nu_in0)								# [s^-1]

	# Effective collision frequencies
	nu_ai.setValue(nu_ii / omega_bi)	# nu_*i
	nu_ae.setValue(nu_ei / omega_be)	# nu_*e (not used as of now)


	## Electron Anomalous Diffusion
	D_an.setValue((aspect)**2*(pi)**(1.0/2.0) / (2*a_m) *\
			(rho_pe * temperature)/B)
	g_n_an.setValue(density*D_an)
	g_T_an.setValue(g_n_an * alpha_an)
	g_Z_an.setValue(g_n_an / rho_pi)

	Gamma_an.setValue((g_n_an*density.grad[0]\
			/density + g_T_an*temperature.grad[0]/temperature\
			+ g_Z_an*Z))										# [m^-2 s^-1]


	## Charge Exchange Friction
	g_n_cx.setValue((-(m_i*n_0*neu_react_rate * density\
			* temperature) / (charge*B_theta**2))\
			* ((B_theta**2 / (aspect*B_phi)**2) + 2.0))
	g_T_cx.setValue(alpha_cx * g_n_cx)
	g_Z_cx.setValue(-g_n_cx / rho_pi)

	Gamma_cx.setValue((g_n_cx * density.grad[0]/density\
			+ g_T_cx * temperature.grad[0]/temperature\
			+ g_Z_cx*Z))										# [m^-2 s^-1]


	## Ion Bulk (Parallel) Viscosity
	bulk_complex_term = 1j*(pi)**(1.0/2.0) * scipy.special.wofz(\
			Z + 1j*nu_ii*aspect*B / (v_Ti*B_theta))
#	print numpy.imag(bulk_complex_term)

	Gamma_bulk.setValue( aspect**2*density*temperature\
			/ (pi**(1.0/2.0)*B*x*charge) * (rho_pi / 0.5 + Z)\
			* numpy.imag(bulk_complex_term) )					# [m^-2 s^-1]
#	print Gamma_bulk


	## Ion Orbit Loss
	g_OL.setValue((charge * density * nu_eff\
			* (aspect)**(1.0/2.0) * rho_pi))

	Gamma_OL.setValue(numerix.exp(-(nu_ai + Z**4)**(1.0/2.0))\
			/ (nu_ai + Z**4)**(1.0/2.0))						# [m^-2 s^-1]

