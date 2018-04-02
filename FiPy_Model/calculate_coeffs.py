"""
	This file contains the $g$ coefficients and plasma
	parameters for use in the alternate $Z$ equation
	developed by Staps.
	It has been written as a function as to be called within
	the solving loop.
"""

from variable_decl import *

import scipy.special		# For the Faddeeva (plasma dispersion) function
import numpy


# ASSUMES density is in m^-3 and temperature is in eV
def calculate_coeffs():
	# Neutrals density in use for CX friction, NEEDS CHANGE!
	n_0.setValue(4.0e17 * a_in0 * (charge*temperature / 100.0)\
			**(3.0/4.0))										# [m^-3]

	# Thermal velocities (most probable)
	v_Ti.setValue(numerix.sqrt(2.0 * charge * temperature / m_i))	# [m/s]
	v_Te.setValue(numerix.sqrt(2.0 * charge * temperature / m_e))	# [m/s]

	# Poloidal gyro-(Larmor) radii
	rho_pi.setValue(m_i * v_Ti / (charge * B_theta))			# [m]
	rho_pe.setValue(m_e * v_Te / (charge * B_theta))			# [m]

	# Transition frequency
	omega_t.setValue(v_Ti / (q*R))

	# Banana orbit bounce frequencies
	omega_bi.setValue(aspect**(3.0/2.0) * omega_t)				# [s^-1]
	omega_be.setValue(aspect**(3.0/2.0) * v_Te / (q * R))		# [s^-1]

	# Banana width
	w_bi.setValue(numerix.sqrt(aspect) * rho_pi)				# [m]

	# Collision frequencies within electrons and ions
	nu_ei.setValue(4.2058e-11*(density)\
			/ (temperature)**(3.0/2.0))							# [s^-1]
	nu_ii.setValue(1.2 * numerix.sqrt(m_e / m_i) * nu_ei)		# [s^-1]

	# Collision frequency of trapped ions and neutrals
	nu_in0.setValue(a_in0 * omega_bi)							# [s^-1]

	# Effective detrapping frequency
	nu_eff.setValue(nu_ii + nu_in0)								# [s^-1]

	# Effective collision frequencies
	nu_ai.setValue(nu_ii / omega_bi)	# nu_*i
	nu_ae.setValue(nu_ei / omega_be)	# nu_*e (not used as of now)


	## Electron Anomalous Diffusion
	D_an.setValue(aspect**2 * numerix.sqrt(pi) * rho_pe * temperature\
			/ (2*a_m * B))
	g_n_an.setValue(charge*density*D_an)						# [A m^-2]
	g_T_an.setValue(g_n_an * alpha_an)							# [A m^-2]
	g_Z_an.setValue(g_n_an / rho_pi)							# [A m^-1]

	Gamma_an.setValue((g_n_an*density.grad[0]\
			/ density + g_T_an*temperature.grad[0]/temperature\
			+ g_Z_an*Z) / charge)								# [m^-2 s^-1]


	## Charge Exchange Friction
	g_n_cx.setValue((-(m_i*n_0*neu_react_rate * density\
			* temperature) / (B_theta**2))\
			* ((B_theta**2 / (aspect*B_phi)**2) + 2.0))			# [A m^-2]
	g_T_cx.setValue(alpha_cx * g_n_cx)							# [A m^-2]
	g_Z_cx.setValue(-g_n_cx / rho_pi)							# [A m^-1]

	Gamma_cx.setValue((g_n_cx * density.grad[0]/density\
			+ g_T_cx * temperature.grad[0]/temperature\
			+ g_Z_cx*Z) / charge)								# [m^-2 s^-1]


	## Ion Bulk (Parallel) Viscosity
	bulk_complex_term = 1j * numerix.sqrt(pi) * scipy.special.wofz(\
			Z + 1j*nu_ii / omega_t)
	D_bulk.setValue(aspect**2 * rho_pi * temperature\
			/ (x * B * numerix.sqrt(pi)))						# [m^2 s^-1

	Gamma_bulk.setValue(density * D_bulk * (L + Z/rho_pi)\
			* numpy.imag(bulk_complex_term) )					# [m^-2 s^-1]


	## Ion Orbit Loss
	g_OL.setValue(charge * density * nu_ii * nu_ai * rho_pi)	# [A m^-2]

	radical_OL = numerix.sqrt(nu_ai + Z**4 + ((x)/w_bi)**4)

	Gamma_OL.setValue(g_OL * numerix.exp(-radical_OL)\
			/ (charge*radical_OL))								# [m^-2 s^-1]

	Z_transient_coeff.setValue(m_i * density * temperature\
			/ (charge* rho_pi * B**2))
	Z_diffusion_coeff.setValue(m_i * mu * density * temperature\
			/ (charge* rho_pi * B_theta**2))

