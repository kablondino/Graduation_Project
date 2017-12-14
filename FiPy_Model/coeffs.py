"""
	This file contains the $g$ coefficients and plasma
	parameters for use in the alternate $Z$ equation
	developed by Staps.
	It NEEDS to be inserted after the declaration of the
	variables and x.
"""

from parameters import *
from variable_decl import *
from fipy.tools import numerix

## Plasma Parameters
# Neutrals density in use for CX friction
n_0 = 4.0e17 * a_in0 * (temperature / 100.0)**(3.0/4.0)

e_p = 1.0 + (m_i*density + m_e*density) / (epsilon_0 * B**2)

v_Ti = (2.0 * charge * temperature / m_i)**(1.0/2.0)

v_Te = (2.0 * charge * temperature / m_e)**(1.0/2.0)

rho_pi = m_i * v_Ti / (charge * B_theta)

rho_pe = m_e * v_Te / (charge * B_theta)

omega_bi = aspect**(3.0/2.0) * v_Ti / (q * R)

omega_be = aspect**(3.0/2.0) * v_Te / (q * R)

## Collision Frequencies
nu_ei = 1.33e5 * (density * 1.0e-20) / (temperature * 1.0e-3)**(3.0/2.0)

nu_ii = 1.2 * (m_e / m_i)**(1.0/2.0) * nu_ei

nu_in0 = a_in0 * omega_bi

nu_eff = nu_ii + nu_in0

nu_ai = nu_ii / omega_bi	# nu_*i

nu_ae = nu_ei / omega_be	# nu_*e (not used as of now)

## Consolidated constant to reduce clutter, listed as N on reference
N = nu_ai * aspect**(3.0/2.0) * nu_ei / nu_ii

## Electron Anomalous Diffusion
D_an = (aspect)**2*(pi)**(1.0/2.0) / (2*a_m) * (rho_pe * temperature)/B
g_n_an = -charge*density*D_an
g_T_an = g_n_an * alpha_an
g_Z_an = g_n_an / rho_pi

## Charge Exchange Friction
g_n_cx = -(m_i*n_0*neu_react_rate) / (B_theta**2)
g_T_cx = alpha_cx * g_n_cx
g_Z_cx = -g_n_cx / rho_pi

## Ion Bulk (Parallel) Viscosity
# xi_p integral
xi_p = (4*N*(27*((N)**2 + Z**2)**2 - 7*((N)**2 - 3*Z**2)*(nu_ai)**(1.0/2))*(nu_ai)**(7.0/4)) / (189*pi*((N)**2 + Z**2)**3)
# xi_t integral
xi_t = (2*N*(135*((N)**2 + Z**2)**2 - 7*(21*(N)**4 + 3*Z**2*(-5 + 7*Z**2) + (N)**2*(5 + 42*Z**2))*(nu_ai)**(1.0/2))*(nu_ai)**(7.0/4)) / (189*pi*((N)**2 + Z**2)**3)

g_n_bulk = aspect**2*(pi)**(1.0/2) / (8*a_m) * density * m_i * rho_pi * (v_Ti)**2 * B_theta * xi_p
g_T_bulk = aspect**2*(pi)**(1.0/2) / (8*a_m) * density * m_i * rho_pi * (v_Ti)**2 * (B_theta*xi_p - B*xi_t)
g_Z_bulk = aspect**2*(pi)**(1.0/2) / (4*a_m) * density * m_i * (v_Ti)**2 * B_theta*xi_p

## Ion Orbit Loss
g_OL = charge * density * nu_eff * (aspect)**(1.0/2.0) * rho_pi
f_OL = g_OL * numerix.exp(-(nu_ai + Z**4)**(1.0/2.0)) / (nu_ai + Z**4)**(1.0/2.0)

