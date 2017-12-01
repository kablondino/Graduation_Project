"""
	This file contains the $g$ coefficients for use in the 
	alternate $Z$ equation developed by Staps.
	It NEEDS to be inserted after the declaration of the
	variables.
"""

## Global parameters and physical constants
pi = 3.141592653589793
charge = 1.602e-19	# Elementary charge
k_B = 8.617e-5		# Boltzmann in eV
m_e = 9.109e-31		# Electron mass
m_i = 1.673e-27		# Ion (H) mass
epsilon_0 = 8.854187817e-12	# Permittivity of free space
mu_0 = 4*pi*1.0e-7			# Permeability of free space

## ASDEX-U specifications
a_v = 0.8	# Vertical minor radius
a_h = 0.5	# Horizontal minor radius
a_m = ( (a_v**2 + a_h**2) / 2.0 )**(1.0/2.0)
R = 1.6		# Major radius
aspect = a_m / R
I_p = 2.0e6
B_t = 3.9
B_p = mu_0 * I_p / ( 2*pi*a_m )
q = aspect * B_t/B_p
B = ( B_t**2 + B_p**2 )**(1.0/2.0)

## Plasma Parameters
def e_p(X): return 1.0 + (m_i*n(X) + m_e*n(X)) / (epsilon_0 * B**2)

def v_Ti(X): return (2.0 * charge * Temp(X) / m_i)**(1.0/2.0)

def v_Te(X): return (2.0 * charge * Temp(X) / m_e)**(1.0/2.0)

def rho_pi(X): return m_i * v_Ti(X) / (charge * B_p)

def rho_pe(X): return m_e * v_Te(X) / (charge * B_p)

def omega_bi(X): return aspect**(3.0/2.0) * v_Ti(X) / (q * R)

def omega_be(X): return aspect**(3.0/2.0) * v_Te(X) / (q * R)

## Collision Frequencies
def nu_ei(X): return 1.33e5 * (n(X) * 1.0e-20) / (Temp(X) * 1.0e-3)**(3.0/2.0)

def nu_ii(X): return 1.2 * (m_e / m_i)**(1.0/2.0) * nu_ei(X)

def nu_in0(X): return a_in0 * omega_bi(X)

def nu_eff(X): return nu_ii(X) + nu_in0(X)

def nu_ai(X): return nu_ii(X) / omega_bi(X)	# nu_*i

def nu_ae(X): return nu_ei(X) / omega_be(X)	# nu_*e (not used as of now)

## Consolidated constant to reduce clutter, listed as N on reference
def N(X): return nu_ai(X) * aspect**(3.0/2.0) * nu_ei(X) / nu_ii(X)

## Electron Anomalous Diffusion
def D_an(X): return (aspect)**2*(pi)**(1.0/2.0) / (2*a_m) * (rho_pe * temperature)/B
def g_n_an(X): return -charge*density*D_an
def g_T_an(X): return g_n_an * alpha_an
def g_Z_an(X): return g_n_an / rho_pi

## Charge Exchange Friction
def g_n_cx(X): return -(m_i*n_0*neu_react_rate) / (B_theta**2)
def g_T_cx(X): return alpha_cx * g_n_cx(X)
def g_Z_cx(X): return -g_n_cx(X) / rho_pi(X)

## Ion Bulk (Parallel) Viscosity
# xi_p integral
def xi_p(X, aZ): return (4*C(X)*(27*((C(X))**2 + aZ**2)**2 - 7*((C(X))**2 - 3*aZ**2)*(nu_ai(X))**(1.0/2))*(nu_ai(X))**(7.0/4)) / (189*pi*((C(X))**2 + aZ**2)**3)
# xi_t integral
def xi_t(X, aZ): return (2*C(X)*(135*((C(X))**2 + aZ**2)**2 - 7*(21*(C(X))**4 + 3*aZ**2*(-5 + 7*aZ**2) + (C(X))**2*(5 + 42*aZ**2))*(nu_ai(X))**(1.0/2))*(nu_ai(X))**(7.0/4)) / (189*pi*((C(X))**2 + aZ**2)**3)

def g_n_bulk(X, aZ): return aspect**2*(pi)**(1.0/2) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * B_p * xi_p(X, aZ)
def g_T_bulk(X, aZ): return aspect**2*(pi)**(1.0/2) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * (B_p*xi_p(X, aZ) - B*xi_t(X, aZ))
def g_Z_bulk(X, aZ): return aspect**2*(pi)**(1.0/2) / (4*a_m) * n(X) * m_i * (v_Ti(X))**2 * B_p*xi_p(X, aZ)

## Ion Orbit Loss
def g_OL(X, aZ) = charge * density * nu_eff * (aspect)**(1.0/2.0) * rho_pi
def f_OL(X, aZ) = g_OL * exp(-(nu_ai + Z**4)**(1.0/2.0)) / (nu_ai + Z**4)**(1.0/2.0)

