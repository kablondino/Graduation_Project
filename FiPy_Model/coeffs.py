"""
	This file contains the $g$ coefficients and plasma
	parameters for use in the alternate $Z$ equation
	developed by Staps.
	It has been written as a function.
"""

from variable_decl import *

## Plasma Parameters
e_p = 1.0 + (m_i*density + m_e*density) / (epsilon_0 * B**2)

# Neutrals density in use for CX friction
n_0 = CellVariable(name=r"$n_0$", mesh=mesh)

# Thermal velocities (most probable)
v_Ti = CellVariable(name=r"$v_{th,i}$", mesh=mesh)
v_Te = CellVariable(name=r"$v_{th,e}$", mesh=mesh)

# Poloidal gyro-(Larmor) radii
rho_pi = CellVariable(name=r"$\rho_{\theta i}$", mesh=mesh)
rho_pe = CellVariable(name=r"$\rho_{\theta e}$", mesh=mesh)

# Banana orbit bounce frequencies
omega_bi = CellVariable(name=r"$\omega_{bi}$", mesh=mesh)
omega_be = CellVariable(name=r"$\omega_{be}$", mesh=mesh)

# Collision Frequencies
nu_ei = CellVariable(name=r"$\nu_{ei}$", mesh=mesh)
nu_ii = CellVariable(name=r"$\nu_{ii}$", mesh=mesh)
nu_in0 = CellVariable(name=r"$\nu_{i0}$", mesh=mesh)
nu_eff = CellVariable(name=r"$\nu_{eff}$", mesh=mesh)
nu_ai = CellVariable(name=r"$\nu_{*i}$", mesh=mesh)
nu_ae = CellVariable(name=r"$\nu_{*e}$", mesh=mesh)

## Electron Anomalous Diffusion
D_an = CellVariable(name=r"$D_{an}$", mesh=mesh)
g_n_an = CellVariable(name=r"$g_n^{an}$", mesh=mesh)
g_T_an = CellVariable(name=r"$g_T^{an}$", mesh=mesh)
g_Z_an = CellVariable(name=r"$g_Z^{an}$", mesh=mesh)
Gamma_an = CellVariable(name=r"$e \Gamma^{an}$", mesh=mesh)

## Charge Exchange Friction
g_n_cx = CellVariable(name=r"$g_n^{cx}$", mesh=mesh)
g_T_cx = CellVariable(name=r"$g_T^{cx}$", mesh=mesh)
g_Z_cx = CellVariable(name=r"$g_Z^{cx}$", mesh=mesh)
Gamma_cx = CellVariable(name=r"$e \Gamma^{cx}$", mesh=mesh)

## Ion Bulk (Parallel) Viscosity
# Consolidated constant to reduce clutter, listed as N on reference
N = CellVariable(name="Consolidated Bulk Viscosity value", mesh=mesh)
# xi_p integral
xi_p = CellVariable(name=r"$\xi_\theta$", mesh=mesh)
# xi_t integral
xi_t = CellVariable(name=r"$\xi_\phi$", mesh=mesh)

g_n_bulk = CellVariable(name=r"$g_n^{\pi\parallel}$", mesh=mesh)
g_T_bulk = CellVariable(name=r"$g_T^{\pi\parallel}$", mesh=mesh)
g_Z_bulk = CellVariable(name=r"$g_Z^{\pi\parallel}$", mesh=mesh)
Gamma_bulk = CellVariable(name=r"$e \Gamma^{\pi\parallel}$", mesh=mesh)

## Ion Orbit Loss
g_OL = CellVariable(name=r"$g^{OL}$", mesh=mesh)
f_OL = CellVariable(name=r"$f^{OL}$", mesh=mesh)

def update_g_coeffs():
	e_p = 1.0 + (m_i*density + m_e*density) / (epsilon_0 * B**2)

	# Neutrals density in use for CX friction
	n_0.setValue(4.0e17 * a_in0 * (temperature / 100.0)**(3.0/4.0))

	# Thermal velocities (most probable)
	v_Ti.setValue((2.0 * charge * temperature / m_i)**(1.0/2.0))
	v_Te.setValue((2.0 * charge * temperature / m_e)**(1.0/2.0))

	# Poloidal gyro-(Larmor) radii
	rho_pi.setValue(m_i * v_Ti / (charge * B_theta))
	rho_pe.setValue(m_e * v_Te / (charge * B_theta))

	# Banana orbit bounce frequencies
	omega_bi.setValue(aspect**(3.0/2.0) * v_Ti / (q * R))
	omega_be.setValue(aspect**(3.0/2.0) * v_Te / (q * R))

	# Collision frequencies within electrons and ions
	nu_ei.setValue((density / temperature)**(3.0/2.0))
	nu_ii.setValue(1.2 * (m_e / m_i)**(1.0/2.0) * nu_ei)

	# Collision frequency of trapped ions and neutrals
	nu_in0.setValue(a_in0 * omega_bi)

	# Effective detrapping frequency
	nu_eff.setValue(nu_ii + nu_in0)

	# Effective collision frequencies
	nu_ai.setValue(nu_ii / omega_bi)	# nu_*i
	nu_ae.setValue(nu_ei / omega_be)	# nu_*e (not used as of now)

	## Consolidated constant to reduce clutter, listed as N on reference
	N.setValue(value = (nu_ai * aspect**(3.0/2.0) * nu_ei / nu_ii))

	## Electron Anomalous Diffusion
	D_an.setValue((aspect)**2*(pi)**(1.0/2.0) / (2*a_m) *\
			(rho_pe * temperature)/B)
	g_n_an.setValue(-charge_true*density*D_an)
	g_T_an.setValue(g_n_an * alpha_an)
	g_Z_an.setValue(g_n_an / rho_pi)
	Gamma_an.setValue((g_n_an*density.grad[0]/density\
			+ g_T_an*temperature.grad[0]/temperature + g_Z_an*Z))

	## Charge Exchange Friction
	g_n_cx.setValue((-(m_i*n_0*neu_react_rate * density * temperature)\
			/ (charge_true*B_theta**2))\
			* ((B_theta**2 / (aspect*B_phi)**2) + 2.0))
	g_T_cx.setValue(alpha_cx * g_n_cx)
	g_Z_cx.setValue(-g_n_cx / rho_pi)
	Gamma_cx.setValue((g_n_cx*density.grad[0]/density\
			+ g_T_cx*temperature.grad[0]/temperature + g_Z_cx*Z))

	## Ion Bulk (Parallel) Viscosity
	# xi_p integral
	xi_p.setValue((4.0*N*(27.0*((N)**2 + Z**2)**2 - 7.0*((N)**2 - 3.0*Z**2)*(nu_ai)**(1.0/2.0))*(nu_ai)**(7.0/4)) / (189.0*pi*((N)**2 + Z**2)**3))
	# xi_t integral
	xi_t.setValue((2.0*N*(135.0*((N)**2 + Z**2)**2 - 7.0*(21.0*(N)**4 + 3.0*Z**2*(-5.0 + 7.0*Z**2) + (N)**2*(5.0 + 42.0*Z**2))*(nu_ai)**(1.0/2))*(nu_ai)**(7.0/4)) / (189.0*pi*((N)**2 + Z**2)**3))

	g_n_bulk.setValue(aspect**2*(pi)**(1.0/2) / (8*a_m) * density * m_i * rho_pi * (v_Ti)**2 * B_theta * xi_p)
	g_T_bulk.setValue(aspect**2*(pi)**(1.0/2) / (8*a_m) * density * m_i * rho_pi * (v_Ti)**2 * (B_theta*xi_p - B*xi_t))
	g_Z_bulk.setValue(aspect**2*(pi)**(1.0/2) / (4*a_m) * density * m_i * (v_Ti)**2 * B_theta*xi_p)
	Gamma_bulk.setValue((1.0/charge_true) * (g_n_bulk*density.grad[0]/density\
			+ g_T_bulk*temperature.grad[0]/temperature + g_Z_bulk*Z))

	## Ion Orbit Loss
	g_OL.setValue((charge * density * nu_eff * (aspect)**(1.0/2.0) * rho_pi))
	f_OL.setValue(numerix.exp(-(nu_ai + Z**4)**(1.0/2.0))/ (nu_ai + Z**4)**(1.0/2.0))

