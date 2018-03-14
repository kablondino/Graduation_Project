"""
	This file contains the $g$ coefficients and plasma
	parameters for use in the alternate $Z$ equation
	developed by Staps.
	It has been written as a function.
"""

import numpy
load("parameters.sage")

assume(x,'real')
assume(x >= 0)

Z(x) = Z_S*(1.0 - tanh((L*x - L) / 2.0))

# ----------------- Diffusivities -------------------------
# Itohs'/Zohm's model
D_Zohm = (D_max + D_min) / 2.0 + ((D_max - D_min)*tanh(Z)) / 2.0
# Stap's Model
alpha_sup = 0.5
D_Staps = D_min + (D_max - D_min) / (1.0 + alpha_sup*(diff(Z,x))^2)
# Flow-Shear Model
a1, a3 = 1.0, 0.5	# ASSUMES a2 = 0
D_Shear = D_min + (D_max - D_min) / (1.0 + a1*(Z)^2 + a3*(diff(Z,x)^2))

density_si_coeff = 1.0e19		# Adjusts to m^-3
temp_si_coeff = 250.0			# Adjusts to eV

## INITAL density and temperature
density(x) = -(density_si_coeff*Gamma_c*lambda_n / D_Staps) * (1.0 + x/lambda_n)

temperature(x) = temp_si_coeff*q_c*((gamma - 1.0) / Gamma_c) * (1.0 - lambda_n / (zeta*lambda_T + lambda_n) * (1.0 + x/lambda_n)^(-zeta))

#state_plots = plot(density, (x,0,L), color='blue', legend_label=r"$n$")\
#		+ plot(temperature, (x,0,L), color='orange', legend_label=r"$T$")\
#		+ plot(Z, (x,0,L), color='green', legend_label=r"$Z$")
#
#show(state_plots + plot(D_Staps, (x,0,L), color='red', legend_label=r"$D$"))

# ----------------- Plasma Parameters ---------------------
e_p = 1.0 + (m_i*density + m_e*density) / (epsilon_0 * B^2)

# Neutrals density in use for CX friction
n_0 = 4.0e17 * a_in0 * (temperature / 100.0)^(3.0/4.0)

# Thermal velocities (most probable)
v_Ti = sqrt(2.0 * charge * temperature / m_i)
v_Te = sqrt(2.0 * charge * temperature / m_e)

# Poloidal gyro-(Larmor) radii
rho_pi = m_i * v_Ti / (charge * B_theta)
rho_pe = m_e * v_Te / (charge * B_theta)

# Banana orbit bounce frequencies
omega_bi = aspect^(3.0/2.0) * v_Ti / (q * R)
omega_be = aspect^(3.0/2.0) * v_Te / (q * R)

# Collision frequencies within electrons and ions
nu_ei = (density / temperature)^(3.0/2.0)
nu_ii = 1.2 * sqrt(m_e / m_i) * nu_ei

# Collision frequency of trapped ions and neutrals
nu_in0 = a_in0 * omega_bi

# Effective detrapping frequency
nu_eff = nu_ii + nu_in0

# Effective collision frequencies
nu_ai = nu_ii / omega_bi	# nu_*i
nu_ae = nu_ei / omega_be	# nu_*e (not used as of now)

## Consolidated constant to reduce clutter, listed as N on reference
#N = value = (nu_ai * aspect^(3.0/2.0) * nu_ei / nu_ii)

# ----------------- Fluxes --------------------------------
## Electron Anomalous Diffusion
D_an = (aspect)^2*sqrt(pi) / (2*a_m) *\
		(rho_pe * temperature) / B
g_n_an = -charge_true*density*D_an
g_T_an = g_n_an * alpha_an
g_Z_an = g_n_an / rho_pi
Gamma_an(x) = g_n_an*diff(density,x)/density + g_T_an*diff(temperature,x)/temperature + g_Z_an*Z

## Charge Exchange Friction
g_n_cx = (-(m_i*n_0*neu_react_rate * density * temperature)\
		/ (charge_true*B_theta^2))\
		* ((B_theta^2 / (aspect*B_phi)^2) + 2.0)
g_T_cx = alpha_cx * g_n_cx
g_Z_cx = -g_n_cx / rho_pi
Gamma_cx(x) = g_n_cx*diff(density,x)/density + g_T_cx*diff(temperature,x)/temperature + g_Z_cx*Z

## Ion Bulk (Parallel) Viscosity ----> OLD!
# xi_p integral
#xi_p = (4.0*N*(27.0*((N)^2 + Z^2)^2 - 7.0*((N)^2 - 3.0*Z^2)*(nu_ai)^(1.0/2.0))*(nu_ai)^(7.0/4)) / (189.0*pi*((N)^2 + Z^2)^3)
## xi_t integral
#xi_t = (2.0*N*(135.0*((N)^2 + Z^2)^2 - 7.0*(21.0*(N)^4 + 3.0*Z^2*(-5.0 + 7.0*Z^2) + (N)^2*(5.0 + 42.0*Z^2))*(nu_ai)^(1.0/2))*(nu_ai)^(7.0/4)) / (189.0*pi*((N)^2 + Z^2)^3)

#g_n_bulk = aspect^2*sqrt(pi) / (8*a_m) * density * m_i * rho_pi * (v_Ti)^2 * B_theta * xi_p
#g_T_bulk = aspect^2*sqrt(pi) / (8*a_m) * density * m_i * rho_pi * (v_Ti)^2 * (B_theta*xi_p - B*xi_t)
#g_Z_bulk = aspect^2*sqrt(pi) / (4*a_m) * density * m_i * (v_Ti)^2 * B_theta*xi_p
#Gamma_bulk(x) = (1.0/charge_true) * (g_n_bulk*diff(density,x)/density + g_T_bulk*diff(temperature,x)/temperature + g_Z_bulk*Z)

var('z')
plasma_disp(z) = I*sqrt(pi)*exp(-z^2) * erfc(-I*z)
bulk_imag_term(x) = plasma_disp(Z + I*nu_ii*aspect*B / (v_Ti * B_theta))

bulk_imag_term_array(x) = [ plasma_disp(Z(k) + I*nu_ii(k)*aspect*B / (v_Ti(k) * B_theta)) for k in numpy.arange(0.0, L, 0.1)]

# SageMath function definition takes WAY too long
#Kobayashi_Gamma_bulk(x) = aspect^2*density*temperature / (sqrt(pi)*B*x) * (rho_pi / 0.5 + Z) * bulk_imag_term.imag()

## Ion Orbit Loss
g_OL = (charge * density * nu_eff * sqrt(aspect) * rho_pi)
Gamma_OL(x) = g_OL * exp(-sqrt(nu_ai + Z^4))/ sqrt(nu_ai + Z^4)

# ----------------- PLOTS ---------------------------------
#plot(density, (x,0,L), color='blue', axes_labels=[r"$x$", r"$10^{-3}$ m"],
#		title="Density").show()
#plot(temperature, (x,0,L), color='orange', axes_labels=[r"$x$", r"$e$V"],
#		title="Temperature").show()
#plot(5.0e-30*Z*temperature / (charge_true*rho_pi), (x,0,L), color='green',
#		axes_labels=[r"$x$", "kV/m"], title="Electric Field").show()

