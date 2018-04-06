load("parameters.sage")

assume(x,'real')
assume(x >= 0)

var('y,Z')

assume(Z, 'real')

Gamma_c = -0.8e20

# STATE!
density(x) = (0.5e19 / L)*x + 0.5e19							# [m^-3]
temperature(x) = 3.0e3*x + 100.0								# [eV]


v_Ti = sqrt(2.0 * charge * temperature / m_i)					# [m/s]
v_Te = sqrt(2.0 * charge * temperature / m_e)					# [m/s]

# Neutrals density in use for CX friction, NEEDS CHANGE!
# The commented version cannot be properly computed
#n_0 = (-0.1*Gamma_c / v_Ti) * numerix.exp(\
#		-sqrt(sigma_ion*cx_rate) * density[0]\
#		* numerix.exp(x * density.grad[0]/density)))			# [m^-3]
# NEED dynamic definition!
n_0 = (-0.1*Gamma_c / v_Ti) / (1.0 + exp(1.0e3*(x-0.01)))		# [m^-3]

# Poloidal gyro-(Larmor) radii
rho_pi = m_i * v_Ti / (charge * B_theta)						# [m]
rho_pe = m_e * v_Te / (charge * B_theta)						# [m]

# Transition frequency
omega_t = v_Ti / (q*R)											# [s^-1]

# Banana orbit bounce frequencies
omega_bi = aspect**(3.0/2.0) * omega_t							# [s^-1]
omega_be = aspect**(3.0/2.0) * v_Te / (q * R)					# [s^-1]

# Banana width
w_bi = sqrt(aspect) * rho_pi									# [m]

# Collision frequencies within electrons and ions
nu_ei = 4.2058e-11*(density)\
		/ (temperature)**(3.0/2.0)								# [s^-1]
nu_ii = 1.2 * sqrt(m_e / m_i) * nu_ei							# [s^-1]

# Collision frequency of trapped ions and neutrals
nu_in0 = a_in0 * omega_bi										# [s^-1]

# Effective detrapping frequency
nu_eff = nu_ii + nu_in0											# [s^-1]

# Effective collision frequencies
nu_ai = nu_ii / omega_bi	# nu_*i
nu_ae = nu_ei / omega_be	# nu_*e (not used as of now)


## Electron Anomalous Diffusion
D_an = aspect**2 * sqrt(pi) * rho_pe * temperature\
		/ (2*a_m * B)
g_n_an = charge*density*D_an									# [A m^-2]
g_T_an = g_n_an * alpha_an										# [A m^-2]
g_Z_an = g_n_an / rho_pi										# [A m^-1]

Gamma_an(Z) = (g_n_an*diff(density, x) / density + g_T_an* diff(temperature, x)/temperature + g_Z_an*Z) / charge									# [m^-2 s^-1]


## Charge Exchange Friction
g_n_cx = (-(m_i*n_0*cx_rate * density\
		* temperature) / (B_theta**2))\
		* ((B_theta**2 / (aspect*B_phi)**2) + 2.0)				# [A m^-2]
g_T_cx = alpha_cx * g_n_cx										# [A m^-2]
g_Z_cx = -g_n_cx / rho_pi										# [A m^-1]

Gamma_cx(Z) = (g_n_cx * diff(density, x)/density + g_T_cx * diff(temperature, x)/temperature + g_Z_cx*Z) / charge									# [m^-2 s^-1]


## Ion Bulk (Parallel) Viscosity
#plasma_disp = (I * sqrt(pi)\
#		* scipy.special.wofz(Z + 1j*nu_ii / omega_t)).imag()
plasma_disp(y) = I*sqrt(pi) * e^(-y^2) * erfc(-I*y)
D_bulk = aspect**2 * rho_pi * temperature\
		/ ((a_m - x) * B * sqrt(pi))							# [m^2 s^-1]

Gamma_bulk(Z) = density * D_bulk * (1.0/L + Z/rho_pi) * (plasma_disp.subs(y=Z + I*nu_ii / omega_t)).imag()	# [m^-2 s^-1]


## Ion Orbit Loss
g_OL = -charge * density * nu_ii * nu_ai * rho_pi				# [A m^-2]

radical_OL(Z) = sqrt(nu_ai + Z**4 + ((x+0.01)/w_bi)**4)

Gamma_OL(Z) = g_OL * exp(-radical_OL) / (charge*radical_OL)		# [m^-2 s^-1]

Z_transient_coeff = m_i * density * temperature\
		/ (charge* rho_pi * B**2)
Z_diffusion_coeff = m_i * mu * density * temperature\
		/ (charge* rho_pi * B_theta**2)

