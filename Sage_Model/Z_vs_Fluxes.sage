reset()

load("parameters.sage")

var('y,Z')

assume(Z, 'real')


# Physical Constants
charge = scipy.constants.e
m_i = scipy.constants.m_p
m_e = scipy.constants.m_e

# ----------------- Input of State ------------------------
## t = 200
# At x ~ 0.002
#x	$n$	$T$	$Z$	$D$	$D_{an}$	$\Gamma_e^{an}$	$n_0$	$\langle\sigma_{cx} v\rangle$	$\Gamma_i^{cx}$	$D_{\pi\parallel}$	$\Gamma_i^{\pi\parallel}$	$\Gamma_i^{ol}$
#0.001875	5.25011171139421e+17	163.482319322465	1.45659916620318	1.95092542714814	0.000688141630100473	2.27690887456374e+17	56501480943441	2.53812246363164e-13	-4.81176944459356e+16	-0.0188249902028863	-1.19168187618349e+18	-557173039481322
#0.002125	5.35434214493223e+17	165.427742871548	1.47166341450206	2.13128525822014	0.000700479137552669	2.30504108694339e+17	56167794581565.9	2.5481649909684e-13	-4.83290085496657e+16	-0.0191697028228131	-1.17180463596314e+18	-466734643459656
#0.002375	5.45017799197483e+17	167.165974949308	1.48552149414179	2.31149230840068	0.000711564705839756	2.33494397297697e+17	55874584307141.2	2.55707178968072e-13	-4.86155187455718e+16	-0.0194804005875539	-1.15395821152485e+18	-377374818968410

# At x ~ 0.03
#0.029875	1.08143075661057e+18	221.574027780129	1.77718267908245	4.94713184277656	0.00108613644485018	6.09230942553967e+17	2496379657.78461	2.80904355987589e-13	-5407156833206.3	-0.0310182685391704	-1.37348458709749e+18	-1.55235354380374e-58
#0.030125	1.08615114035334e+18	221.843028685801	1.77821178953312	4.94743894562273	0.00108811581270749	6.12848344919829e+17	1943024953.06425	2.81018035060206e-13	-4232723648741.13	-0.0310869925862421	-1.37764340545142e+18	-1.15809327369854e-59
#0.030375	1.09087529545641e+18	222.110163544139	1.77923651112516	4.94775127669252	0.00109008262082869	6.16471134861609e+17	1512332115.80433	2.81130833603446e-13	-3313304089634.74	-0.0311554116345023	-1.38179703421831e+18	-8.48298023345065e-61


x_before = 0.001875
x = 0.002125
x_after = 0.002375
density_before = 5.25011171139421e+17
density = 5.35434214493223e+17
density_after = 5.45017799197483e+17
density_grad = (density_after - density_before) / (x_after - x_before)
temperature_before = 163.482319322465
temperature = 165.427742871548
temperature_after = 167.165974949308
temperature_grad = (temperature_after - temperature_before) / (x_after - x_before)

Gamma_c = -1.0e20



# ----------------- Values --------------------------------
# Thermal Velocity
v_Ti = sqrt(2.0 * charge * temperature / m_i)					# [m/s]
v_Te = sqrt(2.0 * charge * temperature / m_e)					# [m/s]

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


# Effective collision frequencies (Collisionalities)
nu_ai = nu_ii / omega_bi	# nu_*i
nu_ae = nu_ei / omega_be	# nu_*e (not used as of now)


## Electron Anomalous Diffusion
D_an = aspect**2 * sqrt(pi) * rho_pe * temperature\
		/ (2*a_m * B)
g_n_an = charge * density * D_an								# [A m^-2]
g_T_an = g_n_an * alpha_an										# [A m^-2]
g_Z_an = g_n_an / rho_pi										# [A m^-1]

Gamma_an(Z) = g_n_an * density_grad / density  +  g_T_an * temperature_grad / temperature  +  g_Z_an * Z


## Charge Exchange Friction
n_0 = (-0.1*Gamma_c / v_Ti) / (1.0 + exp(1.0e3*(x-0.02)))		# [m^-3]
cx_rate = 1.0e-14 * (100*temperature)^(1.0 / 3.0)				# [m^3 s^-1]

g_n_cx = -(m_i * n_0 * cx_rate * density * temperature) / (B_theta^2) * ((B_theta^2 / (aspect * B_phi)^2) + 2.0)			# [A m^-2]
g_T_cx = alpha_cx * g_n_cx										# [A m^-2]
g_Z_cx = g_n_cx / rho_pi										# [A m^-1]

Gamma_cx(Z) = g_n_cx * density_grad / density  +  g_T_cx * temperature_grad / temperature  +  g_Z_an * Z


## Ion Bulk (Parallel) Viscosity
plasma_disp(Z) = imag(I*sqrt(pi) * e^(-(Z + I*nu_ii/omega_t)^2) * erfc(-I*(Z + I*nu_ii/omega_t)))
#plasma_disp(Z) = exp((nu_ii/omega_t)^2 - Z^2) * (sin(2*nu_ii*Z/omega_t) * imag(erfc(i*Z - nu_ii/omega_t)) + cos(2*nu_ii*Z/omega_t) * real(erfc(i*Z - nu_ii/omega_t)))
# Naive plasma dispersion:
#plasma_disp = sqrt(pi) * exp(-Z^2)
D_bulk = aspect^2 * rho_pi * temperature / ((x - a_m) * B * sqrt(pi))	# [m^2 s^-1]

Gamma_bulk(Z) = density * D_bulk * (density_grad / density  +  Z / rho_pi) * plasma_disp	# [m^-2 s^-1]


## Ion Orbit Loss
g_ol = -charge * density * nu_ii * nu_ai * rho_pi				# [A m^-2]
radical_ol = sqrt(nu_ai + (Z)^4 + ((x) / w_bi)^4)

Gamma_ol(Z) = g_ol * exp(-radical_ol) / (charge * radical_ol)	# [m^-2 s^-1]


## Sum of all the fluxes
#Gamma_sum = sage.plot.misc.setup_for_eval_on_grid(Gamma_an + Gamma_cx + Gamma_bulk + Gamma_ol, [(-7,7)], plot_points=10000)

# ----------------- Plotting ------------------------------
Zmin, Zmax = -6, 6
plot_title = "$x = " + str(x)+ "$"
Gamma_an_plot = plot(Gamma_an, (Z, Zmin, Zmax), legend_label=r"$\Gamma_e^{an}$", color='red', gridlines=True)
Gamma_cx_plot = plot(Gamma_cx, (Z, Zmin, Zmax), legend_label=r"$\Gamma_i^{cx}$", color='green')
#Gamma_bulk_plot = plot(Gamma_bulk, (Z, Zmin, Zmax), legend_label=r"$\Gamma_i^{\pi\parallel}$", color='purple')
Gamma_ol_plot = plot(Gamma_ol, (Z, Zmin, Zmax), legend_label=r"$\Gamma_i^{ol}$", color='blue')

Gamma_sum_plot = plot(Gamma_an(Z) + Gamma_cx(Z) + Gamma_ol(Z) + Gamma_bulk(Z), (Z, Zmin, Zmax), legend_label=r"$\sum_k \Gamma^k$", color='magenta', gridlines=True, ymin=-4.0e15, ymax=1.0e14)
Gamma_sum_plot.show()

combined_plot = Gamma_an_plot + Gamma_cx_plot + Gamma_ol_plot + Gamma_bulk_plot
combined_plot.show()

