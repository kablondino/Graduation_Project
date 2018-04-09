"""
	This file contains the simple constant delcarations
	for use in the PDE system. This includes machine
	parameters, length scales, global physical/math
	constants, and non-g model numbers.

"""

from input_handling import *


# ----------------- Constant Parameters -------------------
## Global parameters and physical constants
pi = 3.141592653589793
charge_dummy = 1.0							# 'Natural' units for charge
charge = 1.6021766208e-19					# Elementary charge
k_B_J = 1.38064852e-23						# Boltzmann in J/K
k_B_eV = k_B_J / charge						# Boltzmann in eV/K
m_e = 9.10938356e-31						# Electron mass
m_i = 1.673e-27								# Ion (H) mass
epsilon_0 = 8.854187817e-12					# Permittivity of free space
mu_0 = 4*pi*1.0e-7							# Permeability of free space
c = 1 / (epsilon_0*mu_0)**(1.0/2.0)			# Speed of light
gamma = 5.0/3.0								# Adiabatic index for monoatomic


## ASDEX-U specifications
a_v = 0.8									# Vertical minor radius
a_h = 0.5									# Horizontal minor radius
a_m = ( (a_v**2 + a_h**2) / 2.0 )**(1.0/2.0)# Mean minor radius
R = 1.6										# Major radius
I_p = 2.0e6									# Plasma current
B_phi = 3.9									# Toroidal field
B_theta = mu_0 * I_p / ( 2*pi*a_m )			# Poloidal field
B = ( B_phi**2 + B_theta**2 )**(1.0/2.0) 	# Full field


## ITER specifications
#a_m = 2.0									# Mean minor radius
#R = 6.2										# Major radius
#I_p = 15.0e6								# Plasma current
#B_phi = 5.3									# Toroidal field
#B_theta = mu_0 * I_p / ( 2*pi*a_m )			# Poloidal field
#B = ( B_phi**2 + B_theta**2 )**(1.0/2.0) 	# Full field


aspect = a_m / R							# Aspect Ratio
q = aspect * B_phi/B_theta					# q value

## PRESET parameters for quick calculation, many of which
## are chosen by Staps and Paquay

if config.original_model == True:
	L = 4.0									# in AU
	lambda_n = 5.0/4.0						# Length scales for decay at edge
	lambda_T = 3.0/2.0
	lambda_Z = 5.0/4.0
elif config.original_model == False:
	L = 0.04								# in m
	lambda_n = 0.001						# Length scales for decay at edge
	lambda_T = 0.002
	lambda_Z = 0.001


D_max = 5.0
D_min = 1.0

epsilon = 1.0 / 25.0
mu = 1.0 / 20.0

## Choose set of parameters
# If numerical_choice is not defined, not a string, or set to Paquay:
if config.numerical_choice.lower() == "paquay":
	# Paquay's numbers
	zeta = 1.1
	c_n = 1.1
	c_T = 0.9
	a = -1.5
	b = 1.0
	c = -1.0
	Z_S = 1.4

# Stap's numbers
elif config.numerical_choice.lower() == "staps":
	zeta = 0.5
	c_n = -1.1
	c_T = -0.9
	a = 3.0/2.0
	b = 2.0
	c = -1.0
	Z_S = -3.0/2.0

# Heat diffusivity coefficient choice for Gradient Model
elif config.numerical_choice.lower() == "g_grad"\
		or config.numerical_choice.lower() == "gradient_model":
	zeta = 0.5

print "The numerical parameters are chosen to " + str(config.numerical_choice)

## For use in full flux model
alpha_an = 5.0				# Anomalous loss coefficient
a_in0 = 0.05				# Neutrals coefficient
alpha_cx = 0.9				# Charge exchange coefficient

