"""
	This file contains the simple numerical delcarations
	for use in the PDE system. This includes machine
	parameters, length scales, global physical/math
	constants, and non-g model numbers.

"""
from input_handling import *


# ----------------- Parameters ----------------------------
## Global parameters and physical constants
pi = 3.141592653589793
charge = 1.0					# Dummy charge
charge_true = 1.602e-19			# Elementary charge
k_B = 8.617e-5					# Boltzmann in eV
m_e = 9.109e-31					# Electron mass
m_i = 1.673e-27					# Ion (H) mass
epsilon_0 = 8.854187817e-12		# Permittivity of free space
mu_0 = 4*pi*1.0e-7				# Permeability of free space


## ASDEX-U specifications
a_v = 0.8									# Vertical minor radius
a_h = 0.5									# Horizontal minor radius
a_m = ( (a_v**2 + a_h**2) / 2.0 )**(1.0/2.0)# Mean minor radius
R = 1.6										# Major radius
aspect = a_m / R							# Aspect Ratio
I_p = 2.0e6									# Plasma current
B_phi = 3.9									# Toroidal field
B_theta = mu_0 * I_p / ( 2*pi*a_m )			# Poloidal field
q = aspect * B_phi/B_theta					# q value
B = ( B_phi**2 + B_theta**2 )**(1.0/2.0) 	# Full field


## PRESET parameters for quick calculation, many of which
## are chosen by Staps and Paquay
Gamma_c = -4.0/5.0
q_c = -4.0
gamma = 5.0/3.0

lambda_n = 5.0/4.0		# Length scales for decay at edge
lambda_T = 3.0/2.0
lambda_Z = 5.0/4.0

D_max = 2.0
D_min = 2.0/5.0

epsilon = 1.0 / 25.0
mu = 1.0/20.0

## Choose set of parameters
# If numerical_choice is not defined, not a string, or set to Paquay:
if config.numerical_choice.lower() == "paquay":
	# Paquay's numbers
	zeta = 1.1
	c_n = 1.1
	c_T = 0.9
	a = -1.5
	b = -1.0
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
	zeta = 0.9

print "The numerical parameters are chosen to " + str(config.numerical_choice)


## For use in gradient model
alpha_cx = 0.9				# Charge exchange coefficient
alpha_an = 1.0				# Anomalous loss coefficient
a_in0 = 0.1					# Neutrals coefficient
neu_react_rate = 1.0e-14	# Reaction rate of charge exchange

