"""
	This file contains the simple constant delcarations
	for use in the PDE system. This includes machine
	parameters, length scales, global physical/math
	constants, and non-g model numbers.

"""

from input_handling import *

from fipy.tools.dimensions.physicalField import *

# ----------------- Constant Parameters -------------------
## Global parameters and physical constants
pi = 3.141592653589793
charge = PhysicalField(value=1.0, unit="e")	# Dummy charge
k_B_J = 1.38064852e-23			# Boltzmann in J/K
k_B_eV = k_B_J / charge			# Boltzmann in eV/K
m_e = PhysicalField(value=1.0, unit="me")	# 9.10938356e-31, Electron mass
m_i = PhysicalField(value=1.0, unit="mp")	# 1.67262171e-27, Ion (H) mass
epsilon_0 = PhysicalField(value=1.0, unit="eps0")	# = 8.854187817e-12, Permittivity of free space
mu_0 = PhysicalField(value=1.0, unit="mu0")			# = 4*pi*1.0e-7, Permeability of free space
c = 1 / (epsilon_0*mu_0)**(1.0/2.0)

## ASDEX-U specifications
# Vertical minor radius
a_v = PhysicalField(value=0.8, unit="m")
# Horizontal minor radius
a_h = PhysicalField(value=0.5, unit="m")
# Mean minor radius
a_m = PhysicalField(value=((a_v**2 + a_h**2) / 2.0)**(1.0/2.0), unit="m")
# Major radius
R = PhysicalField(value=1.6, unit="m")
# Aspect Ratio
aspect = PhysicalField(value=a_m / R, unit="")
# Plasma current
I_p = PhysicalField(value=2.0e6, unit="A")
# Toroidal field
B_phi = PhysicalField(value=3.9, unit="T")
# Poloidal field
B_theta = PhysicalField(value=mu_0 * I_p / ( 2*pi*a_m ), unit="T")
# q value
q = PhysicalField(value=aspect * B_phi/B_theta, unit="")
# Full field
B = PhysicalField(value=( B_phi**2 + B_theta**2 )**(1.0/2.0) , unit="T")

## PRESET parameters for quick calculation, many of which
## are chosen by Staps and Paquay
gamma = 5.0/3.0

## Length scales for decay at edge
lambda_n = PhysicalField(value=5.0/4.0, unit="m")
lambda_T = PhysicalField(value=3.0/2.0, unit="m")
lambda_Z = PhysicalField(value=5.0/4.0, unit="m")

D_max = PhysicalField(value=3.0, unit="m**2/s")
D_min = PhysicalField(value=0.2, unit="m**2/s")

epsilon = PhysicalField(value=1.0 / 25.0, unit="s")
mu = PhysicalField(1.0 / 20.0, unit="m**2/s")

## Choose set of parameters
# If numerical_choice is not defined, not a string, or set to Paquay:
if config.numerical_choice.lower() == "paquay":
	# Paquay's numbers
	zeta = 1.1
	c_n = PhysicalField(value=1.1, unit="1.0e19 / (eV*m**3)")
	c_T = PhysicalField(value=0.9, unit="1.0e19 / m**3")
	a = -1.5
	b = 1.0
	c = -1.0
	Z_S = 1.4

# Stap's numbers
elif config.numerical_choice.lower() == "staps":
	zeta = 0.5
	c_n = PhysicalField(value=-1.1, unit="1.0e19 / (eV*m**3)")
	c_T = PhysicalField(value=-0.9, unit="1.0e19 / m**3")
	a = 3.0/2.0
	b = 2.0
	c = -1.0
	Z_S = -3.0/2.0

# Heat diffusivity coefficient choice for Gradient Model
elif config.numerical_choice.lower() == "g_grad"\
		or config.numerical_choice.lower() == "gradient_model":
	zeta = 0.9

print "The numerical parameters are chosen to " + str(config.numerical_choice)

## For use in full flux model
alpha_an = 1.0				# Anomalous loss coefficient
a_in0 = 0.05				# Neutrals coefficient
neu_react_rate = PhysicalField(value=1.0e-6, unit="m**3/s")	# Reaction rate of charge exchange
alpha_cx = 0.9				# Charge exchange coefficient

