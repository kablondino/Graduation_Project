"""
	This file contains the simple numerical delcarations
	for use in the PDE system. This includes machine
	parameters, length scales, global physical/math
	constants, and non-g model numbers.
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
B_phi = 3.9
B_theta = mu_0 * I_p / ( 2*pi*a_m )
q = aspect * B_phi/B_theta
B = ( B_phi**2 + B_theta**2 )**(1.0/2.0)

## PRESET parameters for quick calculation, many of which
## are chosen by Staps and Paquay
#zeta = 0.5		# Stap's
zeta = 1.1		# Paquay's
Gamma_c = -4.0/5.0
q_c = -4.0
gamma = 5.0/3.0

lambda_n = 5.0/4.0
lambda_T = 3.0/2.0
lambda_Z = 5.0/4.0

D_max = 2.0
D_min = 2.0/5.0

epsilon = 1.0/25.0
mu = 1.0/20.0

c_n = 1.1
c_T = 0.9

a = 3.0/2.0
b = 2.0
c = -1.0
Z_S = -3.0/2.0

