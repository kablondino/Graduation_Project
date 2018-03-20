"""
	This file generates the 1D mesh and the 4 cell
	variables needed for the model (including Diffusivity.
	It also sets the initial conditions for each variable.

	A job configuration file is REQUIRED in order to
	properly choose the initial operating mode, number of
	cells, etc.

	In addition, it checks the configuration file's
	variable types to avoid errors.
"""
from fipy import Grid1D, CellVariable, FaceVariable
from fipy.tools import numerix, dump

from parameters import *


# ----------------- Mesh Generation -----------------------
mesh = Grid1D(nx=config.nx, Lx=config.L)

x = mesh.cellCenters[0] # Cell position
X = mesh.faceCenters[0] # Face position, if needed

# ----------------- Variable Declarations -----------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True)

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True)

U = CellVariable(name=r"$U$", mesh=mesh, hasOld=True)

Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True)

Diffusivity = CellVariable(name=r"$D$", mesh=mesh, hasOld=True)

# ----------- Initial Conditions of Z ---------------------
if config.initial_H_mode == False:
	Z.setValue(0.0) # L--mode
elif config.initial_H_mode == True:
	Z.setValue(Z_S*(1.0 - numerix.tanh(\
			(config.L*x - config.L) / 2.0))) # H--mode
else:
	print "The check for initial_H_mode is not working.... stop."

# ----------------- Diffusivities -------------------------
# Itohs'/Zohm's model
D_Zohm = (D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(Z)) / 2.0
# Stap's Model
alpha_sup = 0.5
D_Staps = D_min + (D_max - D_min) / (1.0 + alpha_sup\
		*numerix.dot(Z.grad, Z.grad))
# Flow-Shear Model
a1, a3 = 1.0, 0.5	# ASSUMES a2 = 0
D_Shear = D_min + (D_max - D_min) / (1.0 + a1*(Z)**2 +\
		a3*numerix.dot(Z.grad, Z.grad))

# Set Diffusivity. It defaults to Stap's version
if config.D_choice.lower() == "d_zohm":
	D_choice_local = D_Zohm
elif config.D_choice.lower() == "d_staps":
	D_choice_local = D_Staps
elif config.D_choice.lower() == "d_shear" or\
		config.D_choice.lower() == "d_flow_shear":
	D_choice_local = D_Shear

else:
	print "Something went horribly wrong in choosing the Diffusivity model."

Diffusivity.setValue(D_choice_local)
print "The diffusivity model is set to " + str(config.D_choice)

# --------- Set Initial Conditions for n and T ------------
if config.SI_units == True:					# NEEDS TO BE FIXED!
	density_si_coeff = 0.05		# x 10^20 m^-3
	temp_si_coeff = 0.5			# eV
elif config.SI_units == False:
	density_si_coeff = 1.0		# m^-3
	temp_si_coeff = 1.0			# eV
else:
	print "Something went horribly wrong when choosing the units."

# Initial conditions for L--mode
if config.initial_H_mode == False:
	density.setValue(-(density_si_coeff*config.Gamma_c*lambda_n / Diffusivity)\
			* (1.0 + x/lambda_n))

	temperature.setValue(temp_si_coeff*config.q_c\
			*((gamma - 1.0) / config.Gamma_c)\
			*(1.0 - lambda_n / (zeta*lambda_T + lambda_n)\
			*(1.0 + x/lambda_n)**(-zeta)))

# Initial conditions for H--mode
elif config.initial_H_mode == True:
	density.setValue(-(density_si_coeff*config.Gamma_c*lambda_n\
			/ Diffusivity) * (1.0 + x/lambda_n))

	temperature.setValue(temp_si_coeff*config.q_c*((gamma - 1.0)\
			/ config.Gamma_c) * (1.0 - lambda_n / (zeta*lambda_T + lambda_n)\
			*(1.0 + x/lambda_n)**(-zeta)))

else:
	print "Something went horribly wrong in choosing initial conditions."

U.setValue(density*temperature / (gamma - 1.0))

