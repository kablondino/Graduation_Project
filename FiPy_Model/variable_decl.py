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
from fipy.tools import numerix

from parameters import *


# -------------- Check Configuration Variables ------------
# If saving data is enabled, but not a directory, exit the run.
# It is done early as to not waste time.
if (getattr(config, 'save_directory', None) == None and\
		(getattr(config, 'save_plots', False) == True or\
		getattr(config, 'save_TSVs', False) == True)):
	sys.exit("No directory specified for saving specified files. Exiting...")

if (type(getattr(config, 'nx', None)) != int or\
		getattr(config, 'nx', None) <= 0):
	try:
		config.nx = int(input("nx (Grid number) not properly defined. Enter a positive integer value: "))
		raw_input("nx set to "+str(config.nx))
	except:
		config.nx = 100
		raw_input("nx defaulted to 100.")

if (type(getattr(config, 'L', None)) != float or\
		getattr(config, 'L', None) <= 0.0):
	try:
		config.L = float(input("Length of domain not properly defined. Enter floating-point value: "))
		raw_input("L set to "+str(config.L))
	except:
		config.L = 4.0
		raw_input("L defaulted to 4.0.")


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
if getattr(config, 'initial_H_mode', False) == False:
	config.initial_H_mode = False
	Z.setValue(0.0) # L--mode
	if not hasattr(config, 'initial_H_mode'):
		raw_input("Initial working mode is defaulted to L--mode")
		config.initial_H_mode = False

elif getattr(config, 'initial_H_mode', False) == True:
	Z.setValue(Z_S*(1.0 - numerix.tanh(\
			(config.L*x - config.L) / 2.0))) # H--mode

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
if (type(getattr(config, 'D_choice', None)) != str or\
		getattr(config, 'D_choice', None).lower() not in\	# Not quite working
		["d_zohm", "d_shear"] or getattr(config, 'D_choice', "").lower()\
		== "d_staps"):
	D_choice_local = D_Staps
	config.D_choice = "D_Staps"
	raw_input("Diffusivity model defaulted to Stap's.")

if config.D_choice.lower() == "d_staps":
	D_choice_local = D_Staps
elif config.D_choice.lower() == "d_zohm":
	D_choice_local = D_Zohm
elif config.D_choice.lower() == "d_shear":
	D_choice_local = D_Shear

else:
	raw_input("Something went horribly wrong in choosing the Diffusivity model.")

Diffusivity.setValue(D_choice_local)

# --------- Set Initial Conditions for n and T ------------
# Initial conditions for L--mode
if config.initial_H_mode == False:
	density.setValue(-(Gamma_c*lambda_n / Diffusivity) * (1.0 + x/lambda_n))

	temperature.setValue(q_c*((gamma - 1.0) / Gamma_c) *(1.0 - lambda_n /\
			(zeta*lambda_T + lambda_n)*(1.0 + x/lambda_n)**(-zeta)))

# Initial conditions for H--mode
elif config.initial_H_mode == True:
	density.setValue(-(Gamma_c*lambda_n / Diffusivity) * (1.0 + x/lambda_n))

	temperature.setValue(q_c*((gamma - 1.0) / Gamma_c) *(1.0 - lambda_n /\
			(zeta*lambda_T + lambda_n)*(1.0 + x/lambda_n)**(-zeta)))

else:
	print "Something went horribly wrong in choosing initial conditions."

U.setValue(density*temperature / (gamma - 1.0))

