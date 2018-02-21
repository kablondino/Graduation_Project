"""
	This file generates the 1D mesh and the 4 cell
	variables needed for the model (including Diffusivity.
	It also sets the initial conditions for each variable.

	A job configuration file is REQUIRED in order to
	properly choose the initial operating mode, number of
	cells, etc.
"""
from fipy import Grid1D, CellVariable, FaceVariable
from fipy.tools import numerix

from parameters import *

# ----------------- Mesh Generation -----------------------
L = 5.0
mesh = Grid1D(nx=job_config.nx, Lx=L)

x = mesh.cellCenters[0] # Cell position
X = mesh.faceCenters[0] # Face position, if needed

# ----------------- Variable Declarations -----------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True)

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True)

U = CellVariable(name=r"$U$", mesh=mesh, hasOld=True)

Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True)

Diffusivity = CellVariable(name=r"$D$", mesh=mesh, hasOld=True)

# ----------- Initial Conditions of Z ---------------------
if job_config.Initial_H_mode == True:
	Z.setValue(Z_S*(1.0 - numerix.tanh((L*x - L) / 2.0))) # H--mode
elif job_config.Initial_H_mode == False:
	Z.setValue(0.0) # L--mode
else:
	sys.exit("Initial working mode not properly selected")

# ----------------- Diffusivities -------------------------
# Itohs'/Zohm's model
D_Zohm = (D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(Z)) / 2.0
# Stap's Model
alpha_sup = 0.5
D_Staps = D_min + (D_max - D_min) / (1.0 + alpha_sup*numerix.dot(Z.grad, Z.grad))
# Flow-Shear Model
a1, a3 = 1.0, 0.5	# ASSUMES a2 = 0
D_Shear = D_min + (D_max - D_min) / (1.0 + a1*(Z)**2 + a3*numerix.dot(Z.grad, Z.grad))

if job_config.D_choice.lower() == "d_zohm":
	D_choice_local = D_Zohm
elif job_config.D_choice.lower() == "d_staps":
	D_choice_local = D_Staps
elif job_config.D_choice.lower() == "d_shear":
	D_choice_local = D_Shear
else:
	sys.exit("Error in choosing Diffusivity model")

Diffusivity.setValue(D_choice_local)

# --------- Set Initial Conditions for n and T ------------
# Initial conditions for H--mode
if job_config.Initial_H_mode == True:
	density0H = -(Gamma_c*lambda_n / Diffusivity) * (1.0 + x/lambda_n)
	density.setValue(density0H)

	temp0H = q_c*((gamma - 1.0) / Gamma_c) *(1.0 - lambda_n /\
			(zeta*lambda_T + lambda_n)*(1.0 + x/lambda_n)**(-zeta))
	temperature.setValue(temp0H)

# Initial conditions for L--mode
elif job_config.Initial_H_mode == False:
	density0L = -(Gamma_c*lambda_n / Diffusivity) * (1.0 + x/lambda_n)
	density.setValue(density0L)

	temp0L = q_c*((gamma - 1.0) / Gamma_c) * (1.0 - lambda_n /\
			(zeta*lambda_T + lambda_n)*(1.0 + x/lambda_n)**(-zeta))
	temperature.setValue(temp0L)

else:
	sys.exit("Initial working mode not properly selected")

U.setValue(density*temperature / (gamma - 1.0))

print "Success!"

