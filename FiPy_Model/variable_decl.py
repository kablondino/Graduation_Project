"""
	This file generates the 1D mesh and the 3 cell
	variables needed for the model, with optional
	face variable for particle diffusivity.
"""
from fipy import Grid1D, CellVariable, FaceVariable

from fipy.tools import numerix

from parameters import *

# ----------------- Mesh Generation -----------------------
nx = 100
L = 5.0
mesh = Grid1D(nx=nx, Lx=L)

x = mesh.cellCenters[0]	# Cell position
X = mesh.faceCenters[0] # Face position

# ----------------- Variable Declarations -----------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True)

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True)

Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True)

Diffusivity = CellVariable(name=r"$D$", mesh=mesh, hasOld=True)

# ----------- Initial Conditions of Z ---------------------
Z0L = 0 # L--mode
Z0H = Z_S*(1 - numerix.tanh((L*x - L) / 2.0)) # H--mode
Z.setValue(Z0H)

# ----------------- Diffusivities -------------------------
# Itohs'/Zohm's model
D_Zohm = (D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(Z)) / 2.0
# Stap's Model
alpha_sup = 0.5
D_Staps = D_min + (D_max - D_min) / (1.0 + alpha_sup*numerix.dot(Z.grad, Z.grad))
# Flow-Shear Model
a1, a3 = 1.0, 0.5	# ASSUMES a2 = 0
D_Shear = D_min + (D_max - D_min) / (1.0 + a1*(Z)**2 + a3*numerix.dot(Z.grad, Z.grad))

# CHOOSE DIFFUSIVITY HERE!
D_choice = D_Staps

# If Diffusivity is a Cell/Face variable
Diffusivity.setValue(D_choice)


density0 = CellVariable(name=r"$n_0$", mesh=mesh, value=-(Gamma_c*lambda_n / Diffusivity) * (1.0 + x/lambda_n))

temp0 = CellVariable(name=r"$T_0", mesh=mesh, value = q_c*((gamma - 1.0) / Gamma_c) * (1.0 - lambda_n / (zeta*lambda_T + lambda_n)*(1.0 + x/lambda_n)**(-zeta)))

