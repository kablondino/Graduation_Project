from fipy import TransientTerm, DiffusionTerm, Viewer, TSVViewer

from fipy.solvers import *

# input_handling, parameters.py, variable_decl.py, (in that order)
# and fipy.tools.numerix are included in the following import
from coeffs import *

import os	# For saving files to a specified directory


# ----------------- Boundary Conditions ------------------- #
def set_boundary_values(AGamma_c, Aq_c):
	"""
		Density Boundary Conditions:
		d/dx(n(0)) == n / lambda_n
		d/dx(n(L)) == -Gamma_c / Diffusivity
	"""
	density.faceGrad.constrain(density.faceValue / lambda_n, mesh.facesLeft)
	density.faceGrad.constrain(\
			-AGamma_c / Diffusivity.faceValue, mesh.facesRight)

	"""
		Temperature Boundary Conditions:
		d/dx(T(0)) = T / lambda_T
		d/dx(T(L)) = zeta*(Gamma_c*T - q_c*(gamma - 1)) / (Diffusivity * n)
	"""
	temp_left = temperature.faceValue / lambda_T
	temperature.faceGrad.constrain(temp_left, mesh.facesLeft)
	temp_right =\
	temperature.faceGrad.constrain((zeta * (AGamma_c*temperature.faceValue -\
			Aq_c*(gamma - 1.0))) / (Diffusivity.faceValue *\
			density.faceValue), mesh.facesRight)

	"""
		Paquay considered these Z Boundary Conditions:
		d/dx(Z(0)) == Z / lambda_Z
		mu*D/epsilon * d/dx(Z(L)) == 0
	"""
	Z.faceGrad.constrain(Z.faceValue / lambda_Z, mesh.facesLeft)
	Z.faceGrad.constrain(0.0, mesh.facesRight)

set_boundary_values(Gamma_c, q_c)

# ----------------- PDE Declarations ----------------------
# Density Equation
density.equation = TransientTerm(coeff=1.0, var=density)\
		== DiffusionTerm(coeff=Diffusivity, var=density)

# Energy Equation
temperature.equation = TransientTerm(coeff=density, var=temperature)\
		== DiffusionTerm(coeff=(Diffusivity*density/zeta), var=temperature)\
		+ DiffusionTerm(coeff=Diffusivity*temperature, var=density)

# Z Equation
G = a + b*(Z - Z_S) + c*(Z - Z_S)**3
S_Z = ((c_n*temperature) / density**2) * density.grad[0]\
		+ (c_T / density) * temperature.grad[0] + G
Z.equation = TransientTerm(coeff=epsilon, var=Z)\
		== DiffusionTerm(coeff=mu, var=Z) + S_Z

# Fully-Coupled Equation
full_equation = density.equation & temperature.equation & Z.equation

# ----------------- Choose Solver -------------------------
# Available: LinearPCGSolver (Default), LinearGMRESSolver, LinearLUSolver,
# LinearJORSolver	<-- Not working exactly
PCG_Solver = LinearPCGSolver(iterations=100, tolerance=1.0e-6)
GMRES_Solver = LinearGMRESSolver(iterations=100, tolerance=1.0e-6)
LLU_Solver = LinearLUSolver(iterations=100, tolerance=1.0e-6)

# Initial conditions viewer
#initial_viewer = Viewer((density, temperature, Z, Diffusivity),\
#		xmin=0.0, xmax=config.plotx_max, legend='best')
#raw_input("Pause for Initial Conditions")

timeStep = mu / config.timeStep_denom

if __name__ == '__main__':

	# Initialize viewer
	viewer = Viewer((density, temperature, -Z, Diffusivity),\
			xmin=0.0, xmax=config.plotx_max,\
			datamax=config.ploty_max, legend='best',\
			title = config.plot_title)

	# File writing
	if (hasattr(config, 'save_directory') and\
			(getattr(config, 'save_plots', False) == True or\
			getattr(config, 'save_TSVs', False) == True)):
		if not os.path.exists(os.getcwd() +str("/")+ config.save_directory):
			os.makedirs(os.getcwd() +str("/")+ config.save_directory)
			raw_input("Directory created: " +str(config.save_directory))

	for t in range(config.total_timeSteps):
		# (Re)set residual value
		res_full = 1.0e10

		# Update values
		Diffusivity.setValue(D_choice_local)
		density.updateOld(); temperature.updateOld(); Z.updateOld()

		# Solve the fully coupled equation
		while res_full > config.res_tol:
			print t, res_full
			res_full = full_equation.sweep(dt=timeStep, solver=GMRES_Solver)

		# Plot solution and save, if option is True
		if config.save_plots == True:
			viewer.plot(filename =\
					config.save_directory+"/"+str(t).zfill(4)+".png")
		else:
			viewer.plot()

		# Save TSV's
		if config.save_TSVs == True:
			TSVViewer(vars=(density, temperature, Z, Diffusivity)).plot(\
					filename=config.save_directory+"/"+str(t).zfill(4)+".tsv")


	raw_input(" <================================== End of Program. Press any key to continue. ===================================> ")

