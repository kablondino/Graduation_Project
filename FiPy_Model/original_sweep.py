from fipy import TransientTerm, DiffusionTerm, Viewer, TSVViewer

from fipy.solvers import *

# variable_decl.py, parameters.py, and fipy.tools.numerix are included in the following
from coeffs import *

import os		# For saving files to a specified directory

# ------------- Check time and plot title config ----------
if (type(getattr(config, 'total_time', None)) != int or\
		getattr(config, 'total_time', None) <= 0):
	try:
		config.total_time = int(input("Total number of time steps not properly defined. Enter integer value: "))
		raw_input("Total time set to "+str(config.total_time))
	except:
		raw_input("Total time steps defaulted to 100.")
		config.total_time = 100

if (type(getattr(config, 'timeStep_denom', None)) != float or\
		getattr(config, 'timeStep_denom', None) <= 0.0):
	try:
		config.timeStep_denom = float(input("The denomintor of the time step size is not properly defined. Enter floating-point value: "))
		raw_input("The denominator of the time step size is set to "+str(config.timeStep_denom))
	except:
		config.timeStep_denom = 15.0
		raw_input("The denominator of the time step size is defaulted to 15.0")

if getattr(config, 'total_time', None) == None:
	config.total_time = input("Total number of time steps not properly defined. Enter integer value: ")
if getattr(config, 'timeStep_denom', None) == None:
	config.timeStep_denom = input("The denomintor of the time step size is not properly defined. Enter floating-point value: ")

if getattr(config, 'plot_title', None) == None:
	try:
		config.plot_title = str(input("Set the title of the plots: "))
	except:
		config.plot_title = ""


# ------------- Check file-writing configuration ----------
# Assumes save_directory exists, but not written correctly as a string
if hasattr(config, 'save_directory'):
	if type(config.save_directory) != str:
		config.save_directory = input("Directory for saving data is incorrectly written. Enter string wrapped in quotes: ")

# Assumes save_plots/TSV exist, but incorrectly input-ed
if hasattr(config, 'save_plots'):
	if type(config.save_plots) != bool:
		config.save_plots = input("Save plot images? Enter boolean: ")
if hasattr(config, 'save_TSVs'):
	if type(config.save_TSVs) != bool:
		config.save_TSVs = input("Save numerical data? Enter boolean: ")


# ------------------------- Boundary Conditions ----------------------------- #
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
# Available: LinearPCGSolver (Default), LinearGMRESSOlver, LinearLUSolver,
# LinearJORSolver	<-- Not working exactly
GMRES_Solver = LinearGMRESSolver(iterations=1000, tolerance=1.0e-6)

## Options for viewers
x_min, x_max, y_min, y_max = 0.0, config.L, 0.0, 5.0

#initial_viewer = Viewer((density, temperature, Z, Diffusivity),\
#		xmin=x_min, xmax=x_max, datamin=y_min, legend='best')
#raw_input("Pause for Initial Conditions")

timeStep = epsilon / config.timeStep_denom

res_tol = 1.0e-6

if __name__ == '__main__':

	# Initialize viewer
	viewer = Viewer((density, temperature, Z, Diffusivity),\
			xmin=x_min, xmax=x_max, legend='best',\
			title = config.plot_title)

	# File writing
	if hasattr(config, 'save_directory'):
		if not os.path.exists(os.getcwd() +str("/")+ config.save_directory):
			os.makedirs(os.getcwd() +str("/")+ config.save_directory)
			raw_input("Directory created: "+str(config.save_directory))

	for t in range(config.total_time):
		# (Re)set residual value(s)
		res_D = res_n = res_T = res_Z = res_full = 1.0e10

		# Update values
		Diffusivity.setValue(D_choice_local)
		density.updateOld(); temperature.updateOld(); Z.updateOld()

		# Solve the fully coupled equation
		while res_full > res_tol:
			print t, res_full
			res_full = full_equation.sweep(dt=timeStep, solver=GMRES_Solver)

		# Plot solution and save, if option is True
		if getattr(config, 'save_plots', False) == True:
			viewer.plot(filename =\
					config.save_directory+"/"+str(t).zfill(4)+".png")
		else:
			viewer.plot()

		# Save TSV's
		if getattr(config, 'save_TSVs', False) == True:
			TSVViewer(vars=(density, temperature, Z, Diffusivity)).plot(\
					filename=config.save_directory+"/"+str(t).zfill(4)+".tsv")


	raw_input(" =================================== End of Program. Press any key to continue. ==================================== ")

