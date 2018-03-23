#!/usr/bin/env python
from fipy import TransientTerm, DiffusionTerm, Viewer, TSVViewer
from fipy.solvers import *

# Order of file imports from the following inport: input_handling.py,
# parameters.py, variable_decl.py boundary_init_cond
from calculate_coeffs import *
# fipy.tools.numerix is also imported from the above

import os	# For saving files to a specified directory


# ----------------- PDE Declarations ----------------------
# Density Equation
density.equation = TransientTerm(coeff=1.0, var=density)\
		== DiffusionTerm(coeff=Diffusivity, var=density)

# Energy Equation
temperature.equation = TransientTerm(coeff=density, var=temperature)\
		== DiffusionTerm(coeff=(Diffusivity*density/zeta), var=temperature)\
		+ DiffusionTerm(coeff=Diffusivity*temperature, var=density)

# Z Equation
if config.original_model == False:
	# Full Flux model
	########## WILL BE DEFINED SOON ##########
	print "Flux model chosen."
	G = a + b*(Z - Z_S) + c*(Z - Z_S)**3
	S_Z = ((c_n*temperature) / density**2) * density.grad[0]\
			+ (c_T / density) * temperature.grad[0] + G
	Z.equation = TransientTerm(coeff=epsilon, var=Z)\
			== DiffusionTerm(coeff=mu, var=Z) + S_Z
else:
	# Original numerical, arbitrary units model
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
initial_viewer = Viewer((density, temperature, -Z, Diffusivity),\
		xmin=0.0, xmax=config.plotx_max, legend='best')
raw_input("Pause for Initial Conditions")

timeStep = epsilon / config.timeStep_denom

if __name__ == '__main__':

	# Initialize viewer
	viewer = Viewer((density, temperature, -Z, Diffusivity),\
			xmin=0.0, xmax=config.plotx_max,\
			datamax=config.ploty_max, legend='best',\
			title = config.plot_title)

	# Auxiliary viewers
	if config.aux_plots == True:
		auxiliary1_viewer = Viewer((omega_bi), xmin=0.0,\
				xmax=config.plotx_max, datamin=config.aux1y_min,\
				datamax=config.aux1y_max, legend='best',\
				title = config.aux_title1)
		auxiliary2_viewer = Viewer((omega_be), xmin=0.0,\
				xmax=config.plotx_max, datamin=config.aux2y_min,\
				datamax=config.aux2y_max, legend='best',\
				title = config.aux_title2)

	# File writing
	if (hasattr(config, 'save_directory') and\
			(getattr(config, 'save_plots', False) == True or\
			getattr(config, 'save_TSVs', False) == True)):
		if not os.path.exists(os.getcwd() +str("/")+ config.save_directory):
			os.makedirs(os.getcwd() +str("/")+ config.save_directory)
			raw_input("Directory created: " +str(config.save_directory))
		raw_input("Pause set for writing to file...")

	for t in range(config.total_timeSteps):
		# (Re)set residual value
		res_full = 1.0e100

		# Update values
		Diffusivity.setValue(D_choice_local)
		density.updateOld(); temperature.updateOld(); Z.updateOld()
		update_g_coeffs()

		# Solve the fully coupled equation
		while res_full > config.res_tol:
			print t, res_full
			res_full = full_equation.sweep(dt=timeStep, solver=GMRES_Solver)

		# Plot solution and save, if option is True
		if config.save_plots == True:
			viewer.plot(filename =\
					config.save_directory+"/"+str(t).zfill(4)+".png")
			if config.aux_plots == True:
				auxiliary1_viewer.plot(filename =\
						config.save_directory+"/aux1_"+str(t).zfill(4)+".png")
				auxiliary2_viewer.plot(filename =\
						config.save_directory+"/aux2_"+str(t).zfill(4)+".png")
		else:
			viewer.plot()
			if config.aux_plots == True:
				auxiliary1_viewer.plot(); auxiliary2_viewer.plot()

		# Save TSV's
		if config.save_TSVs == True:
			all_variables = (density, temperature, Z, Diffusivity, v_Ti,\
					v_Te, rho_pi, rho_pe, omega_bi, omega_be, nu_ei, nu_ii,\
					nu_in0, nu_eff, nu_ai, nu_ae, Gamma_an, Gamma_cx,\
					Gamma_bulk, Gamma_OL)
			TSVViewer(vars=all_variables).plot(\
					filename=config.save_directory+"/"+str(t).zfill(4)+".tsv")


	raw_input(" <=============== End of Program. Press any key to continue. ===============> ")

