#!/usr/bin/env python
"""
	This file sets the equation system and does the solving
	and viewing of the original model.
"""

from fipy import TransientTerm, DiffusionTerm, Viewer, TSVViewer
from fipy.solvers import *

# Order of file imports from the following inport: input_handling.py,
# parameters.py, variable_decl.py boundary_init_cond
from boundary_init_cond import *
from calculate_coeffs import *
# fipy.tools.numerix and dump is also imported from the above

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
G = a + b*(Z - Z_S) + c*(Z - Z_S)**3
S_Z = ((c_n*temperature) / density**2) * density.grad[0]\
		+ (c_T / density) * temperature.grad[0] + G
Z.equation = TransientTerm(coeff=epsilon, var=Z)\
		== DiffusionTerm(coeff=mu, var=Z) + S_Z

# Fully-Coupled Equation
full_equation = density.equation & temperature.equation & Z.equation


# Initialize all the coefficients and other variables
#calculate_coeffs()

# ----------------- Choose Solver -------------------------
# Available: LinearPCGSolver (Default), LinearGMRESSolver, LinearLUSolver,
# LinearJORSolver	<-- Not working exactly
PCG_Solver = LinearPCGSolver(iterations=100, tolerance=1.0e-6)
GMRES_Solver = LinearGMRESSolver(iterations=100, tolerance=1.0e-6)
LLU_Solver = LinearLUSolver(iterations=100, tolerance=1.0e-6)


# LOAD pickled H--Mode data
if (__name__ == '__main__' and config.initial_H_mode == True\
		and config.original_model == True):
	H_mode_data = dump.read("./L_start_pickle/state0090.dat")
	density.setValue(H_mode_data['density'])
	temperature.setValue(H_mode_data['temperature'])
	Z.setValue(H_mode_data['Z'])
	Diffusivity.setValue(D_choice_local)




timeStep = epsilon / config.timeStep_denom


if __name__ == '__main__':

	# Declare viewer
	viewer = Viewer((density, temperature, -Z, Diffusivity),\
			xmin=0.0, xmax=L,\
			datamax=config.ploty_max, legend='best',\
			title = config.plot_title)
	if config.show_initial == True:
		raw_input("Pause for Initial Conditions")

	# Auxiliary viewers
	if config.aux_plots == True:
		auxiliary1_viewer = Viewer((Gamma_an), xmin=0.0,\
				xmax=L, datamin=config.aux1y_min,\
				datamax=config.aux1y_max, legend='best',\
				title = config.aux_title1)
		auxiliary2_viewer = Viewer((Gamma_bulk), xmin=0.0,\
				xmax=L, datamin=config.aux2y_min,\
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


	# Set the tolerance for the residual
	res_tol = 1.0e-5


	# ----------------- Time Loop -------------------------
	for t in range(config.total_timeSteps):
		# (Re)set residual values
		current_residual = 1.0e100

		# Update values
		density.updateOld(); temperature.updateOld(); Z.updateOld()
		calculate_coeffs(); Diffusivity.setValue(D_choice_local)


		# --------------- Solving Loop --------------------
		while current_residual > res_tol:
			print t, current_residual
			current_residual = full_equation.sweep(dt=timeStep,\
					solver=GMRES_Solver)


		# Plot solution and save, if option is True
		if config.save_plots == True:
			viewer.plot(filename =\
					config.save_directory+"/"+str(t).zfill(4)+".png")

			# Save auxiliary plots
			if config.aux_plots == True:
				auxiliary1_viewer.plot(filename =\
						config.save_directory+"/aux1_"+str(t).zfill(4)+".png")
				auxiliary2_viewer.plot(filename =\
						config.save_directory+"/aux2_"+str(t).zfill(4)+".png")

		elif config.save_plots == False:
			viewer.plot()

			if config.aux_plots == True:
				auxiliary1_viewer.plot(); auxiliary2_viewer.plot()

		# Save TSV's
		if config.save_TSVs == True:
			TSVViewer(vars=\
					(Gamma_an, Gamma_cx, Gamma_bulk, Gamma_OL,\
					Z_transient_coeff, Z_diffusion_coeff, density.grad[0],\
					temperature.grad[0])).plot(filename=\
					config.save_directory+"/"+str(t).zfill(4)+".tsv")


	raw_input(" <=============== End of Program. Press any key to continue. ===============> ")

