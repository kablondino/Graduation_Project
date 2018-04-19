#!/usr/bin/env python
"""
	This file sets the equation system and does the solving
	and viewing of the full flux model.
"""
# Order of file imports from the following inport: input_handling.py,
# parameters.py, variable_decl.py boundary_init_cond
from boundary_init_cond import *
from calculate_coeffs import *
# fipy.tools.numerix and dump is also imported from the above

if config.view_matplotlib == True:
	import matplotlib.pyplot as plt
	from fipy import MatplotlibViewer

from fipy import TransientTerm, DiffusionTerm, Viewer, TSVViewer
from fipy.solvers import *

import os	# For saving files to a specified directory


# Initialize all the coefficients and other variables
calculate_coeffs()

# ----------------- PDE Declarations ----------------------
# Density Equation
density.equation = TransientTerm(coeff=1.0, var=density)\
		== DiffusionTerm(coeff=Diffusivity, var=density)

# Energy Equation
temperature.equation = TransientTerm(coeff=density, var=temperature)\
		== DiffusionTerm(coeff=(Diffusivity*density/zeta), var=temperature)\
		+ DiffusionTerm(coeff=Diffusivity*temperature, var=density)

# Z Equation
Z.equation = TransientTerm(coeff=Z_transient_coeff, var=Z)\
		== DiffusionTerm(coeff=Z_diffusion_coeff, var=Z)\
		+ (Gamma_an - Gamma_cx  - Gamma_OL)
#Z.equation = TransientTerm(coeff=1.0, var=Z)\
#		== DiffusionTerm(coeff=1.0, var=Z)\
#		+ Flux_coeff * (Gamma_an - Gamma_cx - Gamma_bulk - Gamma_OL)

# Fully-Coupled Equation
full_equation = density.equation & temperature.equation & Z.equation


# LOAD pickled H--Mode data
#if __name__ == '__main__' and config.initial_H_mode == True:
#	H_mode_data = dump.read("./L_start_pickle/state0090.dat")
#	density.setValue(H_mode_data['density'])
#	temperature.setValue(H_mode_data['temperature'])
#	Z.setValue(H_mode_data['Z'])
#	Diffusivity.setValue(D_choice_local)


# ----------------- Choose Solver -------------------------
# Available: LinearPCGSolver (Default), LinearGMRESSolver, LinearLUSolver,
# LinearJORSolver	<-- Not working exactly
PCG_Solver = LinearPCGSolver(iterations=100, tolerance=1.0e-6)
GMRES_Solver = LinearGMRESSolver(iterations=100)
LLU_Solver = LinearLUSolver(iterations=100, tolerance=1.0e-6)


if __name__ == '__main__':

	# Declare viewer
	if config.view_matplotlib == True:
		state_fig = plt.figure()

		density_ax = plt.subplot((221))
		temp_ax = plt.subplot((222))
		Z_ax = plt.subplot((223))
		D_ax = plt.subplot((224))

		density_viewer = MatplotlibViewer(density, xmin=0.0, xmax=L,\
				datamin=0.0, datamax=3.0e19, legend='best', axes=density_ax)
		temp_viewer = MatplotlibViewer(temperature, xmin=0.0, xmax=L,\
				datamin=80.0, datamax=500.0, legend='best', axes=temp_ax)
		Z_viewer = MatplotlibViewer(Z, xmin=0.0, xmax=L,\
				legend='best', axes=Z_ax)
		D_viewer = MatplotlibViewer(Diffusivity, xmin=0.0, xmax=L, datamin=0.0,\
				datamax=D_max + D_max/10.0, legend='best', axes=D_ax)

		plt.suptitle(config.plot_title + diff_title)

	elif config.view_matplotlib == False:
		density_viewer = Viewer(density, xmin=0.0, xmax=L,\
				datamin=0.0, datamax=3.0e19, legend='best')
		temp_viewer = Viewer(temperature, xmin=0.0, xmax=L,\
				datamin=80.0, datamax=500.0, legend='best')
		Z_viewer = Viewer(Z, xmin=0.0, xmax=L,\
				legend='best')
		D_viewer = Viewer(Diffusivity, xmin=0.0, xmax=L, datamin=0.0,\
				datamax=D_max + D_max/10.0, legend='best')

	if config.show_initial == True:
		if config.view_matplotlib == True:
			state_fig.show()
		raw_input("Pause for Viewing Initial Conditions")

	# Auxiliary viewers
	if config.aux_plots == True:
		aux_plot_array = []
		for k in range(len(config.aux_vars)):
			aux_plot_array.append(Viewer(variable_dictionary\
					[config.aux_vars[k]], xmin=0.0, xmax=L,\
					datamin=config.aux_ymin[k], datamax=config.aux_ymax[k],\
					legend='best', title=config.aux_titles[k]))
		raw_input("Pause for Viewing Initial Auxiliary Plots")


	# File writing
	if (hasattr(config, 'save_directory') and\
			(getattr(config, 'save_plots', False) == True or\
			getattr(config, 'save_TSVs', False) == True)):
		if not os.path.exists(os.getcwd() +str("/")+ config.save_directory):
			os.makedirs(os.getcwd() +str("/")+ config.save_directory)
			raw_input("Directory created: " +str(config.save_directory))
		raw_input("Pause set for writing to file...")


	# ----------------- Time Loop -------------------------
	for t in range(config.total_timeSteps):
		# (Re)set residual values
		current_residual = 1.0e100

		# Update values
		density.updateOld(); temperature.updateOld(); Z.updateOld()
		Diffusivity.setValue(D_choice_local)
		calculate_coeffs()

		# --------------- Solving Loop --------------------
		while current_residual > config.res_tol:
			print t, current_residual
			current_residual = full_equation.sweep(dt=config.timeStep,\
					solver=GMRES_Solver)


		# Plot solution and save, if option is True
		if config.save_plots == True:
			density_viewer.plot(filename=config.save_directory + "/n"\
					+str(t).zfill(4)+ ".png")
			temp_viewer.plot(filename=config.save_directory + "/T" \
					+str(t).zfill(4)+ ".png")
			Z_viewer.plot(filename=config.save_directory + "/Z" \
					+str(t).zfill(4)+ ".png")
			D_viewer.plot(filename=config.save_directory + "/D" \
					+str(t).zfill(4)+ ".png")

			# Save auxiliary plots
			if config.aux_plots == True:
				for current_aux in range(len(aux_plot_array)):
					aux_plot_array[current_aux].plot(filename =\
							config.save_directory + "/aux" +str(current_aux)+\
							"_" +str(t).zfill(4)+ ".png")

		# If not set to save
		elif config.save_plots == False:
			density_viewer.plot(); temp_viewer.plot()
			Z_viewer.plot(); D_viewer.plot()

			if config.aux_plots == True:
				for current_aux in aux_plot_array:
					current_aux.plot()

		# Save TSV's
		if config.save_TSVs == True:
			TSVViewer(vars=\
					(Gamma_an, Gamma_cx, Gamma_bulk, Gamma_OL,\
					Z_transient_coeff, Z_diffusion_coeff)).plot(filename=\
					config.save_directory+"/"+str(t).zfill(4)+".tsv")


	raw_input(" <=============== End of Program. Press any key to continue. ===============> ")

