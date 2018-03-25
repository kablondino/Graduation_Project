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
	Z.equation = TransientTerm(coeff=(m_i * density * temperature\
			/ (charge_true**2 * rho_pi * B**2)), var=Z) == DiffusionTerm(\
			coeff=(m_i * mu / (charge_true**2 * rho_pi * B_theta**2)), var=Z)\
			+ Gamma_an - Gamma_cx - Gamma_bulk - Gamma_OL

elif config.original_model == True:
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
if config.original_model == True:
	initial_viewer = Viewer((density, temperature, -Z, Diffusivity),\
			xmin=0.0, xmax=config.plotx_max, legend='best')
	raw_input("Pause for Initial Conditions")
elif config.original_model == False:
	init_density_viewer = Viewer(density, xmin=0.0, xmax=config.plotx_max,\
			legend=None, title="Density")
	init_density_viewer = Viewer(temperature, xmin=0.0,\
			xmax=config.plotx_max, legend=None, title="Temperature")
	init_density_viewer = Viewer((Z, Diffusivity), xmin=0.0,\
			xmax=config.plotx_max, legend='best',\
			title=r"$Z$ and Diffusivity")
	raw_input("Pause for SI Initial Conditions")

timeStep = epsilon / config.timeStep_denom

if __name__ == '__main__':

	# Initialize viewers
	if config.original_model == True:
		viewer = Viewer((density, temperature, -Z, Diffusivity),\
				xmin=0.0, xmax=config.plotx_max,\
				datamax=config.ploty_max, legend='best',\
				title = config.plot_title)
	else:
		density_viewer = Viewer(density, xmin=0.0, xmax=config.plotx_max,\
				datamax=3.0e20, legend='best',\
				title = config.plot_title)
		temp_viewer = Viewer(temperature, xmin=0.0, xmax=config.plotx_max,\
				datamax=2e3, legend='best',\
				title = config.plot_title)
		Z_viewer = Viewer((Z, Diffusivity), xmin=0.0, xmax=config.plotx_max,\
				legend='best',\
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

	# Set the tolerance for the full flux model
	original_res_tol = 1.0e-7
	density_res_tol, temp_res_tol, Z_res_tol = 1.0e12, 1.0e12, 1.0e12

	for t in range(config.total_timeSteps):
		# (Re)set residual values
		original_residual = 1.0e100
		density_residual = 1.0e100
		temp_residual = 1.0e100
		Z_residual = 1.0e100

		# Update values
		Diffusivity.setValue(D_choice_local)
		density.updateOld(); temperature.updateOld(); Z.updateOld()
		update_g_coeffs()

		# Solve the full flux model equations
		if config.original_model == False:
			# Solve density equation
			while density_residual > density_res_tol:
				print t, density_residual, temp_residual, Z_residual
				density_residual = density.equation.sweep(dt=timeStep,\
							solver=GMRES_Solver)

			# Solve temperature equation
			while temp_residual > temp_res_tol:
				print t, density_residual, temp_residual, Z_residual
				temp_residual = temperature.equation.sweep(dt=timeStep,\
							solver=GMRES_Solver)

			# Solve Z equation
			while Z_residual > Z_res_tol:
				print t, density_residual, temp_residual, Z_residual
				Z_residual = Z.equation.sweep(dt=timeStep,\
						solver=GMRES_Solver)

		# Solve the fully coupled, original model equation
		else:
			config.original_model == True
			while original_residual > original_res_tol:
				print t, original_residual
				original_residual = full_equation.sweep(dt=timeStep,\
						solver=GMRES_Solver)

		# Plot solution and save, if option is True
		if config.save_plots == True:
			# Save full flux plots
			if config.original_model == False:
				density_viewer.plot(filename = config.save_directory + "/n"\
						+str(t).zfill(4)+ ".png")
				temp_viewer.plot(filename = config.save_directory + "/T" \
						+str(t).zfill(4)+ ".png")
				Z_viewer.plot(filename = config.save_directory + "/Z" \
						+str(t).zfill(4)+ ".png")
			# Save original model plots
			else:
				viewer.plot(filename =\
						config.save_directory+"/"+str(t).zfill(4)+".png")

			# Save auxiliary plots
			if config.aux_plots == True:
				auxiliary1_viewer.plot(filename =\
						config.save_directory+"/aux1_"+str(t).zfill(4)+".png")
				auxiliary2_viewer.plot(filename =\
						config.save_directory+"/aux2_"+str(t).zfill(4)+".png")

		else:
			if config.original_model == False:
				density_viewer.plot(); temp_viewer.plot(); Z_viewer.plot()
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

