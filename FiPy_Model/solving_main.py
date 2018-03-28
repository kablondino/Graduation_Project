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

## Energy Equation
temperature.equation = TransientTerm(coeff=density, var=temperature)\
		== DiffusionTerm(coeff=(Diffusivity*density/zeta), var=temperature)\
		+ DiffusionTerm(coeff=Diffusivity*temperature, var=density)

# Z Equation
if config.original_model == False:
	# Full Flux model
	Z_transient_coeff = m_i * density * temperature\
			/ (charge**2 * rho_pi * B**2)
	Z_diffusion_coeff = (m_i * mu / (charge**2 * rho_pi * B_theta**2))

	Z.equation = TransientTerm(coeff=Z_transient_coeff, var=Z)\
			== DiffusionTerm(coeff=Z_diffusion_coeff, var=Z)\
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
init_density_viewer = Viewer(density, xmin=0.0, xmax=L,\
		legend=None, title="Density")
if config.original_model == True:
	init_temp_viewer = Viewer(temperature, xmin=0.0,\
			xmax=L, legend=None, title="Temperature")
elif config.original_model == False:
	init_temp_viewer = Viewer(temperature/charge, xmin=0.0,\
			xmax=L, legend=None, title="Temperature")
init_Z_viewer = Viewer((-Z, Diffusivity), xmin=0.0,\
		xmax=L, legend='best',\
		title=r"$Z$ and Diffusivity")
raw_input("Pause for SI Initial Conditions")

timeStep = epsilon / config.timeStep_denom

# Debug
#update_g_coeffs()
print "\n ------------------- Density ---------------------"
print density
print "\n ------------------- Temperature -----------------"
print temperature
#print "\n ------------------- Z, in H--mode ---------------"
#print Z
#print "\nManual E_r calc"
#print (Z * temperature / (charge*rho_pi)).inBaseUnits()
#print "\n ------------------- Diffusivity - ---------------"
#print Diffusivity
#print "\n ------------------- v_Ti ------------------------"
#print v_Ti
#print "\n ------------------- v_Te ------------------------"
#print v_Te
#print "\n ------------------- rho_pi ----------------------"
#print rho_pi
#print "\n ------------------- rho_pe ----------------------"
#print rho_pe
## Frequencies
#print "\n ------------------- omega_t ---------------------"
#print omega_t
#print "\n ------------------- omega_bi --------------------"
#print omega_bi
#print "\n ------------------- omega_be --------------------"
#print omega_be
#print "\n ------------------- w_bi ------------------------"
#print w_bi
#print "\n ------------------- nu_ei -----------------------"
#print nu_ei
#print "\n ------------------- nu_ii -----------------------"
#print nu_ii
#print "\n ------------------- nu_in0 ----------------------"
#print nu_in0
#print "\n ------------------- nu_eff ----------------------"
#print nu_eff
#print "\n ------------------- nu_ai -----------------------"
#print nu_ai
#print "\n ------------------- nu_ae -----------------------"
#print nu_ae
#
## Fluxes
#print "\n ------------------- D_an ------------------------"
#print D_an
#print "\n ------------------- g_n_an ----------------------"
#print g_n_an
#print "\n ------------------- g_T_an ----------------------"
#print g_T_an
#print "\n ------------------- g_Z_an ----------------------"
#print g_Z_an
#print "\n ------------------- Gamma_an --------------------"
#print Gamma_an
#print "\n ------------------- g_n_cx ----------------------"
#print g_n_cx
#print "\n ------------------- g_T_cx ----------------------"
#print g_T_cx
#print "\n ------------------- g_Z_cx ----------------------"
#print g_Z_cx
#print "\n ------------------- Gamma_cx --------------------"
#print Gamma_cx
#print "\n ------------------- Gamma_bulk ------------------"
#print Gamma_bulk
#print "\n ------------------- g_OL ------------------------"
#print g_OL
#print "\n ------------------- Gamma_OL --------------------"
#print Gamma_OL
#
## Transient coefficient
#print "\n ------------------- Transient Coeff -------------"
#print (Z_transient_unit*m_i * density * temperature / (charge**2 * rho_pi * B**2)).inBaseUnits()
## Diffusion coefficient
#print "\n ------------------- Diffusion Coeff -------------"
#print (Z_diffusion_unit*m_i * mu *density*temperature / (charge**2 * rho_pi * B_theta**2)).inBaseUnits()

all_variables = (density, temperature, Z, Diffusivity, v_Ti, v_Te, rho_pi,\
		rho_pe, omega_t, omega_bi, omega_be, w_bi, nu_ei, nu_ii, nu_in0,\
		nu_eff, nu_ai, nu_ae, D_an, g_n_an, g_T_an, g_Z_an, Gamma_an,\
		g_n_cx, g_T_cx, g_Z_cx, Gamma_cx, g_n_bulk, g_T_bulk, g_Z_bulk,\
		Gamma_bulk, g_OL, Gamma_OL)
#TSVViewer(vars=all_variables).plot(filename="all_vars.tsv")

# ----------------- NOTE ----------------------------------
"""
	This file is never meant to solve the system, and
	therefore, the main solving routine has been deleted.
	Refer to the master branch to obtain it.
"""

