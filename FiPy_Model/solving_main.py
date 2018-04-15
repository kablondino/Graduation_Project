#!/usr/bin/env python
from fipy import TransientTerm, DiffusionTerm, Viewer, TSVViewer
from fipy.solvers import *

# Order of file imports from the following inport: input_handling.py,
# parameters.py, variable_decl.py boundary_init_cond
from calculate_coeffs import *
# fipy.tools.numerix is also imported from the above

import os	# For saving files to a specified directory


# ----------------- NOTE ----------------------------------
"""
	This file is never meant to solve the system, and
	therefore, the main solving routine has been deleted.
	Refer to the master branch to obtain it.

	In addition, the 'original model' is also not referenced
	here whatsoever.
"""

# ----------------- PDE Declarations ----------------------
# Density Equation
density.equation = TransientTerm(coeff=1.0, var=density)\
		== DiffusionTerm(coeff=Diffusivity, var=density)

## Energy Equation
temperature.equation = TransientTerm(coeff=density, var=temperature)\
		== DiffusionTerm(coeff=(Diffusivity*density/zeta), var=temperature)\
		+ DiffusionTerm(coeff=Diffusivity*temperature, var=density)

# Z Equation, Full Flux model
Z_transient_coeff = m_i * density * temperature\
		/ (charge**2 * rho_pi * B**2)
Z_transient_coeff.name = r"$\hat{epsilon}$"
Z_diffusion_coeff = (m_i * mu * density * temperature\
		/ (charge**2 * rho_pi * B_theta**2))
Z_diffusion_coeff.name = r"$\hat{mu}$"

Z.equation = TransientTerm(coeff=Z_transient_coeff, var=Z)\
		== DiffusionTerm(coeff=Z_diffusion_coeff, var=Z)\
		+ Gamma_an - Gamma_cx - Gamma_bulk - Gamma_OL

# Fully-Coupled Equation
full_equation = density.equation & temperature.equation & Z.equation

# LOAD pickled H--Mode data for Z
#if __name__ == '__main__' and config.nx == 500:
#	print "Using pickled H--mode data for Z."
#	H_mode_data = dump.read("./L_start_pickle/state0100.dat")
#	Z.setValue(H_mode_data['Z'])
#	Diffusivity.setValue(D_choice_local)

# Initial conditions viewer
init_density_viewer = Viewer(density, xmin=0.0, xmax=L,\
		legend=None, title="Density")
init_temp_viewer = Viewer(temperature/charge, xmin=0.0,\
		xmax=L, legend=None, title="Temperature")
init_Z_viewer = Viewer((-Z, Diffusivity), xmin=0.0,\
		xmax=L, legend='best',\
		title=r"$Z$ and Diffusivity")
raw_input("Pause for SI Initial Conditions")


all_variables = (density, temperature, Z, Diffusivity, v_Ti, v_Te, rho_pi,\
		rho_pe, omega_t, omega_bi, omega_be, w_bi, nu_ei, nu_ii, nu_in0,\
		nu_eff, nu_ai, nu_ae, D_an, g_n_an, g_T_an, g_Z_an, Gamma_an,\
		g_n_cx, g_T_cx, g_Z_cx, Gamma_cx, Gamma_bulk, g_OL, Gamma_OL)
#TSVViewer(vars=all_variables).plot(filename="all_vars.tsv")


# Debug
calculate_coeffs()
print_variables(density, temperature, Z, Diffusivity, v_Ti, v_Te, rho_pi,\
		rho_pe, omega_t, omega_bi, omega_be, w_bi, nu_ei, nu_ii, nu_in0,\
		nu_eff, nu_ai, nu_ae, D_an, g_n_an, g_T_an, g_Z_an, Gamma_an,\
		g_n_cx, g_T_cx, g_Z_cx, Gamma_cx, D_bulk, Gamma_bulk, g_OL, Gamma_OL,\
		Z_transient_coeff, Z_diffusion_coeff, Flux_coeff)


