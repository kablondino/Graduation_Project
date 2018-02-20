from fipy import TransientTerm, DiffusionTerm, Viewer, TSVViewer

from fipy.solvers import *

# variable_decl.py, parameters.py, and fipy.tools.numerix are included in the following
from coeffs import *

# ----------------- Boundary Conditions -------------------
def set_boundary_values(AGamma_c, Aq_c):
	"""
		Density Boundary Conditions:
		d/dx(n(0)) == n / lambda_n
		d/dx(n(L)) == -Gamma_c / Diffusivity
	"""
	density.faceGrad.constrain(density.faceValue / lambda_n, mesh.facesLeft)
	density.faceGrad.constrain(-AGamma_c / Diffusivity.faceValue, mesh.facesRight)

	"""
		Temperature Boundary Conditions:
		d/dx(T(0)) = T / lambda_T
		d/dx(T(L)) = zeta*(Gamma_c*T - q_c*(gamma - 1)) / (Diffusivity * n)
	"""
	temp_left = temperature.faceValue / lambda_T
	temperature.faceGrad.constrain(temp_left, mesh.facesLeft)
	temp_right = (zeta * (AGamma_c*temperature.faceValue - Aq_c*(gamma - 1.0)))\
					/ (Diffusivity.faceValue * density.faceValue)
	temperature.faceGrad.constrain(temp_right, mesh.facesRight)

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
x_min, x_max, y_min, y_max = 0.0, L, 0.0, 5.0

#initial_viewer = Viewer((density, temperature, Z, Diffusivity),\
#		xmin=x_min, xmax=x_max, datamin=y_min, legend='best')
#raw_input("Pause for Initial")

timeStep_denom = 50.0
timeStep = epsilon / timeStep_denom
total_time = 300

res_tol = 1.0e-5


if __name__ == '__main__':
	viewer = Viewer((density, temperature, -Z, Diffusivity), xmin=x_min, xmax=x_max,\
			legend='best',\
			title="GMRES H--Mode Start; $t = $"+str(total_time)+r", $\Delta t = \epsilon / $"+str(timeStep_denom))
	for t in range(total_time):
		# (Re)set residual value(s)
		res_D = res_n = res_T = res_Z = res_full = 1.0e10
		# Update values
		Diffusivity.setValue(D_choice)
		density.updateOld(); temperature.updateOld(); Z.updateOld()

		# Solve the fully coupled equation
		while res_full > res_tol:
			print t, res_full
			res_full = full_equation.sweep(dt=timeStep, solver=GMRES_Solver)

		viewer.plot()
#		if t % 10 == 0:
#		viewer.plot(filename="original_solution/"+str(t).zfill(4)+".png")

	raw_input(" =================================== End of Program. Press any key to continue. ==================================== ")

