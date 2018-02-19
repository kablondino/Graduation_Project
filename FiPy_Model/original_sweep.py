from fipy import TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm, TSVViewer

from fipy.solvers import *

from variable_decl import *
from coeffs import *

# --------- Set Initial Conditions for n and T ------------
# Initial conditions for L--mode
density0L = CellVariable(name=r"$n_{0L}$", mesh=mesh,\
		value=-(Gamma_c*lambda_n / Diffusivity) * (1.0 + x/lambda_n))

temp0L = CellVariable(name=r"$T_{0L}", mesh=mesh,\
		value = q_c*((gamma - 1.0) / Gamma_c) * \
		(1.0 - lambda_n / (zeta*lambda_T + lambda_n)*(1.0 + x/lambda_n)**(-zeta)))

# Initial conditions for H--mode
density0H = CellVariable(name=r"$n_0$", mesh=mesh,\
		value=-(Gamma_c*lambda_n / Diffusivity) * (1.0 + x/lambda_n))

temp0H = CellVariable(name=r"$T_0", mesh=mesh,\
		value = q_c*((gamma - 1.0) / Gamma_c) *\
		(1.0 - lambda_n / (zeta*lambda_T + lambda_n)*(1.0 + x/lambda_n)**(-zeta)))

density.setValue(density0L)
temperature.setValue(temp0L)
U.setValue(density*temperature / (gamma - 1.0))

# ----------------- Boundary Conditions -------------------
"""
	Density Boundary Conditions:
	d/dx(n(0)) == n / lambda_n
	d/dx(n(L)) == -Gamma_c / Diffusivity

	Temperature Boundary Conditions:
	d/dx(T(0)) = T / lambda_T
	d/dx(T(L)) = zeta*(Gamma_c*T - q_c*(gamma - 1)) / (Diffusivity * n)

	Paquay considered these Z Boundary Conditions:
	d/dx(Z(0)) == Z / lambda_Z
	mu*D/epsilon * d/dx(Z(L)) == 0
"""
def set_boundary_values(AGamma_c, Aq_c):
	density.faceGrad.constrain(density.faceValue / lambda_n, mesh.facesLeft)

	density.faceGrad.constrain(-AGamma_c / Diffusivity.faceValue, mesh.facesRight)

	temp_left = temperature.faceValue / lambda_T
	temperature.faceGrad.constrain(temp_left, mesh.facesLeft)
	temp_right = (zeta * (AGamma_c*temperature.faceValue - Aq_c*(gamma - 1.0)))\
					/ (Diffusivity.faceValue * density.faceValue)
	temperature.faceGrad.constrain(temp_right, mesh.facesRight)

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
GMRES_Solver = LinearGMRESSolver(iterations=1000, tolerance=1.0e-4)

## Options for viewers
x_min, x_max, y_min, y_max = 0.0, L, 0.0, 5.0

#initial_viewer = Viewer((density, temperature, Z, Diffusivity),\
#		xmin=x_min, xmax=x_max, datamin=y_min, legend='best')
#raw_input("Pause for Initial")

timeStep_denom = 30.0
timeStep = epsilon / timeStep_denom
total_time = 250

res_tol = 1.0e-4

if __name__ == '__main__':
	viewer = Viewer((density, temperature, Z, Diffusivity), xmin=x_min, xmax=x_max,\
			datamin=y_min, legend='best',\
			title="Default $t = $"+str(total_time)+r", $\Delta t = \epsilon / $"+str(timeStep_denom))
	for t in range(total_time):
		# (Re)set residual value(s)
		res_D = res_n = res_T = res_Z = res_full = 1.0e10
		# Update values
		Diffusivity.setValue(D_choice)
		density.updateOld(); temperature.updateOld(); Z.updateOld()

		# Solve each equation individually
#		while max(res_n, res_T) > res_tol:
#			print t, res_n, res_T, res_Z
#			res_n = density.equation.sweep(var=density, dt=timeStep, solver=GMRES_Solver)
#			res_T = temperature.equation.sweep(var=temperature, dt=timeStep, solver=GMRES_Solver)
#			res_Z = Z.equation.sweep(var=Z, dt=timeStep, solver=GMRES_Solver)
#			viewer.plot()

		# Solve the fully coupled equation
		while res_full > res_tol:
			print t, res_full
			res_full = full_equation.sweep(dt=timeStep)#, solver=GMRES_Solver)

#		viewer.plot()
		viewer.plot(filename="original_solution/"+str(t).zfill(4)+".png")

	raw_input(" =================================== End of Program. Press any key to continue. ==================================== ")

