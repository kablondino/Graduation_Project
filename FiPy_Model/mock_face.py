"""
	The entire point of this file is to experiment with
	the relationship between face and cell variables.
"""

from fipy import TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm

from variable_decl import *
#from coeffs import *


# ----------------- Boundary Conditions -------------------
"""
	Density Boundary Conditions:
	d/dx(n(0)) == n / lambda_n
	d/dx(n(L)) == D / Gamma_c
"""
# I think density.faceValue defaults to the leftmost value
# for mesh.facesLeft, and rightmost for mesh.facesRight
density.faceGrad.constrain(density.faceValue / lambda_n, mesh.facesLeft)

density.faceGrad.constrain(-Gamma_c / Diffusivity_f, mesh.facesRight)

"""
	Temperature Boundary Conditions:
	d/dx(T(0)) = T / lambda_T
	d/dx(T(L)) = (T*Gamma_c - (gamma - 1)*q_c) / (diffusivity * n)
"""
temp_left = temperature.faceValue / lambda_T
temperature.faceGrad.constrain(temp_left, mesh.facesLeft)

temp_right = (zeta / Diffusivity_f)*(temperature.faceValue*Gamma_c - (gamma - 1.0)*q_c) / density.faceValue
temperature.faceGrad.constrain(temp_right, mesh.facesRight)

# ----------------- PDE Declarations ----------------------
# Density Equation
density_equation = TransientTerm(coeff=1.0, var=density) == DiffusionTerm(coeff=Diffusivity_c, var=density)

# Temperature Equation
temp_equation = TransientTerm(coeff=density, var=temperature) == DiffusionTerm(coeff=(Diffusivity_c*density/zeta), var=temperature) + DiffusionTerm(coeff=Diffusivity_c*temperature, var=density)

full_equation = density_equation & temp_equation


## Options for viewers
x_min, x_max, y_min, y_max = 0.0, L, -1.0, 5.0

initial_viewer = Viewer((density, temperature, Diffusivity_c), xmin=x_min, xmax=x_max, datamin=y_min, datamax=y_max, legend='best')
raw_input("Pause for Initial")

timeStep = epsilon / 3.0

if __name__ == '__main__':
	viewer = Viewer((density, temperature, Diffusivity_c), xmin=x_min, xmax=x_max, legend='best')
	for t in range(100):
		density.updateOld(); temperature.updateOld()
		Diffusivity_c.updateOld()
		full_equation.solve(dt=timeStep)
		viewer.plot()

	raw_input("End of Program. <return> to continue...")

