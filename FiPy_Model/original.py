from fipy import TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm, TSVViewer

from variable_decl import *
from coeffs import *

# ----------------- Printing for testing ------------------
def printing():
	print rho_pi
	print nu_ai
	print g_n_an

# --------- Set Initial Conditions for n and T ------------
density.setValue(density0)
temperature.setValue(temp0)

# ----------------- Boundary Conditions -------------------
"""
	Density Boundary Conditions:
	d/dx(n(0)) == n / lambda_n
	d/dx(n(L)) == D / Gamma_c
"""
density.faceGrad.constrain(density.faceValue / lambda_n, mesh.facesLeft)

density.faceGrad.constrain(-Gamma_c / Diffusivity.faceValue, mesh.facesRight)

"""
	Temperature Boundary Conditions:
	d/dx(T(0)) = T / lambda_T
	d/dx(T(L)) = (T*Gamma_c - (gamma - 1)*q_c) / (diffusivity * n)
"""
temp_left = temperature.faceValue / lambda_T
temperature.faceGrad.constrain(temp_left, mesh.facesLeft)

temp_right = (zeta / Diffusivity.faceValue)*(temperature.faceValue*Gamma_c\
				- (gamma - 1.0)*q_c) / density.faceValue
temperature.faceGrad.constrain(temp_right, mesh.facesRight)

"""
	Z Boundary Conditions:
	d/dx(Z(0)) == Z / lambda_Z
	mu*D/epsilon * d/dx(Z(L)) == 0
	d^2/dx^2(Z(0)) == 0
"""
Z.faceGrad.constrain(Z.faceValue / lambda_Z, mesh.facesLeft)
Z.constrain(0.0, mesh.facesRight)

Z.faceGrad.divergence.constrain(0.0, mesh.facesLeft)

"""
	Zohm's N(Z,g) = 0 Boundary:
	g_0 * T/n * ((dn/dx) / n + gamma*(dT/dx) / T)
		= g_1 + alpha*Z - beta*Z^3
"""
#N_boundary = ((density.faceValue)**2 / (temperature.faceValue))\
#		* (a + b*(Z.faceValue) - c*(Z.faceValue)**3) - gamma*density.faceValue\
#		* temperature.faceGrad / temperature.faceValue
#density.faceGrad.constrain(N_boundary, mesh.facesRight)
#density.faceGrad.constrain(N_boundary, mesh.facesLeft)

# ----------------- PDE Declarations ----------------------
# Density Equation
density_equation = TransientTerm(coeff=1.0, var=density)\
		== DiffusionTerm(coeff=Diffusivity, var=density)

# Temperature Equation
temp_equation = TransientTerm(coeff=density, var=temperature)\
		== DiffusionTerm(coeff=(Diffusivity*density/zeta), var=temperature)\
		+ DiffusionTerm(coeff=Diffusivity*temperature, var=density)
#S_T = ((zeta + 1.0)/zeta) * (Diffusivity / density) * numerix.dot(density.grad, temperature.grad) + (Diffusivity*temperature)/density * density.faceGrad.divergence
#temp_equation = TransientTerm(coeff=1.0, var=temperature) == DiffusionTerm(coeff=(Diffusivity/zeta), var=temperature) + S_T

# Z Equation
G = a + b*(Z - Z_S) + c*(Z - Z_S)**3
S_Z = ((c_n*temperature) / density**2) * density.grad[0]\
		+ (c_T / density) * temperature.grad[0] + G
Z_equation = TransientTerm(coeff=epsilon, var=Z) == DiffusionTerm(coeff=mu, var=Z) + S_Z

# Fully-Coupled Equation
full_equation = density_equation & temp_equation & Z_equation

#printing()

## Options for viewers
x_min, x_max, y_min, y_max = 0.0, L, -1.0, 5.0

#initial_viewer = Viewer((density, temperature, Z, Diffusivity), xmin=x_min, xmax=x_max, datamin=y_min, datamax=y_max, legend='best')
#raw_input("Pause for Initial")

timeStep = epsilon / 3.0

if __name__ == '__main__':
	viewer = Viewer((density, temperature, Z), xmin=x_min,\
			xmax=x_max, legend='best')
	for t in range(100):
		density.updateOld(); temperature.updateOld(); Z.updateOld()
		#Diffusivity.updateOld()	# Remanent of when it was a CellVariable
#		TSVViewer(vars=(density, temperature, Z, Diffusivity, rho_pi, n_0, nu_ei, nu_ii, nu_ai)).plot(filename="original_solution/original"+str(t)+".tsv")
		full_equation.solve(dt=timeStep)
		viewer.plot()

	raw_input("End of Program. <return> to continue...")

