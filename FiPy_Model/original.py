from fipy import TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm, TSVViewer

from fipy.tools import numerix

from parameters import *
from variable_decl import *

# ------------- Initial Conditions and Diffusivity---------
Z0 = Z_S*(1 - numerix.tanh((L*x - L) / 2.0))
Z.setValue(Z0)

# Zohm's model
D_Zohm = (D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(Z)) / 2.0
# Stap's Model
alpha_sup = 0.5
D_Staps = D_min + (D_max - D_min) / (1.0 + alpha_sup*numerix.dot(Z.grad, Z.grad))
# Flow-Shear Model
a1, a3 = 1.0, 0.5	# ASSUMES a2 = 0
D_Shear = D_min + (D_max - D_min) / (1.0 + a1*Z**2 + a3*numerix.dot(Z.grad, Z.grad))

# CHOOSE DIFFUSIVITY HERE!
D_choice = D_Staps
Diffusivity = D_choice


density0 = CellVariable(name=r"$n_0$", mesh=mesh, value=-(Gamma_c*lambda_n / Diffusivity) * (1.0 + x/lambda_n))
density.setValue(density0)

temp0 = CellVariable(name=r"$T_0", mesh=mesh, value = q_c*((gamma - 1.0) / Gamma_c) * (1.0 - lambda_n / (zeta*lambda_T + lambda_n)*(1.0 + x/lambda_n)**(-zeta)))
temperature.setValue(temp0)

# ----------------- Printing for testing ------------------
def printing():
	print density
	print temperature
	print density*temperature
	
#	printing_viewer = Viewer((density, density.faceGrad.divergence), legend='best')
#	raw_input("Pause for Test")

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

temp_right = (zeta / Diffusivity.faceValue)*(temperature.faceValue*Gamma_c - (gamma - 1.0)*q_c) / density.faceValue
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

# ----------------- PDE Declarations ----------------------
# Density Equation
density_equation = TransientTerm(var=density) == DiffusionTerm(coeff=Diffusivity, var=density)

# Temperature Equation
chi = Diffusivity / (zeta*(gamma - 1.0))
#temp_equation = TransientTerm(coeff=(density/(gamma - 1.0)), var=temperature) == DiffusionTerm(coeff=chi*density, var=temperature) + DiffusionTerm(coeff=Diffusivity*temperature, var=density)
S_T = ((zeta+1.0)/zeta)*(Diffusivity/density) * numerix.dot(density.grad, temperature.grad) + (Diffusivity/density)*numerix.dot(temperature, density.faceGrad.divergence)
temp_equation = TransientTerm(var=temperature) == DiffusionTerm(coeff=Diffusivity/zeta, var=temperature) + S_T

# Z Equation
G = a + b*(Z - Z_S) + c*(Z - Z_S)**3
S_Z = (c_T / density) * temperature.grad[0] + (c_n*temperature / density**2) * density.grad[0] + G
Z_equation = TransientTerm(coeff=epsilon, var=Z) == DiffusionTerm(coeff=mu*Diffusivity/epsilon, var=Z) + S_Z


# Fully-Coupled Equation
full_equation = density_equation & temp_equation & Z_equation

printing()

## Options for viewers
x_min, x_max, y_min, y_max = 0.0, L, -1.0, 5.0

#initial_viewer = Viewer((density, temperature, Z, Diffusivity), xmin=x_min, xmax=x_max, datamin=y_min, datamax=y_max, legend='best')
#raw_input("Pause for Initial")

timeStep = 0.9*epsilon / 5
print timeStep

#if __name__ == '__main__':
#	viewer = Viewer((density, temperature, -Z, Diffusivity), datamin=y_min, datamax=y_max, xmin=x_min, xmax=x_max, legend='best')
#	for t in range(100):
#		density.updateOld(); temperature.updateOld(); Z.updateOld()
#		#Diffusivity.updateOld()	# Remanent of when it was a CellVariable
##		TSVViewer(vars=(density, temperature, Z, Diffusivity)).plot(filename="full_solution/full"+str(t)+".tsv")
#		full_equation.solve(dt=timeStep)
#		viewer.plot()
#
#	raw_input("End of Program. <return> to continue...")

