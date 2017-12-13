from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm

from fipy.tools import numerix

from parameters import *

nx = 100
L = 5.0
mesh = Grid1D(nx=nx, Lx=L)

x = mesh.cellCenters[0]


# ----------------- Variable Declarations -----------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True)

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True)

Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True)


# ------------- Initial Conditions and Diffusivity---------
Z0 = Z_S*(1 - numerix.tanh((L*x - L) / 2.0))
Z.setValue(Z0)

# Zohm's model
D_Zohm = (D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(Z)) / 2.0
# Stap's Model
alpha_sup = 0.5
D_Staps = D_min + (D_max - D_min) / (1 + alpha_sup*numerix.dot(Z.grad, Z.grad))
# Flow-Shear Model
a1, a3 = 1.0, 0.5	# ASSUMES a2 = 0
D_Shear = D_min + (D_max - D_min) / (1 + a1*Z**2 + a3*numerix.dot(Z.grad, Z.grad))

# CHOOSE DIFFUSIVITY HERE!
Diffusivity = D_Shear


density0 = CellVariable(name=r"$n_0$", mesh=mesh, value=-(Gamma_c*lambda_n / Diffusivity) * (1 + x/lambda_n))
density.setValue(density0)

temp0 = CellVariable(name=r"$T_0", mesh=mesh, value = q_c*((gamma - 1) / Gamma_c) * (1 - lambda_n / (zeta*lambda_T + lambda_n)*(1 + x/lambda_n)**(-zeta)))
temperature.setValue(temp0)

# ----------------- Printing for testing ------------------

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

temp_right = (zeta / Diffusivity.faceValue)*(temperature.faceValue*Gamma_c - (gamma - 1)*q_c) / density.faceValue
temperature.faceGrad.constrain(temp_right, mesh.facesRight)

"""
	Z Boundary Conditions:
	d/dx(Z(0)) == Z / lambda_Z
	mu*D/epsilon * d/dx(Z(L)) == 0
	d^2/dx^2(Z(0)) == 0
"""
Z.faceGrad.constrain(Z.faceValue / lambda_Z, mesh.facesLeft)
Z.constrain(0.0, mesh.facesRight)

Z.grad.faceGrad.constrain(0.0, mesh.facesLeft)


# ----------------- PDE Declarations ----------------------
# Density Equation
density_equation = TransientTerm(var=density) == DiffusionTerm(coeff=Diffusivity, var=density)

# Temperature Equation
S_T = ((zeta+1)/zeta)*(Diffusivity/density) * numerix.dot(density.grad,temperature.grad)
temp_equation = TransientTerm(var=temperature) == DiffusionTerm(coeff=Diffusivity/zeta, var=temperature) + S_T

# Z Equation
G = a + b*(Z - Z_S) + c*(Z - Z_S)**3
S_Z = (c_T / density)*temperature + (c_n*temperature / density**2)*density + G
Z_equation = TransientTerm(coeff=epsilon, var=Z) == DiffusionTerm(coeff=mu*Diffusivity/epsilon, var=Z) + S_Z


# Fully-Coupled Equation
full_equation = density_equation & temp_equation & Z_equation

initial_viewer = Viewer((density, temperature, Z, Diffusivity))
raw_input("Pause for Initial")

if __name__ == '__main__':
	viewer = Viewer((density, temperature, Z, Diffusivity), datamax=3.0, datamin=-1.5)
	for t in range(100):
		density.updateOld(); temperature.updateOld()
		Z.updateOld()#; D.updateOld()
		full_equation.solve(dt=0.01)
		viewer.plot()

	raw_input("End of Program. <return> to continue...")

