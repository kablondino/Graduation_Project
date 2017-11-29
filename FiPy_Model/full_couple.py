from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm

from fipy.tools import numerix

nx = 100
L = 5.0
mesh = Grid1D(nx=nx, Lx=L)

x = mesh.cellCenters[0]

# Parameters
zeta = 0.5
Gamma_c = -4.0/5.0
q_c = -4.0
gamma = 5.0/3.0

lambda_n = 5.0/4.0
lambda_T = 3.0/2.0
lambda_Z = 5.0/4.0

epsilon = 1.0/25.0
mu = 1.0/20.0

c_n = -1.1
c_T = -0.9

a = 3.0/2.0
b = 2.0
c = -1.0
Z_S = -3.0/2.0


# ----------------- Variable Declarations -----------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True)

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True)

Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True)

## Diffusivity
D_max = 2.0
D_min = 2.0/5.0
#D = CellVariable(name=r"$D$", mesh=mesh, hasOld=True)
D = 5.0
#def D(aZ):
#	return (D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(aZ)) / 2.0


# ----------------- Initial Conditions --------------------
density0 = CellVariable(name=r"$n_0$", mesh=mesh, value=-(Gamma_c*lambda_n / D) * (1 + x/0.5))
density.setValue(density0)

temp0 = CellVariable(name=r"$T_0", mesh=mesh, value = q_c*((gamma - 1) / Gamma_c) * (1 - lambda_n / (zeta*lambda_T + lambda_n)*(1 + x/lambda_n)**(-zeta)))
temperature.setValue(temp0)

Z0 = Z_S*(1 - numerix.tanh((L*x - L) / 2))
Z.setValue(Z0)

#D0 = (D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(Z)) / 2.0
#D.setValue(D0)
print D

# ----------------- Boundary Conditions -------------------
# Density Boundary Conditions:
#	d/dx(n(0,t)) == n / lambda_n
#	d/dx(n(L,t)) == D / Gamma_c
density_left = density / lambda_n
density.faceGrad.constrain(density_left, mesh.facesLeft)

density_right = -Gamma_c / D
density.faceGrad.constrain(density_right, mesh.facesRight)

# Temperature Boundary Conditions:
#	d/dx(T(0,t)) = T / lambda_T
#	d/dx(T(L,t)) = (T*Gamma_c - (gamma - 1)*q_c) / (diffusivity * n)
temp_left = temperature / lambda_T
temperature.faceGrad.constrain(temp_left, mesh.facesLeft)
# THE FOLLOWING NEEDS TO BE FIXED, since the index gets messed up
temp_right = (zeta / D)#*(temperature*Gamma_c - (gamma - 1)*q_c) / density
temperature.faceGrad.constrain(temp_right, mesh.facesRight)

# Z Boundary Conditions:
#	d/dx(Z(0,t)) = Z / lambda_Z
#	mu*D/epsilon * d/dx(Z(L,t)) == 0  ---->  Z == 0
Z_left = Z / lambda_Z
Z.faceGrad.constrain(Z_left, mesh.facesLeft)

Z.constrain(0.0, mesh.facesRight)


# ----------------- PDE Declarations ----------------------
# Density Equation
density_equation = TransientTerm(var=density) == DiffusionTerm(coeff=5.0, var=density)

# Temperature Equation
S_T = ((zeta+1)/zeta)*(D/density)*numerix.dot(density.grad,temperature.grad)
temp_equation = TransientTerm(var=temperature) == DiffusionTerm(coeff=5.0/0.5, var=temperature) + S_T

# Z Equation
G = a + b*(Z - Z_S) + c*(Z - Z_S)**3
#S_Z = numerix.dot((c_n*temperature / density**2),density.grad) + numerix.dot((c_T / density), temperature.grad) + G
#S_Z = ConvectionTerm(coeff=((c_n*temperature/density**2),), var=density) + ConvectionTerm(coeff=((c_T / density),), var=temperature) + G
S_Z = G# + -density * temperature.grad
Z_equation = TransientTerm(coeff=epsilon, var=Z) == DiffusionTerm(coeff=mu*D, var=Z) + S_Z


# Fully-Coupled Equation
full_equation = density_equation & temp_equation & Z_equation
print full_equation

if __name__ == '__main__':
	viewer = Viewer((density, temperature, Z))
	initial_view = Viewer((Z))
	raw_input("Pause")

#if __name__ == '__main__':
#	for t in range(100):
#		density.updateOld()
#		temperature.updateOld()
#		full_equation.solve(dt=1.0e-3)
#		viewer.plot()
#
#	raw_input("End of Program. <return> to continue...")

