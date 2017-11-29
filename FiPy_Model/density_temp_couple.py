from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm

from fipy.tools import numerix

nx = 100
L = 5.0
mesh = Grid1D(nx=nx, Lx=L)

x = mesh.cellCenters[0]

# Parameters
D = 1.0
zeta = 0.5
Gamma_c = -4.0/5.0
q_c = -4.0
gamma = 5.0/3.0
lambda_n = 5.0/4.0
lambda_T = 3.0/2.0

# Initial Conditions
density0 = CellVariable(name=r"$n_0$", mesh=mesh, value=-(Gamma_c*lambda_n / D) * (1 + x/0.5))

temp0 = CellVariable(name=r"$T_0", mesh=mesh, value = q_c*((gamma - 1) / Gamma_c) * (1 - lambda_n / (zeta*lambda_T + lambda_n)*(1 + x/lambda_n)**(-zeta)))

# Variable Declarations
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True, value=density0)

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True, value=temp0)

## Density Boundary Conditions:
#	d/dx(n(0,t)) == n / lambda_n
#	d/dx(n(L,t)) == D / Gamma_c
density_left = density / lambda_n
density.faceGrad.constrain(density_left, mesh.facesLeft)
density_right = -Gamma_c / D
density.faceGrad.constrain(density_right, mesh.facesRight)

## Temperature Boundary Conditions:
#	d/dx(T(0,t)) = T / lambda_T
#	d/dx(T(L,t)) = (T*Gamma_c - (gamma - 1)*q_c) / (diffusivity * n)
temp_left = temperature / lambda_T
temperature.faceGrad.constrain(temp_left, mesh.facesLeft)
# THE FOLLOWING NEEDS TO BE FIXED, since the index gets messed up
temp_right = (zeta / D)#*(temperature*Gamma_c - (gamma - 1)*q_c) / density
temperature.faceGrad.constrain(temp_right, mesh.facesRight)

# Density Equation
density_equation = TransientTerm(var=density) == DiffusionTerm(coeff=D, var=density)

# Temperature Equation
S_T = ((zeta+1)/zeta)*(D/density)*numerix.dot(density.grad,temperature.grad)
temp_equation = TransientTerm(var=temperature) == DiffusionTerm(coeff=D/zeta, var=temperature) + S_T


# Coupled Equation
full_equation = density_equation & temp_equation
#print full_equation

viewer = Viewer((density, temperature), datamin=0.0)

if __name__ == '__main__':
	for t in range(100):
		density.updateOld()
		temperature.updateOld()
		full_equation.solve(dt=1.0e-3)
		viewer.plot()

	raw_input("End of Program. <return> to continue...")

