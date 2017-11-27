## This file calculates the density with naive diffusivity

from fipy import Variable, CellVariable, Grid1D, TransientTerm, DiffusionTerm, Viewer

from fipy.tools import numerix

L = 5.0
nx = 50
dx = L / nx
mesh = Grid1D(nx=nx, dx=dx)

# Parameters from core, SOL, etc.
Gamma_c = -4.0/5.0
lambda_n = 5.0/4.0

# Diffusivity function(s) OR constant
diffusivity = 5

timeStepDuration = 9 * dx**2 / (2*diffusivity)
steps = 50

x = mesh.cellCenters[0]
t = timeStepDuration * steps

# Solution variable
density = CellVariable(name=r"$n$", mesh=mesh)

# Initial Condition
density_initial = CellVariable(name=r"$n_0$", mesh=mesh, value=-Gamma_c*lambda_n / diffusivity * (1 + x/lambda_n))
density.setValue(density_initial)

# Boundary conditions:
#	D*n(0,t) / lambda_n + D*d/dx(n(0,t)) == 0
#	Gamma_c(t) + D*d/dx(n(L,t)) == 0
density_left_neumann = density / lambda_n
density.faceGrad.constrain(density_left_neumann, mesh.facesLeft)

density_right_neumann = -Gamma_c / diffusivity
density.constrain(density_right_neumann, mesh.facesRight)


# Set the equation
density_equation = TransientTerm(var=density) == DiffusionTerm(coeff=diffusivity, var=density)
 
if __name__ == '__main__':
	viewer = Viewer(vars=(density, density_initial), datamin=0.0, datamax=max(density_initial))

for step in range(100):
	density_equation.solve(var=density, dt=timeStepDuration)
	if __name__ == '__main__':
		viewer.plot()

if __name__ == '__main__':
	raw_input("Pause. Press <return> to continue.")

