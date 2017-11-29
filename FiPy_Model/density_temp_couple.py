## This file calculates the density AND temperature with naive diffusivity

from fipy import Variable, CellVariable, Grid1D, TransientTerm, DiffusionTerm, ExplicitDiffusionTerm, ConvectionTerm, ImplicitSourceTerm, Viewer

from fipy.tools import numerix

L = 5.0
nx = 100
dx = L / nx
mesh = Grid1D(nx=nx, dx=dx)

# Parameters from core, SOL, etc.
gamma = 5.0/3.0
Gamma_c = -4.0/5.0
q_c = -4.0
zeta = 0.5
lambda_n = 5.0/4.0
lambda_T = 3.0/2.0

# Diffusivity function(s) OR constant
diffusivity = 5

timeStepDuration = 9 * dx**2 / (2*diffusivity)
steps = 50

x = mesh.cellCenters[0]

# Solution variable
density = CellVariable(name=r"$n$", mesh=mesh)
temperature = CellVariable(name=r"$T$", mesh=mesh)

# Initial Conditions
density_initial = CellVariable(name="Initial Density", mesh=mesh, value=-Gamma_c*lambda_n / diffusivity * (1 + x/lambda_n))
density.setValue(density_initial)

temp_initial = CellVariable(name="Initial Temperature", mesh=mesh, value= q_c * (gamma - 1) / Gamma_c * (1 - (lambda_n / (zeta*lambda_T + lambda_n)) * (1 + x/lambda_n)**(-zeta)))
temperature.setValue(temp_initial)


## Density Boundary Conditions:
#	d/dx(n(0,t)) == n(0,t) / lambda_n
#	d/dx(n(L,t)) == -Gamma_c(t) / D
density_left_neumann = density / lambda_n
density.faceGrad.constrain(density_left_neumann, mesh.facesLeft)
density_right_neumann = -Gamma_c / diffusivity
density.constrain(density_right_neumann, mesh.facesRight)

## Temperature Boundary Conditions:
#	d/dx(T(0,t)) = T / (lambda_T*(1 + (zeta+1)*(1/n)*dn/dx))
#	d/dx(T(L,t)) = (T*Gamma_c - (gamma - 1)*q_c) / (diffusivity * n)
temp_left_neumann = temperature / (lambda_T*(1 + (zeta+1)*(1/density)*density.grad))		# FIX!
temperature.faceGrad.constrain(temp_left_neumann, mesh.facesLeft)
#print temp_left_neumann		# How to represent derivative, not grad!

temp_right_neumann = (zeta*(temperature*Gamma_c - (gamma-1)*q_c)) / (diffusivity*density)	# FIX!
temperature.faceGrad.constrain(temp_right_neumann, mesh.facesLeft)


# Set the equations
density_equation = TransientTerm(var=density) == DiffusionTerm(coeff=diffusivity, var=density)

S_T = (zeta + 1) / zeta * diffusivity / density * (density.grad.mag * temperature.grad.mag)
#print S_T

#print S_T * ConvectionTerm((zeta,0))# * ConvectionTerm((1.0,), var=temperature)
temp_equation = TransientTerm(var=temperature) == ExplicitDiffusionTerm(coeff=diffusivity / zeta, var=temperature) + S_T

semifull_eq = density_equation & temp_equation
print semifull_eq
 
#if __name__ == '__main__':
#	viewer = Viewer((density, temperature), datamin=0.0, datamax=max(temp_initial))
#
#for step in range(steps):
#	semifull_eq.solve(dt=timeStepDuration)
#	if __name__ == '__main__':
#		viewer.plot()

if __name__ == '__main__':
	raw_input("Pause. Press <return> to continue.")
