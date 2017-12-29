from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm

from fipy.tools import numerix

# Import parameter file
from parameters import *

nx = 100
L = 5.0
mesh = Grid1D(nx=nx, Lx=L)

x = mesh.cellCenters[0]

# ----------------- Variable Declaration ------------------
# u[0] = density
# u[1] = temperature
# u[2] = Z
u = CellVariable(mesh=mesh, hasOld=True, elementshape=(3,))

# Initial Condition for Z
Z0 = Z_S*(1 - numerix.tanh((L*x - L) / 2.0))

# ----------------- Diffusivity models --------------------
# Zohm's model
D_Zohm = (D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(u[2])) / 2.0
# Stap's Model
alpha_sup = 0.5
D_Staps = D_min + (D_max - D_min) / (1 + alpha_sup*numerix.dot(u[2].grad, u[2].grad))
## Flow-Shear Model
a1, a2, a3 = 0.7, 1.25, 0.5
D_Shear = D_min + (D_max - D_min) / (1 + a1*(u[2])**2 + a2*u[2]*(u[2].grad) + a3*numerix.dot(u[2].grad, u[2].grad))

# CHOOSE DIFFUSIVITY HERE!
D_choice = D_Staps

# Initial conditions for density and temperature
density0 = -(Gamma_c*lambda_n / D_choice) * (1 + x/lambda_n)
temp0 = q_c*((gamma - 1.0) / Gamma_c) * (1.0 - lambda_n / (zeta*lambda_T + lambda_n)*(1.0 + x/lambda_n)**(-zeta))

u.setValue((density0, temp0, Z0))

#initial_viewer = Viewer((u[0], u[1], u[2], D_choice))
#raw_input("Pause for Initial")

# ----------------- Boundary Conditions -------------------
"""
	Density Boundary Conditions:
	d/dx(n(0)) == n / lambda_n
	d/dx(n(L)) == D / Gamma_c

	Temperature Boundary Conditions:
	d/dx(T(0)) = T / lambda_T
	d/dx(T(L)) = (T*Gamma_c - (gamma - 1)*q_c) / (diffusivity * n)

	Z Boundary Conditions:
	d/dx(Z(0)) == Z / lambda_Z
	mu*D/epsilon * d/dx(Z(L)) == 0
	d^2/dx^2(Z(0)) == 0
"""
densityL = u[0].faceValue / lambda_n
densityR = -Gamma_c / D_choice.faceValue

tempL = u[1].faceValue / lambda_T
tempR = (zeta / D_choice.faceValue)*(u[1].faceValue*Gamma_c - (gamma - 1)*q_c) / u[0].faceValue

ZL = u[2].faceValue / lambda_Z
ZR = 0.0

u.faceGrad.constrain([densityL, tempL, ZL], mesh.facesLeft)
u.faceGrad.constrain([densityR, tempR, ZR], mesh.facesRight)

#u[2].grad.faceGrad.constrain(0.0, mesh.facesLeft)

# ----------------- Equation ------------------------------
eqn = TransientTerm([1.0, 1.0, epsilon], var=u) == DiffusionTerm([D_choice, D_choice/zeta, mu*D_choice], var=u)

#print eqn

if __name__ == '__main__':
	viewer = Viewer((u[0], u[1], u[2], D_choice), datamax=3.0, datamin=-1.5)
	for t in range(100):
		u.updateOld()
		eqn.solve(dt=0.01)
		viewer.plot()

	raw_input("End of Program. <return> to continue...")

