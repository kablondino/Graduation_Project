from fipy import CellVariable, Grid1D, DiffusionTerm, Viewer, ImplicitSourceTerm
from fipy.tools import numerix

from fipy.solvers import *

from parameters import *

# ----------------- Solver --------------------------------
mySolver = LinearGMRESSolver(iterations=1000, tolerance=1.0e-6)

# ----------------- Mesh Generation -----------------------
nx = 100
L = 5.0
mesh = Grid1D(nx=nx, Lx=L)

x = mesh.cellCenters[0] # Cell position

# ----------------- Variable Declarations -----------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True)
temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True)
Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True)

Diffusivity = CellVariable(name=r"$D$", mesh=mesh, hasOld=True)

# ----------- Initial Conditions of Z ---------------------
Z0L = 1.0 # L--mode
Z0H = Z_S*(1.0 - numerix.tanh((L*x - L) / 2.0)) # H--mode
Z.setValue(Z0H)

# ----------------- Diffusivities -------------------------
# Itohs'/Zohm's model
D_Zohm = (D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(Z)) / 2.0
# Stap's Model
alpha_sup = 0.5
D_Staps = D_min + (D_max - D_min) / (1.0 + alpha_sup*numerix.dot(Z.grad, Z.grad))
# Flow-Shear Model
a1, a3 = 1.0, 0.5	# ASSUMES a2 = 0
D_Shear = D_min + (D_max - D_min) / (1.0 + a1*(Z)**2 + a3*numerix.dot(Z.grad, Z.grad))

# CHOOSE DIFFUSIVITY HERE!
D_choice = D_Staps

# If Diffusivity is a Cell/Face variable
Diffusivity.setValue(D_choice)
Diffusivity.equation = (ImplicitSourceTerm(1.0))

# ----------------- Boundary Conditions -------------------
"""
	Density Boundary Conditions:
	d/dx(n(0)) == n / lambda_n
	d/dx(n(L)) == -Gamma_c / Diffusivity
"""
density.faceGrad.constrain(density.faceValue / lambda_n, mesh.facesLeft)

density.faceGrad.constrain(-Gamma_c / Diffusivity.faceValue, mesh.facesRight)

"""
	Temperature Boundary Conditions:
	d/dx(T(0)) = T / lambda_T
	d/dx(T(L)) = zeta*(Gamma_c*T - q_c*(gamma - 1)) / (Diffusivity * n)
"""
temp_left = temperature.faceValue / lambda_T
temperature.faceGrad.constrain(temp_left, mesh.facesLeft)

temp_right = (zeta * (Gamma_c*temperature.faceValue - q_c*(gamma - 1.0)))\
				/ (Diffusivity.faceValue * density.faceValue)
temperature.faceGrad.constrain(temp_right, mesh.facesRight)

"""
	Paquay considered these Z Boundary Conditions:
	d/dx(Z(0)) == Z / lambda_Z
	mu*D/epsilon * d/dx(Z(L)) == 0
"""
Z.faceGrad.constrain(Z.faceValue / lambda_Z, mesh.facesLeft)
Z.faceGrad.constrain(0.0, mesh.facesRight)

# ----------------- PDE Declarations ----------------------
density.equation = (DiffusionTerm(coeff=Diffusivity, var=density) == 0)

temperature.equation = (0 ==\
		DiffusionTerm(coeff=(Diffusivity*density/zeta), var=temperature)\
		 + DiffusionTerm(coeff=-Diffusivity*temperature, var=density))

G = a + b*(Z - Z_S) + c*(Z - Z_S)**3
S_Z = ((c_n*temperature) / density**2) * density.grad[0]\
		+ (c_T / density) * temperature.grad[0] + G
Z.equation = (0 == DiffusionTerm(coeff=mu, var=Z) + S_Z)

full_equation = density.equation & temperature.equation & Z.equation

viewer = Viewer((density, temperature, Z, Diffusivity), xmin=0.0, xmax=L)

restol = 1e-1

res_D = res_n = res_T = res_Z = 1.0e10

if __name__ == '__main__':
	while res_n > restol:
		res_n = density.equation.sweep(var=density, solver=mySolver)
		viewer.plot()
		raw_input("n_res = %f" % res_n)
	while res_T > restol:
		res_T = temperature.equation.sweep(var=temperature, solver=mySolver)
		viewer.plot()
		raw_input("T_res = %f" % res_T)
#	while res_Z > restol:
#		res_Z = Z.equation.sweep(var=Z, solver=mySolver)
#		Diffusivity.setValue(D_choice)
#		viewer.plot()
#		raw_input("Z_res = %f" % res_Z)

#	density.equation.solve(var=density, solver=mySolver)
#	viewer.plot()
	raw_input("Pause")

#print density
