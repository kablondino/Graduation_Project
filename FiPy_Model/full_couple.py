from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm

from fipy.tools import numerix

nx = 1000
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

D_max = 2.0
D_min = 2.0/5.0

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
D = CellVariable(name=r"$D$", mesh=mesh, hasOld=True)


# ------------- Initial Conditions and Diffusivity---------
Z0 = Z_S*(1 - numerix.tanh((L*x - L) / 2))
Z.setValue(Z0)

# Zohm's model
#D.setValue((D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(Z)) / 2.0)
# Stap's Model
alpha_sup = 1.0
#D = (D_min + (D_max - D_min) / (1 + alpha_sup*(Z.grad.mag)**2))
# Flow-Shear Model
a1, a2, a3 = 1.0, 1.0, 1.0
D.setValue(D_min + (D_max - D_min) / (1 + a1*Z**2 + a2*Z*(Z.grad) + a3*(Z.grad)**2))

density0 = CellVariable(name=r"$n_0$", mesh=mesh, value=-(Gamma_c*lambda_n / D) * (1 + x/lambda_n))
density.setValue(density0)

temp0 = CellVariable(name=r"$T_0", mesh=mesh, value = q_c*((gamma - 1) / Gamma_c) * (1 - lambda_n / (zeta*lambda_T + lambda_n)*(1 + x/lambda_n)**(-zeta)))
temperature.setValue(temp0)


# ----------------- Boundary Conditions -------------------
# Density Boundary Conditions:
#	d/dx(n(0,t)) == n / lambda_n
#	d/dx(n(L,t)) == D / Gamma_c
density.faceGrad.constrain(density.faceValue / lambda_n, mesh.facesLeft)

density.faceGrad.constrain(-Gamma_c / D.faceValue, mesh.facesRight)

# Temperature Boundary Conditions:
#	d/dx(T(0,t)) = T / lambda_T
#	d/dx(T(L,t)) = (T*Gamma_c - (gamma - 1)*q_c) / (diffusivity * n)
temp_left = temperature[0] / lambda_T
temperature.faceGrad.constrain(temp_left, mesh.facesLeft)
# THE FOLLOWING NEEDS TO BE FIXED, since the index gets messed up
temp_right = (zeta / D[nx-1])*(temperature[nx-1]*Gamma_c - (gamma - 1)*q_c) / density[nx-1]
temperature.faceGrad.constrain(temp_right, mesh.facesRight)

# Z Boundary Conditions:
#	d/dx(Z(0,t)) == Z / lambda_Z
#	mu*D/epsilon * d/dx(Z(L,t)) == 0
#	d^2/dx^2(Z(0,t)) == 0
Z.faceGrad.constrain(Z.faceValue / lambda_Z, mesh.facesLeft)
Z.constrain(0.0, mesh.facesRight)

Z.grad.faceGrad.constrain(0.0, mesh.facesLeft)


# ----------------- PDE Declarations ----------------------
# Density Equation
density_equation = TransientTerm(var=density) == DiffusionTerm(coeff=D, var=density)

# Temperature Equation
#S_T = ((zeta+1)/zeta)*(D/density) * numerix.dot(density.grad,temperature.grad)
S_T = ConvectionTerm(coeff=((zeta+1)/zeta)*(D/density)*density.grad, var=temperature)
temp_equation = TransientTerm(var=temperature) == DiffusionTerm(coeff=5.0/0.5, var=temperature) + S_T

# Z Equation
G = a + b*(Z - Z_S) + c*(Z - Z_S)**3
S_Z = G + (c_T / density)*temperature.grad.mag + (c_T*temperature / density**2)*density.grad.mag
Z_equation = TransientTerm(coeff=epsilon, var=Z) == DiffusionTerm(coeff=mu*D/epsilon, var=Z) + S_Z

## ALTERNATE Z Equation
#S_Za = B_theta**2*((g_n_an - g_n_cx - g_n_bulk)*density.grad.mag/density + (g_T_an - g_T_cx - g_T_bulk)*temperature.grad.mag/temperature + (g_Z_an - g_Z_cx - g_Z_bulk)*Z - f_OL)
#Z_equation = TransientTerm(coeff=, var=Z) == DiffusionTerm(coeff=, var=Z) + S_Za


# Fully-Coupled Equation
full_equation = density_equation & temp_equation & Z_equation


initial_viewer = Viewer((density, temperature, Z, D))
raw_input("Pause for Initial")

if __name__ == '__main__':
	viewer = Viewer((density, temperature, Z, D), datamin=-3.0, datamax=3.0)
	for t in range(100):
		density.updateOld(); temperature.updateOld()
		Z.updateOld(); D.updateOld()
		full_equation.solve(dt=1.0e-2)
		viewer.plot()

	raw_input("End of Program. <return> to continue...")

