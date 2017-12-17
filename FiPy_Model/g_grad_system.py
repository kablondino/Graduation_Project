from fipy import TransientTerm, DiffusionTerm, Viewer, ConvectionTerm, ImplicitSourceTerm

from parameters import *
from coeffs import *
#q_c = -1000.0		# HOT SWAP heat and particle
#Gamma_c = -10.0	# fluxes from the core.

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
D_choice = D_Staps
Diffusivity.setValue(D_choice)


density0 = -(Gamma_c*lambda_n / Diffusivity) * (1 + x/lambda_n)
density.setValue(1.0e19*density0)

temp0 = q_c*((gamma - 1) / Gamma_c) * (1 - lambda_n / (zeta*lambda_T + lambda_n)*(1 + x/lambda_n)**(-zeta))
temperature.setValue(10.0*temp0)

# ----------------- Printing for testing ------------------
#print numerix.dot([1.0/density,], density.grad)

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
# Diffusivity Equation
diffusivity_equation = ImplicitSourceTerm(coeff=1.0, var=Diffusivity) == D_choice

# Density Equation
density_equation = TransientTerm(var=density) == DiffusionTerm(coeff=Diffusivity, var=density)

# Temperature Equation
S_T = ((zeta+1)/zeta)*(Diffusivity/density) * numerix.dot(density.grad,temperature.grad)
temp_equation = TransientTerm(var=temperature) == DiffusionTerm(coeff=Diffusivity/zeta, var=temperature) + S_T

# Z Equation
S_Z = B_theta**2 * ( (g_n_an - g_n_cx - g_n_bulk)*numerix.dot([1.0/density,], density.grad) + (g_T_an - g_T_cx - g_T_bulk)*numerix.dot([1.0/temperature,], temperature.grad) + (g_Z_an - g_Z_cx - g_Z_bulk)*Z - f_OL )
Z_equation = TransientTerm(coeff=(m_i * density*temperature/ (charge)*(B_theta/B)**2), var=Z) == DiffusionTerm(coeff=m_i*mu * density*temperature/ (charge), var=Z)# + S_Z
Z_equation = TransientTerm(coeff=epsilon, var=Z) == DiffusionTerm(coeff=mu*Diffusivity, var=Z)

#print Z_equation

# Fully-Coupled Equation
full_equation = density_equation & temp_equation & Z_equation & diffusivity_equation

#initial_viewer = Viewer((density, temperature, Z, Diffusivity))
#raw_input("Pause for Initial")
print density, temperature
if __name__ == '__main__':
	viewer = Viewer((density / 1.0e19, temperature / 10.0, -Z, Diffusivity), legend='best')#, ylog=True, datamin = 1.0e-4)
#	viewer_rho_pi = Viewer((rho_pi))
#	viewer_nu_ai = Viewer((nu_ai), name="TEST", title=None, legend=None)
	for t in range(100):
#		print t, q_c, Gamma_c
#		q_c = q_c - 1.0
		Gamma_c = Gamma_c - 1.0
		density.updateOld(); temperature.updateOld()
		Z.updateOld(); Diffusivity.updateOld()
		full_equation.solve(dt=0.01)
		viewer.plot()

	raw_input("End of Program. <return> to continue...")

