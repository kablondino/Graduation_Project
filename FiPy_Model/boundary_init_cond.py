"""
	Declares the initial conditions, including whether it
	starts in L-- or H--mode.

	Also declares the boundary conditions of the state
	variables, as a function.
"""

from variable_decl import *


# ---------------- Set Initial Conditions -----------------
# Initial conditions for L--mode
if config.initial_H_mode == False:
	density.setValue(1.5*x / 4.0 + 0.5)

	temperature.setValue(250.0*(0.0125*x**2 + 0.2*x + 1.2))

#	Z.setValue(0.0)
	Z.setValue(-3.0 / (1.0 + numerix.exp(12.0*(x - 1.75))))

# Initial conditions for H--mode
elif config.initial_H_mode == True:
	density.setValue((0.5/1.5)*x + 0.5)
	density.setValue(3.0*x - 3.5, where = x > 1.5)
	density.setValue((0.5/1.5)*x + (11.0/6.0),\
			where = x > 2.0)

	temperature.setValue(0.2*x + 1.2)
	temperature.setValue(1.8*x - 1.2, where = x > 1.5)
	temperature.setValue(0.2*x + 2.0, where = x > 2.0)

	Z.setValue(-3.0 / (1.0 + numerix.exp(12.0*(x - 1.75))))

else:
	print "Something went horribly wrong in choosing initial conditions."


# ----------------- Set Diffusivity Model -----------------
## It defaults to Stap's version
# Itohs'/Zohm's model
if config.D_choice.lower() == "d_zohm":
	D_choice_local = (D_max + D_min) / 2.0\
			+ ((D_max - D_min)*numerix.tanh(Z)) / 2.0

# Stap's Model
elif config.D_choice.lower() == "d_staps":
	D_choice_local = D_min + (D_max - D_min) / (1.0 + config.alpha_sup\
			* numerix.sign(Z.grad[0])*(abs(Z.grad[0]))**config.beta)

# Flow-Shear Model
elif config.D_choice.lower() == "d_shear" or\
		config.D_choice.lower() == "d_flow_shear":
	a1, a2, a3 = 1.0, 0.0, 1.0
	D_choice_local = D_min + (D_max - D_min)\
			/ (1.0 + a1*(Z)**2 + a2*Z*Z.grad[0]\
			+ a3*numerix.dot(Z.grad, Z.grad))

else:
	print "Something went horribly wrong in choosing the Diffusivity model."

Diffusivity.setValue(D_choice_local)
print "The diffusivity model is set to " + str(config.D_choice)


# ------ Old definitions, which requires Diffusivity to ALREADY be set -------
#	density.setValue(-(config.Gamma_c*lambda_n / Diffusivity)\
#			* (1.0 + x/lambda_n))
#	temperature.setValue(config.q_c\
#			*((gamma - 1.0) / config.Gamma_c)\
#			*(1.0 - lambda_n / (zeta*lambda_T + lambda_n)\
#			*(1.0 + x/lambda_n)**(-zeta)))

#	Z.setValue(Z_S*(1.0 - numerix.tanh(\
#			(config.L*x - config.L) / 2.0)))	# OLD, by Staps


# ----------------- Boundary Conditions ------------------- #
def set_boundary_values(AGamma_c, Aq_c):
	the_Gamma_c = PhysicalField(value=AGamma_c, unit="1 / (m**2 * s)")
	the_q_c = PhysicalField(value=Aq_c, unit="eV / (m**2 * s)")

	"""
		Density Boundary Conditions:
		d/dx(n(0)) == n / lambda_n
		d/dx(n(L)) == -Gamma_c / Diffusivity
	"""
	density.faceGrad.constrain(density.faceValue / lambda_n, mesh.facesLeft)
	density.faceGrad.constrain(\
			-the_Gamma_c / Diffusivity.faceValue, mesh.facesRight)

	"""
		Temperature Boundary Conditions:
		d/dx(T(0)) = T / lambda_T
		d/dx(T(L)) = zeta*(Gamma_c*T - q_c*(gamma - 1)) / (Diffusivity * n)
	"""
	temp_left = temperature.faceValue / lambda_T
	temperature.faceGrad.constrain(temp_left, mesh.facesLeft)
	temp_right = temperature.faceGrad.constrain(\
			(zeta * (the_Gamma_c*temperature.faceValue -\
			the_q_c*(gamma - 1.0))) / (Diffusivity.faceValue *\
			density.faceValue), mesh.facesRight)

	"""
		Paquay considered these Z Boundary Conditions:
		d/dx(Z(0)) == Z / lambda_Z
		mu*D/epsilon * d/dx(Z(L)) == 0
	"""
#	Z.faceGrad.constrain(Z.faceValue / lambda_Z, mesh.facesLeft)
	Z.faceGrad.constrain(0.0, mesh.facesRight)


set_boundary_values(config.Gamma_c, config.q_c)

