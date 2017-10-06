# Import file with parameters, etc.
load("/home/kabv/Documents/Masters/Graduation_Project/Model_Source/parameters.sage")

## ------- MATHEMATICA definitions of xi functions --------
# Result of integral from Taylor expanded integrand of x^2*exp(-x)*arctan(...)
xi_p_math(T,Z) = (4*C*(27*(C^2 + Z^2)^2 - 7*(C^2 - 3*Z^2)*sqrt(nu_ai))*nu_ai^(7/4)) / (189*pi*(C^2 + Z^2)^3)

xi_t_math(T,Z) = (2*C*(135*(C^2 + Z^2)^2 - 7*(21*C^4 + 3*Z^2*(-5 + 7*Z^2) + C^2*(5 + 42*Z^2))*sqrt(nu_ai))*nu_ai^(7/4)) / (189*pi*(C^2 + Z^2)^3)

# Taylor expand of the above to fit the target model's form
xi_p_math_exp(T,Z) = (4*(27*C^2*nu_ai^(7/4) - 7*nu_ai^(9/4))) / (189*pi*C^3) - (4*(9*C^2*nu_ai^(7/4) - 14*nu_ai^(9/4)) * Z^2) / (63*pi*C^5)

xi_t_math_exp(T,Z) = (2*(135*C^2*nu_ai^(7/4) - 35*nu_ai^(9/4) - 147*C^2*nu_ai^(9/4))) / (189*pi*C^3) + (2*(-45*C^2*nu_ai(7/4) + 70*nu_ai^(9/4) + 49*C^2*nu_ai^(9/4)) * Z^2) / (63*pi*C^5)


## ------- Sage Numerical Integral ------------------------
p_integrand(x, T, Z) = x^2 * exp(-x) * arctan((2*C^2*sqrt(x)) / (C^2 + Z^2 - x))
t_integrand(x, T, Z) = p_integrand * (5/2 * x)


## ------- Dump to file -----------------------------------
with open('/home/kabv/Documents/Masters/Graduation_Project/Model_Source/bulk_visc.dat', 'w') as the_file:
# Print labels
	the_file.write("T Z xi_p_math xi_t_math xi_p_math_expand xi_t_math_expand p_numerical p_numerical_error t_numerical t_numerical_error\n")
	for i in numpy.arange(1e2, 2e3 + 25, 25):
		for j in numpy.arange(-3.0, 3.1, 0.1):
			T = i; Z = j
			the_file.write("%s %s %s %s %s %s %s %s\n" % (i, j, N(xi_p_math(i, j)), N(xi_t_math(i, j)), N(xi_p_math_exp(i, j)), N(xi_t_math_exp(i, j)), numerical_integral(N(1/pi)*p_integrand(x), 0, sqrt(nu_ai(i)), params=[T,Z]), numerical_integral(N(1/pi)*t_integrand(x), 0, sqrt(nu_ai(i)), params=[T,Z])))

the_file.close()

