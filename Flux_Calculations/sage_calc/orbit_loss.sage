# Import file with parameters, etc.
load("/home/kabv/Documents/Masters/Graduation_Project/Flux_Calculations/parameters.sage")

g_ol(T) = e*n*nu_eff*sqrt(aspect)*rho_pi
orbit_flux(T,Z) = g_ol*exp(-sqrt(nu_ai + Z^4)) / sqrt(nu_ai + Z^4)

#plot3d(orbit_flux, (T,1e2,2e3), (Z,-3,3), adaptive=True).show()

orbit_flux_exp(T,Z) = exp(-sqrt(nu_ai)) / sqrt(nu_ai)# - (exp(-sqrt(nu_ai))*g_ol*(1 + sqrt(nu_ai)))*Z^4 / (2*nu_ai^(3/2))
#plot3d(orbit_flux_exp, (T,1e2,2e3), (Z,-3,3), adaptive=True).show()

orbit_flux_exp2(T,Z) = exp(-sqrt(nu_ai)) / sqrt(nu_ai) - (exp(-sqrt(nu_ai))*g_ol*(1 + sqrt(nu_ai)))*Z^4 / (2*nu_ai^(3/2))

with open('/home/kabv/Documents/Masters/Graduation_Project/Flux_Calculations/orbit_loss.dat', 'w') as the_file:
	# Print labels
	the_file.write("T Z g_ol orbit_flux orbit_flux_exp\n")
	
	for i in numpy.arange(1e2,2e3 + 25, 25):
		for j in numpy.arange(-3.0, 3.1, 0.1):
			the_file.write("%s %s %s %s %s %s\n" % (i, j, N(g_ol(i)), N(orbit_flux(i, j)), N(orbit_flux_exp(i, j)), N(orbit_flux_exp2(i, j))))

the_file.close()

