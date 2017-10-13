#from __future__ import print_function
from fenics import *
from dolfin import *

# Final time, number of steps, and step size
total_T = 2.0
num_steps = 100
dt = total_T / num_steps

# Create a mesh
nx = 100
mesh = IntervalMesh(nx, 0, 1.0)
V = FunctionSpace(mesh, 'P', 1)
x = SpatialCoordinate(mesh)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, 0, boundary)

# VALUES for const terms go here
from parameters import *

# Xi expressions for ion bulk viscosity
#xi_pol = (108*aspect**3*nu_ai**(3.0/2)*nu_ei**2*nu_ii - 28*nu_ii**3) / (189*pi*aspect**(9.0/2)*nu_ai**(3.0/4)*nu_ei**3) + ((-36*aspect**2*nu_ai**(3.0/2.0)*nu_ei**2*nu_ii**3 + 56*nu_ii**5) / (63*pi*aspect**(15.0/2.0)*nu_ai**(11.0/4.0)*nu_ei**5)) * Z**2

#xi_tor = (270*aspect**2*nu_ai**(3.0/2.0)*nu_ei**2*nu_ii - 294*aspect**3*nu_ai**2*nu_ei**2*nu_ii - 70*nu_ii**3) / (189*pi*aspect) + ((-90*aspect**3*nu_ai**(3.0/2.0)*nu_ei**2*nu_ii**3 + 98*aspect**3*nu_ai**2*nu_ei**2*nu_ii**3 + 140*nu_ii**5) / (63*pi*aspect**(15.0/2.0)*nu_ai**(11.0/4.0)*nu_ei**5)) * Z**2

## Anomalous diffusion coefficient
D_an = lambda X: (aspect**2 * sqrt(pi)) / (2*a_m) * (rho_pe(X) * T(X)) / B

## ----------------- g_n coefficients -----------
#g_n_pi_paral = aspect**2*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * (v_T_i)**2 * B_pol * xi_pol
g_n_an = lambda X: -e*n(X)*D_an(X)
g_n_cx = lambda X: -m_i*n_0(X)*R_cx *n(X)*T(X) / B_p**2

## ----------------- g_T coefficients -----------
#g_T_pi_paral = aspect**2*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * (v_T_i)**2 * (B_pol*xi_pol - B*xi_tor)
g_T_cx = lambda X: a_cx*g_n_cx(X)
g_T_an = lambda X: g_n_cx(X)*a_an

## ----------------- g_Z coefficients -----------
#g_Z_pi_paral = aspect**2*sqrt(pi) / (4*a_m) * n * m_i * (v_T_i)**2 * B_pol*xi_pol
g_Z_cx = lambda X: m_i*n_0(X)*R_cx *n(X)*T(X) / (rho_pi(X)*B_p**2)
g_Z_an = lambda X: -e*n(X)*D_an(X) / rho_pi(X)

test = 0.5

print D_an(test), g_n_an(test), g_n_cx(test), g_T_an(test), g_T_cx(test), g_Z_an(test), g_Z_cx(test)

"""
# Coefficients and the G function(?) ARE THESE NECESSARY? MAYBE NOT, as I might write out the entire thing with g-coefficients
dielectric_const = (n * T / v_phi) * (B_pol / B**2)

viscosity = n * mu_i * T / (v_pi * B_pol)

c_n = -n**2 * ((m_i * cx_rate * n_0) / (B_pol**2) + (e * aspect**2 * sqrt(pi) * rho_pe) / (2 * a_m * B))

c_T = -n**2 / T * ((alpha_cx * m_i * cx_rate * n_0) / (B_pol**2) - (alpha_an * e * aspect**2 * sqrt(pi) * rho_pe) / (2 * a_m * B))

G = n*T*Z * ((m_i * cx_rate) / (rho_pi * B_pol**2) - (e * aspect**2 * rho_pe * sqrt(pi)) / (2 * a_m * rho_pi * B))
"""

# TEMPORARY mu_i
mu_i = 1

# MODEL PDE?
Z = Function(V)
v = TestFunction(V)
F = -(m_i * mu_i / e) * n * T / rho_pi * dot(grad(Z), grad(v))*dx + B_p**2*((g_n_an + g_n_cx)*L_n + (g_T_an + g_T_cx)*L_T + (g_Z_an + g_Z_cx)*Z)*v*dx

# Create VTK file
vtkfile = File('first_try/solution.pvd')

# Time step
t = 0

for n in range(num_steps):
	# Update time
	t += dt

	# Compute solution
	solve(F == 0, u, bc) # FIX!

	# Save to file and plot
	vtkfile << (Z, t) #
	plot(u)

	# Update previous solution
	Z_n = interpolate( INITIAL_Z, V ) # FIX!

