#from __future__ import print_function
from fenics import *
#from dolfin import *

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


"""
# Coefficients and the G function(?) ARE THESE NECESSARY? MAYBE NOT, as I might write out the entire thing with g-coefficients
dielectric_const = (n * T / v_phi) * (B_pol / B**2)

viscosity = n * mu_i * T / (v_pi * B_pol)

c_n = -n**2 * ((m_i * cx_rate * n_0) / (B_pol**2) + (e * aspect**2 * sqrt(pi) * rho_pe) / (2 * a_m * B))

c_T = -n**2 / T * ((alpha_cx * m_i * cx_rate * n_0) / (B_pol**2) - (alpha_an * e * aspect**2 * sqrt(pi) * rho_pe) / (2 * a_m * B))

G = n*T*Z * ((m_i * cx_rate) / (rho_pi * B_pol**2) - (e * aspect**2 * rho_pe * sqrt(pi)) / (2 * a_m * rho_pi * B))

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
"""

