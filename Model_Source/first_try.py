from __future__ import print_function
from fenics import *

# Final time, number of steps, and step size
total_T = 2.0
num_steps = 100
dt = total_T / num_steps

# Create a mesh
nx = 50
mesh = IntervalMesh(nx, 0, 1)
V = FunctionSpace(mesh, 'P', 1)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, 0, boundary)

# VALUES for const terms go here
e = 1.602e-19
B = 
B_pol = 
v_phi = 
mu_i = 
m_i = 
cx_rate = 
n_0 = 
aspect = 
a_m = 
rho_pe = 
rho_pi = 

# Functions for n and T go here


# Coefficients and the G function(?)
dielectric_const = (n_i * T / v_phi) * (B_pol / B**2)

viscosity = n_i * mu_i * T / (v_pi * B_pol)

c_n = -n**2 * ((m_i * cx_rate * n_0) / (B_pol**2) + (e * aspect**2 * sqrt(pi) * rho_pe) / (2 * a_m * B))

c_T = -n**2 / T * ((alpha_cx * m_i * cx_rate * n_0) / (B_pol**2) - (alpha_an * e * aspect**2 * sqrt(pi) * rho_pe) / (2 * a_m * B))

G = n*T*Z * ((m_i * cx_rate) / (rho_pi * B_pol**2) - (e * aspect**2 * rho_pe * sqrt(pi)) / (2 * a_m * rho_pi * B))


# Create VTK file
vtkfile = File('first_try/solution.pvd')

# Time step
Z = Function(V)
t = 0

for n in range(num_steps):
	# Update time
	t += dt

	# Compute solution
	solve(a == L, u, bc) # FIX!

	# Save to file and plot
	vtkfile << (Z, t) #
	plot(u)

	# Update previous solution
	Z_n = interpolate( INITIAL_Z, V) # FIX!

# Hold plot, if available
interactive()
