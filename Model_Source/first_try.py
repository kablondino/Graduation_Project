from __future__ import print_function
from fenics import *

# Final time, number of steps, and step size
total_T = 2.0
num_steps = 100
dt = total_T / num_steps

# Create a mesh
nx = 50
L = 1 # Length of domain
mesh = IntervalMesh(nx, 0, L)
V = FunctionSpace(mesh, 'P', 1)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, 0, boundary)

# VALUES for const terms go here
e = 1.602e-19
B = 1
B_pol = 1
mu_i = 1
m_i = 1
cx_rate = 1
n_0 = 1
aspect = 1
a_m = 1
rho_pe = 1
rho_pi = 1

# Functions for n and T go here
n = 5*x**-2
T = 10*x**-3

# Coefficients and the G function(?)
dielectric_const = (n * T / v_phi) * (B_pol / B**2)

viscosity = n * mu_i * T / (v_pi * B_pol)

c_n = -n**2 * ((m_i * cx_rate * n_0) / (B_pol**2) + (e * aspect**2 * sqrt(pi) * rho_pe) / (2 * a_m * B))

c_T = -n**2 / T * ((alpha_cx * m_i * cx_rate * n_0) / (B_pol**2) - (alpha_an * e * aspect**2 * sqrt(pi) * rho_pe) / (2 * a_m * B))

G = n*T*Z * ((m_i * cx_rate) / (rho_pi * B_pol**2) - (e * aspect**2 * rho_pe * sqrt(pi)) / (2 * a_m * rho_pi * B))


# MODEL PDE?

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
