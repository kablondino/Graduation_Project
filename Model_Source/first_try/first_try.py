from fenics import *

# Final time, number of steps, and step size
total_T = 2.0
num_steps = 100
dt = total_T / num_steps

# Create a mesh
nx = 100
mesh = IntervalMesh(nx, 0, 1.0)
V = FunctionSpace(mesh, 'P', 1)
x = SpatialCoordinate(mesh)

# MODEL PDE
Z = Function(V)
v = TestFunction(V)


def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, 0, boundary)

# Temporary mu_i mysterious value
mu_i = 1

import sys
sys.path.append("/home/kabv/Documents/Masters/Graduation_Project/Model_Source/")
from coefficients import *

Z_0 = Constant(0)
Z_n = interpolate(Z_0, V)

# Consolidate even more fluxes
flux_sum = Expression('(g_n_an + g_n_cx + g_n_bulk)*L_n + (g_T_an + g_T_cx + g_T_bulk)*L_T + (g_Z_an + g_Z_cx + g_Z_bulk)*Z + f_OL', degree=1, g_n_an=g_n_an, g_n_cx=g_n_cx, g_n_bulk=g_n_bulk, L_n=L_n, g_T_an=g_T_an, g_T_cx=g_T_cx, g_T_bulk=g_T_bulk, L_T=L_T, g_Z_an=g_Z_an, g_Z_cx=g_Z_cx, g_Z_bulk=g_Z_bulk, f_OL=f_OL)
print flux_sum(0.5)

F = (-(m_i*B_p**2/(e*B**2)) * n*Temp*Z*v + dt*(m_i*mu_i/(e*rho_pi) * n*Temp*dot(grad(Z), grad(v))))*dx + (m_i*B_p**2/(e*B**2)*Z_n + dt*B_p**2)*v*dx

"""
solve(F == 0, Z, bc)

# Create VTK file
vtkfile = File('first_try/solution.pvd')

# Time step
t = 0
for n in range(num_steps):
	# Update time
	t += dt

	# Compute solution
	solve(F == 0, Z, bc) # FIX!

	# Save to file and plot
	vtkfile << (Z, t) #

	# Update previous solution
	Z_n = interpolate( Z_0, V ) # FIX!
"""
