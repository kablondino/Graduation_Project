from fenics import *
from dolfin import *

# Create mesh
nx = 100
mesh = IntervalMesh(nx, 0, 1)
V = FunctionSpace(mesh, 'P', 1)
x = SpatialCoordinate(mesh)

# Time dependence
Ttotal = 2.0
num_steps = 50
dt = Ttotal / num_steps

# Boundary condition
Z_D = Constant(0)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, Z_D, boundary)

# Initial value
#Z_0 = Expression('pow(x,2)')
Z_0 = Constant(0)
Z_n = interpolate(Z_0, V)

# Define Constants
epsilon = 1
mu = 1
A = 1
B = 1
C = 1
D = 1

# Define variational problem
v = TestFunction(V)
Z = Function(V)
f = A + B*Z + C*Z**2 + D*Z**3
F = (epsilon*Z*v + dt*mu*dot(grad(Z), grad(v)))*dx - (Z_n + dt*f)*v*dx

a, L = lhs(F), rhs(F)

vtkfile = File('data/naive_time_dep.pvd')

t = 0
for n in range(num_steps):
	# Update time
	t += dt

	# Compute Solution
	solve(F == 0, Z, bc, solver_parameters={"newton_solver":{"maximum_iterations":1000}})

	# Write to vtu file
	vtkfile << (Z, t)

	# Print out raw data to output
	Z_nodal_values = Z.vector()
	Z_array = Z_nodal_values.array()

	for i in range(len(Z_array)):
		print '%g %g %g' % (t, mesh.coordinates()[i], Z_array[i])

	plot(Z)

	# Reassign old value to new
	Z_n.assign(Z)

