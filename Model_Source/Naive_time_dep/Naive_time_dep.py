from fenics import *
from dolfin import *

# Create mesh
nx = 100
mesh = IntervalMesh(nx, 0, 1)
V = FunctionSpace(mesh, 'P', 1)
x = SpatialCoordinate(mesh)

# Boundary condition
Z_D = Constant(0)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, Z_D, boundary)

# Initial value
Z_0 = Constant(0)

# Define Constants
mu = 1
A = 1
B = 1
C = 1
D = 1

# Define variational problem
v = TestFunction(V)
Z = Function(V)
f = A + B*Z + C*Z**2 + D*Z**3
F = (mu*dot(grad(Z), grad(v)))*dx - f*v*dx

vtkfile = File('naive_time_ind.pvd')

# Compute Solution
solve(F == 0, Z, bc)

vtkfile << Z

plot(Z)

# Print out raw data to output
Z_nodal_values = Z.vector()
Z_array = Z_nodal_values.array()

for i in range(len(Z_array)):
	    print '%g %g' % (mesh.coordinates()[i], Z_array[i])

