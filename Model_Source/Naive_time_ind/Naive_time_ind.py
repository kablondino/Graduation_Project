from fenics import *
import numpy

def solver(Z_D, Nx, A, B, C, D, mu, degree=1):
	# Create mesh
	mesh = IntervalMesh(Nx, 0, 1)
	V = FunctionSpace(mesh, 'P', degree)
	x = SpatialCoordinate(mesh)

	def boundary(x, on_boundary):
		return on_boundary

	bc = DirichletBC(V, Z_D, boundary)

	# Define variational problem
	v = TestFunction(V)
	Z = Function(V)
	F = (mu*dot(grad(Z), grad(v)))*dx - (A + B*Z + C*pow(Z,2) + D*pow(Z,3))*v*dx

	# Solve
	solve(F == 0, Z, bc, solver_parameters={"newton_solver":{"maximum_iterations":10000}})

	return Z


def run_solver():
	"""Run solver to compute and post-process it"""

	Z_D = Expression('1 + tanh((1.0/2.0)*(x[0] - 1))',degree=1)
	Z = solver(Z_D, 100, 0, 0, 0, 1, 1, 1)

	plot(Z)
	plot(Z.function_space().mesh())

	vtkfile = File('naive_time_ind.pvd')
	vtkfile << Z


def test_solver():
	"""Test somehow"""

	tol = 1e-10
	Z_D = Expression('1 + tanh((1.0/2.0)*(x[0] - 1))',degree=1)

	# Iterate over mesh sizes and degrees
	for Nx in [3, 10, 20, 50, 100, 500]:
		for degree in 1, 2, 3:
			print 'Solving on a 1 x (%d) mesh with P%d elements.' % (Nx, degree)

			Z = solver(Z_D, Nx, degree)

			# Extract the mesh
			mesh = Z.function_space().mesh()
			vertex_values_Z_D = Z_D.compute_vertex_values(mesh)
			vertex_values_Z = Z.compute_vertex_values(mesh)
			error_max = numpy.max(numpy.abs(vertex_values_Z_D - vertex_values_Z))

			# Check maximum error
			msg = 'Error_max = %g' % error_max
			assert error_max < tol, msg


if __name__ == '__main__':
	run_solver()
	interactive()


# Print out raw data to output
#Z_nodal_values = Z.vector()
#Z_array = Z_nodal_values.array()
#
#for i in range(len(Z_array)):
#	print '%g %g' % (mesh.coordinates()[i], Z_array[i])

