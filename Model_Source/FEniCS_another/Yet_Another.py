from fenics import *
import sys
sys.path.append("/home/kabv/Documents/Masters/Graduation_Project/Model_Source/")
from parameters import *
import coefficients as coef
#from coefficients import *

# Create a mesh
nx = 100
L = 1.0
mesh = IntervalMesh(nx, 0, L)
V = FunctionSpace(mesh, 'P', 1)

# Set up independent and dependent variables
x = SpatialCoordinate(mesh)
z = TrialFunction(V)
v = TestFunction(V)
Z = Function(V)

Z_D = Expression('1 - tanh(L/2.0*(x[0] - 1))', degree=1, L=L)

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, Z_D, boundary)

# Write the coefficients over from the coefficients file
g_n_an, g_T_an, g_Z_an = coef.g_n_an, coef.g_T_an, coef.g_Z_an
g_n_cx, g_T_cx, g_Z_cx = coef.g_n_cx, coef.g_T_cx, coef.g_Z_cx
g_n_bulk, g_T_bulk, g_Z_bulk = coef.g_n_bulk, coef.g_T_bulk, coef.g_Z_bulk
f_OL = coef.f_OL

#f = Expression('pow(B_p,2)*((g_n_an - g_n_cx)*L_n + (g_T_an - g_T_cx)*L_T)', degree=2, B_p=B_p, g_n_an=g_n_an, g_n_cx=g_n_cx, L_n=L_n, g_T_an=g_T_an, g_T_cx=g_T_cx, L_T=L_T)
f = Expression('(g_n_an - g_n_cx)*L_n + (g_T_an - g_T_cx)*L_T', degree=1, g_n_an=g_n_an, g_n_cx=g_n_cx, L_n=L_n, g_T_an=g_T_an, g_T_cx=g_T_cx, L_T=L_T)
g = Expression('(-g_Z_an + g_Z_cx + g_Z_bulk)*Z + g_n_bulk*L_n + g_T_bulk*L_T - f_OL', degree=1, g_n_bulk=g_n_bulk, g_T_bulk=g_T_bulk, g_Z_bulk=g_Z_bulk, L_n=L_n, L_T=L_T, g_Z_an=g_Z_an, g_Z_cx=g_Z_cx, f_OL=f_OL, Z=Z)
print f

#Z = interpolate(Expression('1 - tanh(L/2.0*(x[0] - 1))', degree=1, L=L), V)
F = Z*v*dx + (-m_i / (charge*B_p**2))*rho_pi**(-1) * dot(grad(Z), grad(v))*dx + f*v*dx

J = derivative(F, Z, v)
print J

#problem = NonlinearVariationalProblem(F, Z, bc, J)

#solver = NonlinearVariationalSolver(problem)

#solver.solve()

#vtkfile = File('solution.pvd')

solve(F == 0, Z, bc, solver_parameters={"newton_solver":{"maximum_iterations":1000}})

#vtkfile << Z

