from fenics import *

# Final time, number of steps, and step size
total_T = 2.0
num_steps = 100
dt = total_T / num_steps

# Create a mesh
nx = 100
L = 1.0
mesh = IntervalMesh(nx, 0, L)
V = FunctionSpace(mesh, 'P', 1)
x = SpatialCoordinate(mesh)

# MODEL PDE
Z = Function(V)
v = TestFunction(V)

# Properly make the boundary conditions (all of them!)
def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, 0, boundary)

# Temporary mu_i mysterious value
mu_i = 1.0

import sys
sys.path.append("/home/kabv/Documents/Masters/Graduation_Project/Model_Source/")
from coefficients import *

Z_0 = Expression('1.0 - (exp(L*x[0] - L) - 1) / (exp(L*x[0] - L) + 1)', degree=1, L=L)
#Z_0 = Expression('1.0 - tanh((L / 2.0)*(x - 1))', degree=1, L=L, x=x)
Z_n = interpolate(Z_0, V)

## ----------------- Z-dependent Fluxes -------------------
# Xi expressions for ion bulk viscosity
xi_p = Expression('(4*C*(27*pow((pow(C,2) + pow(Z,2)),2) - 7*(pow(C,2) - 3*pow(Z,2))*sqrt(nu_ai))*pow(nu_ai, 7.0/4.0)) / (189*pi*pow((pow(C,2) + pow(Z,2)),3))', degree=1, C=C, Z=Z, nu_ai=nu_ai)

xi_t = Expression('(2*C*(135*pow((pow(C,2) + pow(Z,2)),2) - 7*(21*pow(C,4) + 3*pow(Z,2)*(-5 + 7*pow(Z,2)) + pow(C,2)*(5 + 42*pow(Z,2)))*sqrt(nu_ai))*pow(nu_ai, 7.0/4.0)) / (189*pi*pow((pow(C,2) + pow(Z,2)),3))', degree=1, C=C, Z=Z, nu_ai=nu_ai)

# Bulk viscosity g's
g_n_bulk = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * B_p * xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p)

g_T_bulk = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * (B_p*xi_p - B*xi_t)', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, B=B, xi_p=xi_p, xi_t=xi_t)

g_Z_bulk = Expression('pow(aspect,2)*sqrt(pi) / (4*a_m) * n * m_i * pow(v_Ti,2) * B_p*xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p)

# Ion orbit loss
g_OL = Expression('charge*n*nu_eff*sqrt(aspect)*rho_pi', degree=1, charge=charge, n=n, nu_eff=nu_eff, aspect=aspect, rho_pi=rho_pi)

f_OL = Expression('g_OL*exp(-sqrt(nu_ai + pow(Z,4))) / sqrt(nu_ai + pow(Z,4))', degree=1, g_OL=g_OL, nu_ai=nu_ai, Z=Z)

## --------------------------------------------------------
# Consolidate even more fluxes
#flux_sum = Expression('(g_n_an + g_n_cx + g_n_bulk)*L_n + (g_T_an + g_T_cx + g_T_bulk)*L_T + (g_Z_an + g_Z_cx + g_Z_bulk)*Z + f_OL', degree=1, Z=Z, g_n_an=g_n_an, g_n_cx=g_n_cx, g_n_bulk=g_n_bulk, L_n=L_n, g_T_an=g_T_an, g_T_cx=g_T_cx, g_T_bulk=g_T_bulk, L_T=L_T, g_Z_an=g_Z_an, g_Z_cx=g_Z_cx, g_Z_bulk=g_Z_bulk, f_OL=f_OL)

flux_sum = Expression('(g_n_an + g_n_cx)*L_n + (g_T_an + g_T_cx)*L_T + (g_Z_an + g_Z_cx)*Z', degree=1, g_n_an=g_n_an, g_n_cx=g_n_cx, L_n=L_n, g_T_an=g_T_an, g_T_cx=g_T_cx, L_T=L_T, g_Z_an=g_Z_an, g_Z_cx=g_Z_cx, Z=Z)
print flux_sum(0.5)

F = (-(m_i*B_p**2/(e*B**2)) * n*Temp*Z*v + dt*(m_i*mu_i/(e*rho_pi) * n*Temp*dot(grad(Z), grad(v))))*dx + (m_i*B_p**2/(e*B**2)*Z_n + dt*B_p**2*flux_sum)*v*dx
print F

# Create VTK file
vtkfile = File('first_try_answer/solution.pvd')

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
	Z_n = interpolate(Z, V) # FIX!

