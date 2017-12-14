from fenics import *

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


## ----------------- Z-dependent Fluxes -------------------
# Xi expressions for ion bulk viscosity
#xi_p = Expression('(4*C*(27*pow((pow(C,2) + pow(Z,2)),2) - 7*(pow(C,2) - 3*pow(Z,2))*sqrt(nu_ai))*pow(nu_ai, 7.0/4.0)) / (189*pi*pow((pow(C,2) + pow(Z,2)),3))', degree=1, C=C, Z=Z, nu_ai=nu_ai)
xi_p_taylor = Expression('(4*(27*pow(C,2)*pow(nu_ai,7.0/4.0))) / (189*pi*pow(C,3)) - (4*(9*pow(C,2)*pow(nu_ai, 7.0/4.0) - 14*pow(nu_ai, 9.0/4.0))*pow(Z,2)) / (63*pi*pow(C,5))', degree=1, C=C, Z=Z, nu_ai=nu_ai)

#xi_t = Expression('(2*C*(135*pow((pow(C,2) + pow(Z,2)),2) - 7*(21*pow(C,4) + 3*pow(Z,2)*(-5 + 7*pow(Z,2)) + pow(C,2)*(5 + 42*pow(Z,2)))*sqrt(nu_ai))*pow(nu_ai, 7.0/4.0)) / (189*pi*pow((pow(C,2) + pow(Z,2)),3))', degree=1, C=C, Z=Z, nu_ai=nu_ai)
xi_t_taylor = Expression('-(2*(-135*pow(C,2)*pow(nu_ai, 7.0/4.0) + 35*pow(nu_ai, 9.0/4.0) + 147*pow(C,2)*pow(nu_ai, 9.0/4.0))) / (189*pi*pow(C,3)) + (2*(-45*pow(C,2)*pow(nu_ai, 7.0/4.0) + 70*pow(nu_ai, 9.0/4.0) + 49*pow(C,2)*pow(nu_ai, 9.0/4.0))*pow(Z,2)) / (63*pi*pow(C,5))', degree=1, C=C, Z=Z, nu_ai=nu_ai)

# Bulk viscosity g's
g_n_bulk = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * B_p * xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p_taylor)

g_T_bulk = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * (B_p*xi_p - B*xi_t)', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, B=B, xi_p=xi_p_taylor, xi_t=xi_t_taylor)

g_Z_bulk = Expression('pow(aspect,2)*sqrt(pi) / (4*a_m) * n * m_i * pow(v_Ti,2) * B_p*xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p_taylor)

# Ion orbit loss
#f_OL = Expression('g_OL*exp(-sqrt(nu_ai + pow(Z,4))) / sqrt(nu_ai + pow(Z,4))', degree=1, g_OL=g_OL, nu_ai=nu_ai, Z=Z)
f_OL_taylor = Expression('exp(-sqrt(nu_ai))*g_OL / sqrt(nu_ai) - (exp(-sqrt(nu_ai))*g_OL * (1 + sqrt(nu_ai)))*pow(Z,4) / (2*pow(nu_ai, 3.0/2.0))', degree=1, nu_ai=nu_ai, Z=Z, g_OL=g_OL)

## --------------------------------------------------------
# Consolidate even more fluxes
flux_sum = Expression('(g_n_an + g_n_cx + g_n_bulk)*L_n + (g_T_an + g_T_cx + g_T_bulk)*L_T + (g_Z_an + g_Z_cx + g_Z_bulk)*Z + f_OL', degree=1, Z=Z, g_n_an=g_n_an, g_n_cx=g_n_cx, g_n_bulk=g_n_bulk, L_n=L_n, g_T_an=g_T_an, g_T_cx=g_T_cx, g_T_bulk=g_T_bulk, L_T=L_T, g_Z_an=g_Z_an, g_Z_cx=g_Z_cx, g_Z_bulk=g_Z_bulk, f_OL=f_OL_taylor)

F = (m_i*mu_i / (charge*rho_pi) * dot(grad(Z), grad(v)))*dx + B_p**2*flux_sum*v*dx

# Create VTK file
vtkfile = File('first_try_answer_ind/solution.pvd')
vtkfile << Z

solve(F == 0, Z, bc)

