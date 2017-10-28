"""
	This file imports the numerical parameters (as well as the density and
	temperature profiles defined in parameters.py) in order to calculated the
	coefficients used, including the Taylor expanded versions (in Z -> 0).

	In addition, it defines the mesh (size and shape), spatial coordinate,
	V-function space, and the dependent variable Z.
"""

from parameters import *

# Create a mesh
nx = 100
L = 1.0
mesh = IntervalMesh(nx, 0, L)
#ny = 600
#mesh = RectangleMesh(Point(0.0,-3.0), Point(L,3.0), nx, ny)
V = FunctionSpace(mesh, 'P', 1)

x = SpatialCoordinate(mesh)
Z = Function(V)

# Xi expressions for ion bulk viscosity
def python_xi_p(X, aZ): return (4*C(X)*(27*((C(X))**2 + aZ**2)**2 - 7*((C(X))**2 - 3*aZ**2)*(nu_ai(X))**(1.0/2))*(nu_ai(X))**(7.0/4)) / (189*pi*((C(X))**2 + aZ**2)**3)

xi_p = Expression('(4*C*(27*pow((pow(C,2) + pow(x[1],2)),2) - 7*(pow(C,2) - 3*pow(x[1],2))*sqrt(nu_ai))*pow(nu_ai, 7.0/4.0)) / (189*pi*pow((pow(C,2) + pow(x[1],2)),3))', degree=1, C=C, nu_ai=nu_ai)

## Taylor expanded versions
def python_xi_p_taylor(X, aZ): return (4*(27*(C(X))**2*(nu_ai(X))**(7.0/4.0))) / (189*pi*(C(X))**3) - (4*(9*(C(X))**2*(nu_ai(X))**(7.0/4.0) - 14*(nu_ai(X))**(9.0/4.0))*aZ**2) / (63*pi*(C(X))**5)
xi_p_taylor = Expression('(4*(27*pow(C,2)*pow(nu_ai,7.0/4.0))) / (189*pi*pow(C,3)) - (4*(9*pow(C,2)*pow(nu_ai, 7.0/4.0) - 14*pow(nu_ai, 9.0/4.0))*pow(x[1],2)) / (63*pi*pow(C,5))', degree=1, C=C, nu_ai=nu_ai)

def python_xi_t(X, aZ): return (2*C(X)*(135*((C(X))**2 + aZ**2)**2 - 7*(21*(C(X))**4 + 3*aZ**2*(-5 + 7*aZ**2) + (C(X))**2*(5 + 42*aZ**2))*(nu_ai(X))**(1.0/2))*(nu_ai(X))**(7.0/4)) / (189*pi*((C(X))**2 + aZ**2)**3)
xi_t = Expression('(2*C*(135*pow((pow(C,2) + pow(x[1],2)),2) - 7*(21*pow(C,4) + 3*pow(x[1],2)*(-5 + 7*pow(x[1],2)) + pow(C,2)*(5 + 42*pow(x[1],2)))*sqrt(nu_ai))*pow(nu_ai, 7.0/4.0)) / (189*pi*pow((pow(C,2) + pow(x[1],2)),3))', degree=1, C=C, nu_ai=nu_ai)

## Taylor Expanded versions
def python_xi_t_taylor(X, aZ): return -(2*(-135*(C(X)**2)*(nu_ai(X))**(7.0/4.0) + 35*(nu_ai(X))**(9.0/4.0) + 147*(C(X))**2*(nu_ai(X))**(9.0/4.0))) / (189*pi*(C(X))**3) + (2*(-45*(C(X))**2*(nu_ai(X))**(7.0/4.0) + 70*(nu_ai(X))**(9.0/4.0) + 49*(C(X))**2*(nu_ai(X))**(9.0/4.0))*aZ**2) / (63*pi*(C(X))**5)
xi_t_taylor = Expression('-(2*(-135*pow(C,2)*pow(nu_ai, 7.0/4.0) + 35*pow(nu_ai, 9.0/4.0) + 147*pow(C,2)*pow(nu_ai, 9.0/4.0))) / (189*pi*pow(C,3)) + (2*(-45*pow(C,2)*pow(nu_ai, 7.0/4.0) + 70*pow(nu_ai, 9.0/4.0) + 49*pow(C,2)*pow(nu_ai, 9.0/4.0))*pow(x[1],2)) / (63*pi*pow(C,5))', degree=1, C=C, nu_ai=nu_ai)

## Anomalous diffusion coefficient
def python_D_an(X): return (aspect**2 * pi**(1.0/2)) / (2*a_m) * (rho_pe(X) * Temp(X)) / B
D_an = Expression('pow(aspect,2)*sqrt(pi) / (2*a_m) * (rho_pe * Temp) / B', degree=1, aspect=aspect, a_m=a_m, rho_pe=rho_pe, Temp=Temp, B=B)


## ----------------- g_n coefficients ---------------------
def python_g_n_bulk(X, aZ): return aspect**2*(pi)**(1.0/2) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * B_p * xi_p(X, aZ)
g_n_bulk = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * B_p * xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p)

# TAYLOR EXPANDED VERSIONS
def python_g_n_bulk_taylor(X, aZ): return aspect**2*pi**(1.0/2.0) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * B_p * python_xi_p_taylor(X, aZ)
g_n_bulk_taylor = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * B_p * xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p_taylor)

def python_g_n_an(X): return -charge*n(X)*D_an(X)
g_n_an = Expression('-charge*n*D_an', degree=1, charge=charge, n=n, D_an=D_an)

def python_g_n_cx(X): return -m_i*n_0(X)*R_cx *n(X)*Temp(X) / B_p**2
g_n_cx = Expression('-m_i*n_0*R_cx *n*Temp / pow(B_p,2)', degree=1, m_i=m_i, n_0=n_0, R_cx=R_cx, n=n, Temp=Temp, B_p=B_p)


## ----------------- g_T coefficients ---------------------
def python_g_T_bulk(X, aZ): return aspect**2*(pi)**(1.0/2) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * (B_p*xi_p(X, aZ) - B*xi_t(X, aZ))
g_T_bulk = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * (B_p*xi_p - B*xi_t)', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, B=B, xi_p=xi_p, xi_t=xi_t)

# TAYLOR EXPANDED VERSIONS
def python_g_T_bulk_taylor(X, aZ): return aspect**2*pi**(1.0/2.0) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * (B_p*python_xi_p_taylor(X, aZ) - B*python_xi_t_taylor(X, aZ))
g_T_bulk_taylor = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * (B_p*xi_p - B*xi_t)', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, B=B, xi_p=xi_p_taylor, xi_t=xi_t_taylor)

def python_g_T_an(X): return g_n_cx(X)*a_an
g_T_an = Expression('g_n_an*a_an', degree=1, g_n_an=g_n_an, a_an=a_an)

def python_g_T_cx(X): return a_cx*g_n_cx(X)
g_T_cx = Expression('a_cx*g_n_cx', degree=1, a_cx=a_cx, g_n_cx=g_n_cx)


## ----------------- g_Z coefficients ---------------------
def python_g_Z_bulk(X, aZ): return aspect**2*(pi)**(1.0/2) / (4*a_m) * n(X) * m_i * (v_Ti(X))**2 * B_p*xi_p(X, aZ)
g_Z_bulk = Expression('pow(aspect,2)*sqrt(pi) / (4*a_m) * n * m_i * pow(v_Ti,2) * B_p*xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p)

# TAYLOR EXPANDED VERSIONS
def python_g_Z_bulk_taylor(X, aZ): return aspect**2*pi**(1.0/2.0) / (4*a_m) * n(X) * m_i * (v_Ti(X))**2 * B_p*python_xi_p_taylor(X, aZ)
g_Z_bulk_taylor = Expression('pow(aspect,2)*sqrt(pi) / (4*a_m) * n * m_i * pow(v_Ti,2) * B_p*xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p_taylor)

def python_g_Z_an(X): return -charge*n(X)*D_an(X) / rho_pi(X)
g_Z_an = Expression('-charge*n*D_an / rho_pi', degree=1, charge=charge, n=n, D_an=D_an, rho_pi=rho_pi)

def python_g_Z_cx(X): return m_i*n_0(X)*R_cx *n(X)*Temp(X) / (rho_pi(X)*B_p**2)
g_Z_cx = Expression('m_i*n_0*R_cx * n*Temp / (rho_pi*pow(B_p,2))', degree=1, m_i=m_i, n_0=n_0, R_cx=R_cx, n=n, Temp=Temp, rho_pi=rho_pi, B_p=B_p)


## ----------------- f_OL -------------------- ------------
def python_g_OL(X): return charge*n(X)*nu_eff(X)*(aspect)**(1.0/2)*rho_pi(X)
g_OL = Expression('charge*n*nu_eff*sqrt(aspect)*rho_pi', degree=1, charge=charge, n=n, nu_eff=nu_eff, aspect=aspect, rho_pi=rho_pi)

def python_f_OL(X, aZ): return g_OL(X)*exp(-(nu_ai(X) + aZ**4)**(1.0/2)) / (nu_ai(X) + aZ**4)**(1.0/2)
f_OL = Expression('g_OL*exp(-sqrt(nu_ai + pow(x[1],4))) / sqrt(nu_ai + pow(x[1],4))', degree=1, g_OL=g_OL, nu_ai=nu_ai)
# Taylor expanded f_OL's
def python_f_OL_taylor(X, aZ): return exp(-(nu_ai(X))**(1.0/2.0))*python_g_OL(X) / (nu_ai(X))**(1.0/2.0) - (exp(-(nu_ai(X))**(1.0/2.0))*python_g_OL(X) * (1 + (nu_ai(X))**(1.0/2.0)))*aZ**4 / (2*(nu_ai(X))**(3.0/2.0))
f_OL_taylor = Expression('exp(-sqrt(nu_ai))*g_OL / sqrt(nu_ai) - (exp(-sqrt(nu_ai))*g_OL * (1 + sqrt(nu_ai)))*pow(x[1],4) / (2*pow(nu_ai, 3.0/2.0))', degree=1, nu_ai=nu_ai, g_OL=g_OL)


import numpy

## ----------------- PRINT NON-Z Dependent Terms ----------
#table_head_array = ['x', 'n', 'Temp', 'g_n_an', 'g_n_cx', 'g_T_an', 'g_T_cx', 'g_Z_an', 'g_Z_cx']
#
#for i in range(len(table_head_array)):
#	print "#" + str(i+1) + "\t" + str(table_head_array[i])
#
#for j in numpy.arange(0.0, 1.0, 0.001):		# x loop
#	print j, n(j), Temp(j), g_n_an(j), g_n_cx(j), g_T_an(j), g_T_cx(j), g_Z_an(j), g_Z_cx(j)

## ----------------- PRINT Z-Dependent Terms --------------
#table_head_array_Z = ['x', 'Z', 'n', 'Temp', 'g_n_bulk', 'g_T_bulk', 'g_Z_bulk', 'f_OL']
#
#for i in range(len(table_head_array_Z)):
#	print "#" + str(i+1) + "\t" + str(table_head_array_Z[i])
#
#Z_values = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
#
#for k in Z_values:
#	for j in numpy.arange(0.0, 1.0, 0.001):			# x loop
#		print j, k, n(j), Temp(j), g_n_bulk(j,k), g_T_bulk(j,k), g_Z_bulk(j,k), f_OL(j,k)

## ----------------- Python Comparison --------------------
#table_head_array_Z = ['Z', 'Python g_n_bulk', 'g_n_bulk', 'Python g_n_bulk_taylor', 'g_n_bulk_taylor', 'Python g_T_bulk', 'g_T_bulk', 'Python g_T_bulk_taylor', 'g_T_bulk_taylor', 'Python g_Z_bulk', 'g_Z_bulk', 'Python g_Z_bulk_taylor', 'g_Z_bulk_taylor', 'Python f_OL', 'f_OL', 'Python f_OL_taylor', 'f_OL_taylor']
#
#print '\n'
#for i in range(len(table_head_array_Z)):
#	print "#" + str(i+1) + "\t" + str(table_head_array_Z[i])
#
#for k in numpy.arange(-3.0, 3.0, 0.01):
#		print k, python_g_n_bulk(j,k), g_n_bulk(j,k), python_g_n_bulk_taylor(j,k), g_n_bulk_taylor(j,k), python_g_T_bulk(j,k), g_T_bulk(j,k), python_g_T_bulk_taylor(j,k), g_T_bulk_taylor(j,k), python_g_Z_bulk(j,k), g_Z_bulk(j,k), python_g_Z_bulk_taylor(j,k), g_Z_bulk_taylor(j,k), python_f_OL(j,k), f_OL(j,k), python_f_OL_taylor(j,k), f_OL_taylor(j,k)

