from dolfin import *	# Required for 'Expression(...)'

# Global parameters and physical constants
pi = 3.141592653589793
a_cx = 1.0
a_an = 1.0
a_in0 = 0.1
charge = 1.602e-19
k_B = 8.617e-5		# Boltzmann in eV
m_e = 9.109e-31
m_i = 1.673e-27
c = 3.0e8
epsilon_0 = 8.854187817e-12
mu_0 = 4*pi*1.0e-7
R_cx = 1.0e-14

## Plasma profiles n_e, T_i, n_0
# Plasma edge is at x = 1
#n = Expression('(x[0] <= 0) ? (1.0e+19) : ((x[0] >= 1) ? (0) :  1.0e+19 * (1 - x[0]))', degree=1)
#
#Temp = Expression('(x[0] <= 0) ? (2000.0) : ((x[0] >= 1) ? (0) : 2000.0 * (1 - x[0]))', degree=1)

# Plasma edge is at x = 0
n = Expression('((x[0] <= 0) ? (1.0e+19) : ((x[0] >= 1) ? (0) : (-1.0e+19*x[0] + 1.0e+19)))', degree=1)

Temp = Expression('((x[0] <= 0) ? (2000.0) : ((x[0] >= 1) ? (0) : (-2000.0*x[0] + 2000.0)))', degree=1)

#def n_0(x): return 4.0e17 * a_in0 * ( Temp(x) / 100 )**(3.0/4.0)
n_0 = Expression('4.0e17 * a_in0 * pow((Temp / 100.0), 3.0/4.0)', degree=1, a_in0=a_in0, Temp=Temp)

# ASDEX-U specifications
a_v = 0.8	# Vertical minor radius
a_h = 0.5	# Horizontal minor radius
a_m = ( (a_v**2 + a_h**2) / 2.0 )**(1.0/2.0)
R = 1.6		# Major radius
aspect = a_m / R
I_p = 2.0e6
B_t = 3.9
B_p = mu_0 * I_p / ( 2*pi*a_m )
q = aspect * B_t/B_p
B = ( B_t**2 + B_p**2 )**(1.0/2.0)

# Length scales
L_n = 0.02
L_T = 0.015
W_ped = 0.04
L_nn = L_n / W_ped
L_Tn = L_T / W_ped

# Parameters
#def e_p(X): return 1.0 + (m_i*n(X) + m_e*n(X)) / (epsilon_0 * B**2)
e_p = Expression('1.0 + (m_i*n / (epsilon_0 * pow(B,2)))', degree=1, m_i=m_i, epsilon_0=epsilon_0, B=B, n=n)

#def v_Ti(X): return (2.0 * charge * Temp(X) / m_i)**(1.0/2.0)
v_Ti = Expression('sqrt(2.0 * charge * Temp / m_i)', degree=1, charge=charge, m_i=m_i, Temp=Temp)

#def v_Te(X): return (2.0 * charge * Temp(X) / m_e)**(1.0/2.0)
v_Te = Expression('sqrt(2.0 * charge * Temp / m_e)', degree=1, charge=charge, Temp=Temp, m_e=m_e)

#def rho_pi(X): return m_i * v_Ti(X) / (charge * B_p)
rho_pi = Expression('m_i * v_Ti / (charge * B_p)', degree=1, m_i=m_i, v_Ti=v_Ti, charge=charge, B_p=B_p)

#def rho_pe(X): return m_e * v_Te(X) / (charge * B_p)
rho_pe = Expression('m_e * v_Te / (charge * B_p)', degree=1, m_e=m_e, v_Te=v_Te, charge=charge, B_p=B_p)

#def omega_bi(X): return aspect**(3.0/2.0) * v_Ti(X) / (q * R)
omega_bi = Expression('pow(aspect, 3.0/2.0) * v_Ti / (q * R)', degree=1, aspect=aspect, v_Ti=v_Ti, q=q, R=R)

#def omega_be(X): return aspect**(3.0/2.0) * v_Te(X) / (q * R)
omega_be = Expression('pow(aspect, 3.0/2.0) * v_Te / (q * R)', degree=1, aspect=aspect, v_Te=v_Te, q=q, R=R)

## Collision Frequencies
#def nu_ei(X): return 1.33e5 * (n(X) * 1.0e-20) / (Temp(X) * 1.0e-3)**(3.0/2.0)
nu_ei = Expression('1.33e5 * (n * 1.0e-20) / pow((Temp * 1.0e-3), 3.0/2.0)', degree=1, n=n, Temp=Temp)

#def nu_ii(X): return 1.2 * (m_e / m_i)**(1.0/2.0) * nu_ei(X)
nu_ii = Expression('1.2 * sqrt(m_e / m_i) * nu_ei', degree=1, m_e=m_e, m_i=m_i, nu_ei=nu_ei)

#def nu_in0(X): return a_in0 * omega_bi(X)
nu_in0 = Expression('a_in0 * omega_bi', degree=1, a_in0=a_in0, omega_bi=omega_bi)

#def nu_eff(X): return nu_ii(X) + nu_in0(X)
nu_eff = Expression('nu_ii + nu_in0', degree=1, nu_ii=nu_ii, nu_in0=nu_in0)

#def nu_ai(X): return nu_ii(X) / omega_bi(X)	# nu_*i
nu_ai = Expression('nu_ii / omega_bi', degree=1, nu_ii=nu_ii, omega_bi=omega_bi)

#def nu_ae(X): return nu_ei(X) / omega_be(X)	# nu_*e (not used as of now)
nu_ae = Expression('nu_ei / omega_be', degree=1, nu_ei=nu_ei, omega_be=omega_be)

## Consolidated constant to reduce clutter
#def C(X): return nu_ai(X) * aspect**(3.0/2.0) * nu_ei(X) / nu_ii(X)
C = Expression('nu_ai * pow(aspect, 3.0/2.0) * nu_ei / nu_ii', degree=1, nu_ai=nu_ai, aspect=aspect, nu_ei=nu_ei, nu_ii=nu_ii)

