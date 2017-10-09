import numpy
var('T,Z')

# Global parameters and physical constants
a_cx = 1.0
a_an = 1.0
a_in0 = 0.1
e = 1.602e-19
k_B = 8.617e-5 # Boltzmann in eV
m_e = 9.109e-31
m_i = 1.673e-27
c = 3.0e8
epsilon_0 = 8.854187817e-12
mu_0 = 4*pi*1.0e-7
R_cx = 1.0e-14

# Plasma variables n_e, T_i, n_0
n = 1.0e19
np = 1.0e19
Tp = 1.0e3
n_0(T) =  4.0e17 * a_in0 * ( T / 100 )^(3.0/4)

# JET specifications
a_v = 2.10
a_h = 1.25
a_m = sqrt( (a_v^2 + a_h^2) / 2.0 )
R = 1.3
aspect = a_m / R
I_p = 3.2e6
B_t = 3.45
B_p = mu_0 * I_p / ( 2*pi*a_m )
q = aspect * B_t/B_p
B = sqrt( B_t^2 + B_p^2 )

# Length scales
L_n = 0.02
L_T = 0.015
W_ped = 0.04
L_nn = L_n / W_ped
L_Tn = L_T / W_ped

# Parameters
e_p = 1 + (m_i*n + m_e*n) / (epsilon_0 * B^2)
v_Ti(T) =  (2 * e * T / m_i)^(1.0/2)
v_Te(T) =  (2 * e * T / m_e)^(1.0/2)
rho_pi(T) =  m_i * v_Ti(T) / (e * B_p)
rho_pe(T) =  m_e * v_Te(T) / (e * B_p)
omega_bi(T) =  aspect^(3.0/2) * v_Ti(T) / (q * R)
omega_be(T) =  aspect^(3.0/2) * v_Te(T) / (q * R)

# Frequencies
nu_ei(T) = 1.33e5 * (n * 1.0e-20) / (T * 1.0e-3)^(3.0/2)
nu_ii(T) =  1.2 * (m_e / m_i)^(1.0/2) * nu_ei(T)
nu_in0(T) =  a_in0 * omega_bi(T)
nu_eff(T) =  nu_ii(T) + nu_in0(T)
nu_ai(T) =  nu_ii(T) / omega_bi(T)	# nu_*i
nu_ae(T) =  nu_ei(T) / omega_be(T)	# nu_*e (not used as of now)


# Consolidated constant to reduce clutter
C = nu_ai * aspect^(3.0/2) * nu_ei / nu_ii

print C(1e2)

