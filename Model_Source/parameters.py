# Global parameters and physical constants
pi = 3.141592653589793
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
n = lambda x: (1 - x) * 1.0e19
T = lambda x: (1 - x) * 2.0e3
np = 1.0e19
Tp = 1.0e3
n_0 = lambda x: 4.0e17 * a_in0 * ( T(x) / 100 )**(3.0/4)

# ASDEX-U specifications
a_v = 0.8	# Vertical minor radius
a_h = 0.5	# Horizontal minor radius
a_m = ( (a_v**2 + a_h**2) / 2.0 )**(1.0/2)
R = 1.6		# Major radius
aspect = a_m / R
I_p = 2.0e6
B_t = 3.9
B_p = mu_0 * I_p / ( 2*pi*a_m )
q = aspect * B_t/B_p
B = ( B_t**2 + B_p**2 )**(1.0/2)

# Length scales
L_n = 0.02
L_T = 0.015
W_ped = 0.04
L_nn = L_n / W_ped
L_Tn = L_T / W_ped

# Parameters
e_p = lambda X: 1.0 + (m_i*n(X) + m_e*n(X)) / (epsilon_0 * B**2)

v_Ti = lambda X: (2 * e * T(X) / m_i)**(1.0/2)

v_Te = lambda X: (2 * e * T(X) / m_e)**(1.0/2)

rho_pi = lambda X: m_i * v_Ti(X) / (e * B_p)

rho_pe = lambda X: m_e * v_Te(X) / (e * B_p)

omega_bi = lambda X: aspect**(3.0/2) * v_Ti(X) / (q * R)

omega_be = lambda X: aspect**(3.0/2) * v_Te(X) / (q * R)

# Frequencies
nu_ei = lambda X: 1.33e5 * (n(X) * 1.0e-20) / (T(X) * 1.0e-3)**(3.0/2)

nu_ii = lambda X: 1.2 * (m_e / m_i)**(1.0/2) * nu_ei(X)

nu_in0 = lambda X: a_in0 * omega_bi(X)

nu_eff = lambda X: nu_ii(X) + nu_in0(X)

nu_ai = lambda X: nu_ii(X) / omega_bi(X)	# nu_*i

nu_ae = lambda X: nu_ei(X) / omega_be(X)	# nu_*e (not used as of now)


# Consolidated constant to reduce clutter
C = lambda X: nu_ai(X) * aspect**(3.0/2) * nu_ei(X) / nu_ii(X)

#print C(0.5)

