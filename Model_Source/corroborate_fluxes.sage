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


## ------------ Ion bulk (parallel) viscosity -------------------
# Xi functions
#xi_p1 = (108*aspect^3*nu_ai^(3.0/2)*nu_ei^2*nu_ii - 28*nu_ii^3) / (189*pi*aspect^(9.0/2)*nu_ai^(3.0/4)*nu_ei^3) + (-36*aspect^2*nu_ai^(3.0/2)*nu_ei^2*nu_ii^3 + 56*nu_ii^5) / (63*pi*aspect^(15.0/2)*nu_ai^(11.0/4)*nu_ei^5) * Z^2
xi_p = (4*(27*C^2*nu_ai^(7/4) - 7*nu_ai^(9/4))) / (189*pi*C^3) + (-4*(9*C^2*nu_ai^(7/4) - 14*nu_ai^(9/4))) / (63*pi*C^5)*Z^2

# NUMERICAL INTEGRAL
#var('x,y')
#Fp_BV = (1/pi)*        x^2*exp(-x)*(C/sqrt(x)) / ((y + Z/sqrt(x))^2 + (C/sqrt(x))^2)
#Ft_BV = (1/pi)*(5/2-x)*x^2*exp(-x)*(C/sqrt(x)) / ((y + Z/sqrt(x))^2 + (C/sqrt(x))^2)
#I_p = numerical_integral(lambda y: Fp_BV, -1, 1)

#xi_t = (270*aspect^2*nu_ai^(3.0/2)*nu_ei^2*nu_ii - 294*aspect^3*nu_ai^2*nu_ei^2*nu_ii - 70*nu_ii^3) / (189*pi*aspect) + (-90*aspect^3*nu_ai^(3.0/2)*nu_ei^2*nu_ii^3 + 98*aspect^3*nu_ai^2*nu_ei^2*nu_ii^3 + 140*nu_ii^5) / (63*pi*aspect^(15.0/2)*nu_ai^(11.0/4)*nu_ei^5) * Z^2
xi_t = (2*(-135*C^2*nu_ai^(7/4) + 35*nu_ai^(9/4) * 147*C^2*nu_ai^(9/4))) / (189*pi*C^3) + (2*(-45*C^2*nu_ai^(7/4) + 70*nu_ai^(9/4) + 49*C^2*nu_ai(9/4))) / (63*pi*C^5) * Z^2

#g_n_bulk(T,Z) = (aspect^2*sqrt(pi)) / (8*a_m) * n*m_i*rho_pi*v_Ti^2 * B_p * xi_p
#g_T_bulk(T,Z) = (aspect^2*sqrt(pi)) / (8*a_m) * n*m_i*rho_pi*v_Ti^2 * (B_p*xi_p - B*xi_t)
#g_Z_bulk(T,Z) = (aspect^2*sqrt(pi)) / (4*a_m) * n*m_i*v_Ti^2 * B_p * xi_p

g_pol_bulk = (aspect^2*sqrt(pi)) / (4*a_m) * n*m_i*v_Ti * B * xi_p
g_tor_bulk = (aspect^2*sqrt(pi)) / (4*a_m) * n*m_i*v_Ti * B * xi_t

up0(T,Z) = -rho_pi*v_Ti / (2*L_T)
up = (v_Ti*B_p / B) * (Z + rho_pi/2 * (1/L_n + 1/L_T))

G_bulk(T,Z) = g_tor_bulk*up + g_pol_bulk*up0

#plot3d(G_bulk, (T,1e2,2e3), (Z,-3,3), adaptive=True).show()


## ------------ Anomalous loss calculation ----------------------
# Sage definition
D_an(T) = (aspect^2 * sqrt(pi)) / (2*a_m) * (rho_pe * T) / B
g_an_n(T) = -e*n*D_an
g_an_T(T) = -e*n*D_an*a_an
g_an_Z(T) = -e*n*D_an / rho_pi
G_an(T, Z) = g_an_n/L_n + g_an_T/L_T + g_an_Z*Z


#plot3d(G_an, (T,1e2,2e3), (Z,-3,3), adaptive=True).show()


## ------------ Charge exchange friction ------------------------
g_cx_n(T) = -m_i*n_0*R_cx *n*T / B_p^2
g_cx_T(T) = a_cx*g_cx_n
g_cx_Z(T) = m_i*n_0*R_cx *n*T / (rho_pi*B_p^2)

G_cx(T,Z) = g_cx_n/L_n + g_cx_T/L_T + g_cx_Z*Z

#plot3d(G_cx, (T,1e2,2e3), (Z,-3,3), adaptive=True).show()


## ------------ Ion Orbit Loss ----------------------------------
g_ol = e*n*nu_eff*sqrt(aspect)*rho_pi
G_ol(T,Z) = g_ol*exp(-sqrt(nu_ai + Z^4)) / sqrt(nu_ai + Z^4)
#plot3d(G_ol, (T,1e2,2e3), (Z,-3,3), adaptive=True).show()


## ----------- Spit out numerical data --------------------------
# Open file, if necessary
with open('/home/kabv/Documents/Masters/Graduation_Project/Model_Source/fluxes.dat', 'w') as the_file:
	# Print labels
	the_file.write("T Z n_0 rho_pi rho_pe v_Ti v_Te omega_bi omega_be nu_ii nu_ei bulk_visc anom_loss charge_ex orbit_loss\n")

	# Loop over EVERYTHING!
	for i in numpy.arange(1e2, 2e3 + 25, 25):	# i is Temperature
		for j in numpy.arange(-3, 3.1, 0.1):	# j is Z (electric field)
			the_file.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % (i, j, N(n_0(i)), N(rho_pi(i)), N(rho_pe(i)), N(v_Ti(i)), N(v_Te(i)), N(omega_bi(i)), N(omega_be(i)), N(nu_ii(i)), N(nu_ei(i)), N(G_bulk(i,j)), N(G_an(i,j)), N(G_cx(i,j)), N(G_ol(i,j))))

the_file.close()

