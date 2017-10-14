from parameters import *

# Xi expressions for ion bulk viscosity
#xi_pol = (108*aspect**3*nu_ai**(3.0/2)*nu_ei**2*nu_ii - 28*nu_ii**3) / (189*pi*aspect**(9.0/2)*nu_ai**(3.0/4)*nu_ei**3) + ((-36*aspect**2*nu_ai**(3.0/2.0)*nu_ei**2*nu_ii**3 + 56*nu_ii**5) / (63*pi*aspect**(15.0/2.0)*nu_ai**(11.0/4.0)*nu_ei**5)) * Z**2
xi_p = lambda X: (4*C(X)*(27*(C(X)**2 + Z**2)**2 - 7*(C(X)**2 - 3*Z**2)*nu_ai(X)**(1.0/2))*nu_ai(X)**(7.0/4)) / (189*pi*(C(X)**2 + Z**2)**3)

#xi_tor = (270*aspect**2*nu_ai**(3.0/2.0)*nu_ei**2*nu_ii - 294*aspect**3*nu_ai**2*nu_ei**2*nu_ii - 70*nu_ii**3) / (189*pi*aspect) + ((-90*aspect**3*nu_ai**(3.0/2.0)*nu_ei**2*nu_ii**3 + 98*aspect**3*nu_ai**2*nu_ei**2*nu_ii**3 + 140*nu_ii**5) / (63*pi*aspect**(15.0/2.0)*nu_ai**(11.0/4.0)*nu_ei**5)) * Z**2
xi_t = lambda X: (2*C(X)*(135*(C(X)**2 + Z**2)**2 - 7*(21*C(X)**4 + 3*Z**2*(-5 + 7*Z**2) + C(X)**2*(5 + 42*Z**2))*nu_ai(X)**(1.0/2))*nu_ai(X)**(7.0/4)) / (189*pi*(C(X)**2 + Z**2)**3)

Z = 1

## Anomalous diffusion coefficient
D_an = lambda X: (aspect**2 * pi**(1.0/2)) / (2*a_m) * (rho_pe(X) * T(X)) / B

## ----------------- g_n coefficients -----------
g_n_bulk = lambda X: aspect**2*(pi)**(1.0/2) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * B_p * xi_p(X)
g_n_an = lambda X: -e*n(X)*D_an(X)
g_n_cx = lambda X: -m_i*n_0(X)*R_cx *n(X)*T(X) / B_p**2

## ----------------- g_T coefficients -----------
g_T_bulk = lambda X: aspect**2*(pi)**(1.0/2) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * (B_p*xi_p(X) - B*xi_t(X))
g_T_cx = lambda X: a_cx*g_n_cx(X)
g_T_an = lambda X: g_n_cx(X)*a_an

## ----------------- g_Z coefficients -----------
g_Z_bulk = lambda X: aspect**2*(pi)**(1.0/2) / (4*a_m) * n(X) * m_i * (v_Ti(X))**2 * B_p*xi_p(X)
g_Z_cx = lambda X: m_i*n_0(X)*R_cx *n(X)*T(X) / (rho_pi(X)*B_p**2)
g_Z_an = lambda X: -e*n(X)*D_an(X) / rho_pi(X)

import numpy

print "#x n T D_an g_n_an g_n_cx g_T_an g_T_cx g_Z_an g_Z_cx"
for j in numpy.arange(0.0, 1.0, 0.001):
	print j, n(j), T(j), D_an(j), g_n_bulk(j), g_n_an(j), g_n_cx(j), g_T_bulk(j), g_T_an(j), g_T_cx(j), g_Z_bulk(j), g_Z_an(j), g_Z_cx(j)

