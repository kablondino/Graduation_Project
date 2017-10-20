from parameters import *

# Xi expressions for ion bulk viscosity
#def xi_p(X, aZ): return (4*C(X)*(27*(C(X)**2 + aZ**2)**2 - 7*(C(X)**2 - 3*aZ**2)*nu_ai(X)**(1.0/2))*nu_ai(X)**(7.0/4)) / (189*pi*(C(X)**2 + aZ**2)**3)
#xi_p = Expression('(4*C*(27*pow((pow(C,2) + pow(Z,2)),2) - 7*(pow(C,2) - 3*pow(Z,2))*sqrt(nu_ai))*pow(nu_ai, 7.0/4.0)) / (189*pi*pow((pow(C,2) + pow(Z,2)),3))', degree=1, C=C, Z=Z, nu_ai=nu_ai)
## Taylor expanded version
#xi_p_taylor = Expression('(4*(27*pow(C,2)*pow(nu_ai,7.0/4.0))) / (189*pi*pow(C,3)) - (4*(9*pow(C,2)*pow(nu_ai, 7.0/4.0) - 14*pow(nu_ai, 9.0/4.0))*pow(Z,2)) / (63*pi*pow(C,5))', degree=1, C=C, Z=Z, nu_ai=nu_ai)

#def xi_t(X, aZ): return (2*C(X)*(135*(C(X)**2 + aZ**2)**2 - 7*(21*C(X)**4 + 3*aZ**2*(-5 + 7*aZ**2) + C(X)**2*(5 + 42*aZ**2))*nu_ai(X)**(1.0/2))*nu_ai(X)**(7.0/4)) / (189*pi*(C(X)**2 + aZ**2)**3)
#xi_t = Expression('(2*C*(135*pow((pow(C,2) + pow(Z,2)),2) - 7*(21*pow(C,4) + 3*pow(Z,2)*(-5 + 7*pow(Z,2)) + pow(C,2)*(5 + 42*pow(Z,2)))*sqrt(nu_ai))*pow(nu_ai, 7.0/4.0)) / (189*pi*pow((pow(C,2) + pow(Z,2)),3))', degree=1, C=C, Z=Z, nu_ai=nu_ai)
## Taylor Expanded version
#xi_t_taylor = Expression('-(2*(-135*pow(C,2)*pow(nu_ai, 7.0/4.0) + 35*pow(nu_ai, 9.0/4.0) + 147*pow(C,2)*pow(nu_ai, 9.0/4.0))) / (189*pi*pow(C,3)) + (2*(-45*pow(C,2)*pow(nu_ai, 7.0/4.0) + 70*pow(nu_ai, 9.0/4.0) + 49*pow(C,2)*pow(nu_ai, 9.0/4.0))*pow(Z,2)) / (63*pi*pow(C,5))', degree=1, C=C, Z=Z, nu_ai=nu_ai)

## Anomalous diffusion coefficient
#def D_an(X): return (aspect**2 * pi**(1.0/2)) / (2*a_m) * (rho_pe(X) * Temp(X)) / B
D_an = Expression('pow(aspect,2)*sqrt(pi) / (2*a_m) * (rho_pe * Temp) / B', degree=1, aspect=aspect, a_m=a_m, rho_pe=rho_pe, Temp=Temp, B=B)

## ----------------- g_n coefficients -----------
#def g_n_bulk(X, aZ): return aspect**2*(pi)**(1.0/2) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * B_p * xi_p(X, aZ)
#g_n_bulk = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * B_p * xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p)

#def g_n_an(X): return -charge*n(X)*D_an(X)
g_n_an = Expression('-charge*n*D_an', degree=1, charge=charge, n=n, D_an=D_an)

#def g_n_cx(X): return -m_i*n_0(X)*R_cx *n(X)*Temp(X) / B_p**2
g_n_cx = Expression('-m_i*n_0*R_cx *n*Temp / pow(B_p,2)', degree=1, m_i=m_i, n_0=n_0, R_cx=R_cx, n=n, Temp=Temp, B_p=B_p)


## ----------------- g_T coefficients -----------
#def g_T_bulk(X, aZ): return aspect**2*(pi)**(1.0/2) / (8*a_m) * n(X) * m_i * rho_pi(X) * (v_Ti(X))**2 * (B_p*xi_p(X, aZ) - B*xi_t(X, aZ))
#g_T_bulk = Expression('pow(aspect,2)*sqrt(pi) / (8*a_m) * n * m_i * rho_pi * pow(v_Ti,2) * (B_p*xi_p - B*xi_t)', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, rho_pi=rho_pi, v_Ti=v_Ti, B_p=B_p, B=B, xi_p=xi_p, xi_t=xi_t)

#def g_T_an(X): return g_n_cx(X)*a_an
g_T_an = Expression('g_n_an*a_an', degree=1, g_n_an=g_n_an, a_an=a_an)

#def g_T_cx(X): return a_cx*g_n_cx(X)
g_T_cx = Expression('a_cx*g_n_cx', degree=1, a_cx=a_cx, g_n_cx=g_n_cx)


## ----------------- g_Z coefficients -----------
#def g_Z_bulk(X, aZ): return aspect**2*(pi)**(1.0/2) / (4*a_m) * n(X) * m_i * (v_Ti(X))**2 * B_p*xi_p(X, aZ)
#g_Z_bulk = Expression('pow(aspect,2)*sqrt(pi) / (4*a_m) * n * m_i * pow(v_Ti,2) * B_p*xi_p', degree=1, aspect=aspect, a_m=a_m, n=n, m_i=m_i, v_Ti=v_Ti, B_p=B_p, xi_p=xi_p)

#def g_Z_an(X): return -charge*n(X)*D_an(X) / rho_pi(X)
g_Z_an = Expression('-charge*n*D_an / rho_pi', degree=1, charge=charge, n=n, D_an=D_an, rho_pi=rho_pi)

#def g_Z_cx(X): return m_i*n_0(X)*R_cx *n(X)*Temp(X) / (rho_pi(X)*B_p**2)
g_Z_cx = Expression('m_i*n_0*R_cx * n*Temp / (rho_pi*pow(B_p,2))', degree=1, m_i=m_i, n_0=n_0, R_cx=R_cx, n=n, Temp=Temp, rho_pi=rho_pi, B_p=B_p)

## ----------------- g_OL coefficient -----------
#def g_OL(X): return charge*n(X)*nu_eff(X)*(aspect)**(1.0/2)*rho_pi(X)
g_OL = Expression('charge*n*nu_eff*sqrt(aspect)*rho_pi', degree=1, charge=charge, n=n, nu_eff=nu_eff, aspect=aspect, rho_pi=rho_pi)

#def f_OL(X, aZ): return g_OL(X)*exp(-(nu_ai(X) + aZ**4)**(1.0/2)) / (nu_ai(X) + aZ**4)**(1.0/2)
#f_OL = Expression('g_OL*exp(-sqrt(nu_ai + pow(Z,4))) / sqrt(nu_ai + pow(Z,4))', degree=1, g_OL=g_OL, nu_ai=nu_ai, Z=Z)
# Taylor expanded f_OL
#f_OL_taylor = Expression('exp(-sqrt(nu_ai))*g_OL / sqrt(nu_ai) - (exp(-sqrt(nu_ai))*g_OL * (1 + sqrt(nu_ai)))*pow(Z,4) / (2*pow(nu_ai, 3.0/2.0))', degree=1, nu_ai=nu_ai, Z=Z, g_OL=g_OL)


#import numpy
#
#print "#x n Temp Z D_an g_n_an g_n_cx g_T_an g_T_cx g_Z_an g_Z_cx g_OL f_OL g_n_bulk g_T_bulk g_Z_bulk"
#for j in numpy.arange(0.0, 1.0, 0.005):
#	for k in numpy.arange(-3.0, 3.0, 0.01):
#		print j, n(j), Temp(j), k, D_an(j), g_n_an(j), g_n_cx(j), g_T_an(j), g_T_cx(j), g_Z_an(j), g_Z_cx(j), g_OL(j), f_OL(j,k), g_n_bulk(j,k), g_T_bulk(j,k), g_Z_bulk(j,k)

