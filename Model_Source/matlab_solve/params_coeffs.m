clear all; close all; clc

% Coordinate
syms x

%% Global parameters
% Constant degrees of freedom
global a_cx
a_cx        = 1;                                            % constant ~ 1
global a_an
a_an        = 1;                                            % constant ~ 1
global a_in0
a_in0       = .01;                                          % constant ~ 1

% Plasma variables: n{e}, T{i}, n{0}
%n = piecewise(x < 0, 1.0e19, 0 < x < 1, 1.0e19 * (1 - x), x > 1, 0)
global n
n = symfun( 1.0e19 * (1 - x), x )
%Temp = piecewise(x < 0, 2000.0, 0 < x < 1, 2000.0 * (1 - x), x > 1, 0)
global Temp
Temp = symfun( 2000.0 * (1 - x), x )
global n_0
n_0  = symfun( 4e17*a_in0*(Temp(x)/100)^(3/4), x );  % [1/m^3]

% Physical constants
global charge
charge      = 1.602e-19;                                    % [J/eV]
global k
k           = 8.617e-5;                                     % [eV/K]
global m_e
m_e         = 9.109e-31;                                    % [kg]
global m_i
m_i         = 1.673e-27;                                    % [kg]
global c
c           = 3e8;                                          % [m/s]
global epsilon_0
epsilon_0   = 8.854187817e-12;                              % [F/m]
global mu_0
mu_0        = 4*pi*1e-7;                                    % [N/A^2]
global R_cx
R_cx        = 1e-14;                                        % [m^3/s]

% ASDEX-U tokamak
global a_v
a_v         = 0.8;                                          % [m]
global a_h
a_h         = 0.5;                                          % [m]
global a_m
a_m         = sqrt((a_v^2 + a_h^2)/2);                      % [m]
global R
R           = 1.6;                                          % [m]
global aspect
aspect      = a_m / R;                                      % [-]
global I_p
I_p         = 2.0e6                                         % [A]
global B_t
B_t         = 3.45;                                         % [T]
global B_p
B_p         = mu_0*I_p / (2*pi*a_m);                        % [T]
global q
q           = aspect*B_t / B_p;                             % [-]
global B
B           = sqrt(B_t^2 + B_p^2);                          % [T]

% Length scales
global L_n
L_n         = 0.020;                                        % [m]
global L_T
L_T         = 0.015;                                        % [m]
global W_ped
W_ped       = 0.040;                                        % [m]
global L_nn
L_nn        = L_n / W_ped;                                  % [-]
global L_Tn
L_Tn        = L_T / W_ped;                                  % [-]

% Parameters
global e_p
e_p         = 1 + (m_i*n + m_e*n) / (epsilon_0*B^2);        % [-]
global v_Ti
v_Ti        = symfun( sqrt(2*charge*Temp(x) / m_i), x );    % [m/s]
global v_Te
v_Te        = symfun( sqrt(2*charge*Temp(x) / m_e), x );    % [m/s]
global rho_pi
rho_pi      = symfun( m_i*v_Ti(x) / (charge*B_p), x );      % [m]
global rho_pe
rho_pe      = symfun( m_e*v_Te(x) / (charge*B_p), x );      % [m]
global omega_bi
omega_bi    = symfun( (aspect)^(3/2)*v_Ti(x) / (q*R), x );  % [1/s]
global omega_be
omega_be    = symfun( (aspect)^(3/2)*v_Te(x) / (q*R), x );  % [1/s]

% Frequencies
global nu_ei
nu_ei       = symfun( 1.33e5*(n(x)*1e-20)/(Temp(x)*1e-3)^(3/2), x );% [1/s]
global nu_ii
nu_ii       = symfun( 1.2*sqrt(m_e/m_i)*nu_ei(x), x );      % [1/s]
global nu_in0
nu_in0      = symfun( a_in0*omega_bi(x), x );               % [1/s]
global nu_eff
nu_eff      = symfun( nu_ii(x) + nu_in0(x), x );            % [1/s]
global nu_ai
nu_ai       = symfun( nu_ii(x) / omega_bi(x), x );          % [-]
global nu_ae
nu_ae       = symfun( nu_ei(x) / omega_be(x), x );          % [-]

% Consolidated constant to reduce clutter in bulk definitions
global C
C = symfun( nu_ai(x) * aspect^(3.0/2.0) * nu_ei(x) / nu_ii(x), x )

%% Anomalous diffusion coefficient D
global D_an
D_an = symfun( aspect^2*sqrt(pi) / (2*a_m)*(rho_pe(x) * Temp(x)) / B, x )

syms Z
%% Xi integral results (both Taylor and non-Taylor expanded versions
global xi_p
xi_p = symfun( (4*C(x)*(27*(C(x)^2 + Z^2)^2 - 7*(C(x)^2 - 3*Z^2)*sqrt(nu_ai(x)))*(nu_ai(x))^(7.0/4)) / (189*pi*(C(x)^2 + Z^2)^3), [x Z] )

global xi_p_taylor
xi_p_taylor = symfun( (4*(27*C^2*(nu_ai(x))^(7.0/4.0))) / (189*pi*C^3) - (4*(9*C^2*(nu_ai(x))^(7.0/4.0) - 14*(nu_ai(x)^(9.0/4.0)))*Z^2) / (63*pi*C^5), [x Z] )

global xi_t
xi_t = symfun( (2*C(x)*(135*(C(x)^2 + Z^2)^2 - 7*(21*C(x)^4 + 3*Z^2*(-5 + 7*Z^2) + C(x)^2*(5 + 42*Z^2))*sqrt(nu_ai(x)))*(nu_ai(x))^(7.0/4)) / (189*pi*(C(x)^2 + Z^2)^3), [x Z] )

global xi_t_taylor
xi_t_taylor = symfun( -(2*(-135*C(x)^2*(nu_ai(x))^(7.0/4.0) + 35*(nu_ai(x))^(9.0/4.0) + 147*C^2*(nu_ai(x))^(9.0/4.0))) / (189*pi*C^3) + (2*(-45*C^2*(nu_ai(x))^(7.0/4.0) + 70*(nu_ai(x))^(9.0/4.0) + 49*C^2*(nu_ai(x))^(9.0/4.0))*Z^2) / (63*pi*C^5), [x Z] )


%% ---------------- g_n coefficients ------------------------------

global g_n_an
g_n_an = symfun( -charge* n(x) * D_an(x), x )

global g_n_cx
g_n_cx = symfun( -m_i*n_0(x)*R_cx * n(x)*Temp(x) / B_p^2, x )

global g_n_bulk
g_n_bulk = symfun( aspect^2*sqrt(pi) / (8*a_m) * n(x) * m_i * rho_pi(x) * (v_Ti(x))^2 * B_p * xi_p(x, Z), [x Z] )

global g_n_sum
g_n_sum = symfun( g_n_an(x) + g_n_cx(x) + g_n_bulk(x, Z), [x Z] )

%% ---------------- g_T coefficients ------------------------------

global g_T_an
g_T_an = symfun( g_n_cx * a_an, x )

global g_T_cx
g_T_cx = symfun( a_cx * g_n_cx, x )

global g_T_bulk
g_T_bulk = symfun( aspect^2*sqrt(pi) / (8*a_m) * n(x) * m_i * rho_pi(x) * (v_Ti(x))^2 * (B_p*xi_p(x, Z) - B*xi_t(x, Z)), [x Z] )

global g_T_sum
g_T_sum = symfun( g_T_an(x) + g_T_cx(x) + g_T_bulk(x, Z), [x Z] )

%% ---------------- g_Z coefficients ------------------------------

global g_Z_an
g_Z_an = symfun( -charge * n(x) * D_an(x) / rho_pi(x), x )

global g_Z_cx
g_Z_cx = symfun( m_i*n_0(x)*R_cx * n(x)*Temp(x) / (rho_pi(x) * B_p^2), x )

global g_Z_bulk
g_Z_bulk = symfun( aspect^2*sqrt(pi) / (4*a_m) * n(x) * m_i * (v_Ti(x))^2 * B_p * xi_p(x, Z), [x Z] )

global g_Z_sum
g_Z_sum = symfun( g_Z_an(x) + g_Z_cx(x) + g_Z_bulk(x, Z), [x Z] )


%% ---------------- Orbit loss terms ------------------------------

global g_OL
g_OL = symfun( charge * n(x) * nu_eff(x) * sqrt(aspect) * rho_pi(x), x )

global f_OL
f_OL = symfun( g_OL(x) * exp(-sqrt(nu_ai(x) + Z^4)) / sqrt(nu_ai(x) + Z^4), [x Z] )

