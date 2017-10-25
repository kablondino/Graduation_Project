% Coordinate
syms x

%% Global parameters
% Constant degrees of freedom
a_cx        = 1;                                            % constant ~ 1
a_an        = 1;                                            % constant ~ 1
a_in0       = .01;                                          % constant ~ 1

% Plasma variables: n{e}, T{i}, n{0}
%n = piecewise(x < 0, 1.0e19, 0 < x < 1, 1.0e19 * (1 - x), x > 1, 0)
n = symfun( 1.0e19 * x, x )
%Temp = piecewise(x < 0, 2000.0, 0 < x < 1, 2000.0 * (1 - x), x > 1, 0)
Temp = symfun( 2000.0 * x, x )
n_0  = symfun( 4e17*a_in0*(Temp(x) / 100)^(3/4), x );  % [1/m^3]

% Physical constants
charge      = 1.602e-19;                                    % [J/eV]
k           = 8.617e-5;                                     % [eV/K]
m_e         = 9.109e-31;                                    % [kg]
m_i         = 1.673e-27;                                    % [kg]
epsilon_0   = 8.854187817e-12;                              % [F/m]
mu_0        = 4*pi*1e-7;                                    % [N/A^2]
R_cx        = 1e-14;                                        % [m^3/s]

% ASDEX-U tokamak
a_v         = 0.8;                                          % [m]
a_h         = 0.5;                                          % [m]
a_m         = sqrt((a_v^2 + a_h^2)/2);                      % [m]
R           = 1.6;                                          % [m]
aspect      = a_m / R;                                      % [-]
I_p         = 2.0e6                                         % [A]
B_t         = 3.45;                                         % [T]
B_p         = mu_0*I_p / (2*pi*a_m);                        % [T]
q           = aspect*B_t / B_p;                             % [-]
B           = sqrt(B_t^2 + B_p^2);                          % [T]

% Length scales
L_n         = 0.020;                                        % [m]
L_T         = 0.015;                                        % [m]
W_ped       = 0.040;                                        % [m]
L_nn        = L_n / W_ped;                                  % [-]
L_Tn        = L_T / W_ped;                                  % [-]

% Parameters
e_p         = 1 + (m_i*n + m_e*n) / (epsilon_0*B^2);        % [-]
v_Ti        = symfun( sqrt(2*charge*Temp(x) / m_i), x );    % [m/s]
v_Te        = symfun( sqrt(2*charge*Temp(x) / m_e), x );    % [m/s]
rho_pi      = symfun( m_i*v_Ti(x) / (charge*B_p), x );      % [m]
rho_pe      = symfun( m_e*v_Te(x) / (charge*B_p), x );      % [m]
omega_bi    = symfun( (aspect)^(3/2)*v_Ti(x) / (q*R), x );  % [1/s]
omega_be    = symfun( (aspect)^(3/2)*v_Te(x) / (q*R), x );  % [1/s]

% Frequencies
nu_ei       = symfun( 1.33e5*(n(x)*1e-20)/(Temp(x)*1e-3)^(3/2), x );% [1/s]
nu_ii       = symfun( 1.2*sqrt(m_e/m_i)*nu_ei(x), x );      % [1/s]
nu_in0      = symfun( a_in0*omega_bi(x), x );               % [1/s]
nu_eff      = symfun( nu_ii(x) + nu_in0(x), x );            % [1/s]
nu_ai       = symfun( nu_ii(x) / omega_bi(x), x );          % [-]
nu_ae       = symfun( nu_ei(x) / omega_be(x), x );          % [-]

% Consolidated constant to reduce clutter in bulk definitions
C = symfun( nu_ai(x) * aspect^(3.0/2.0) * nu_ei(x) / nu_ii(x), x )

%% Anomalous diffusion coefficient D
D_an = symfun( aspect^2*sqrt(pi) / (2*a_m)*(rho_pe(x) * Temp(x)) / B, x )

syms Z
%% Xi integral results (both Taylor and non-Taylor expanded versions
xi_p = symfun( (4*C(x)*(27*(C(x)^2 + Z^2)^2 - 7*(C(x)^2 - 3*Z^2)*sqrt(nu_ai(x)))*(nu_ai(x))^(7.0/4)) / (189*pi*(C(x)^2 + Z^2)^3), [x Z] )

xi_p_taylor = symfun( (4*(27*C^2*(nu_ai(x))^(7.0/4.0))) / (189*pi*C^3) - (4*(9*C^2*(nu_ai(x))^(7.0/4.0) - 14*(nu_ai(x)^(9.0/4.0)))*Z^2) / (63*pi*C^5), [x Z] )

xi_t = symfun( (2*C(x)*(135*(C(x)^2 + Z^2)^2 - 7*(21*C(x)^4 + 3*Z^2*(-5 + 7*Z^2) + C(x)^2*(5 + 42*Z^2))*sqrt(nu_ai(x)))*(nu_ai(x))^(7.0/4)) / (189*pi*(C(x)^2 + Z^2)^3), [x Z] )

xi_t_taylor = symfun( -(2*(-135*C(x)^2*(nu_ai(x))^(7.0/4.0) + 35*(nu_ai(x))^(9.0/4.0) + 147*C^2*(nu_ai(x))^(9.0/4.0))) / (189*pi*C^3) + (2*(-45*C^2*(nu_ai(x))^(7.0/4.0) + 70*(nu_ai(x))^(9.0/4.0) + 49*C^2*(nu_ai(x))^(9.0/4.0))*Z^2) / (63*pi*C^5), [x Z] )


%% ---------------- g_n coefficients ------------------------------

g_n_an = symfun( -charge* n(x) * D_an(x), x )

g_n_cx = symfun( -m_i*n_0(x)*R_cx * n(x)*Temp(x) / B_p^2, x )

g_n_bulk = symfun( aspect^2*sqrt(pi) / (8*a_m) * n(x) * m_i * rho_pi(x) * (v_Ti(x))^2 * B_p * xi_p(x, Z), [x Z] )

g_n_sum = symfun( g_n_an(x) + g_n_cx(x) + g_n_bulk(x, Z), [x Z] )

%% ---------------- g_T coefficients ------------------------------

g_T_an = symfun( g_n_cx * a_an, x )

g_T_cx = symfun( a_cx * g_n_cx, x )

g_T_bulk = symfun( aspect^2*sqrt(pi) / (8*a_m) * n(x) * m_i * rho_pi(x) * (v_Ti(x))^2 * (B_p*xi_p(x, Z) - B*xi_t(x, Z)), [x Z] )

g_T_sum = symfun( g_T_an(x) + g_T_cx(x) + g_T_bulk(x, Z), [x Z] )

%% ---------------- g_Z coefficients ------------------------------

g_Z_an = symfun( -charge * n(x) * D_an(x) / rho_pi(x), x )

g_Z_cx = symfun( m_i*n_0(x)*R_cx * n(x)*Temp(x) / (rho_pi(x) * B_p^2), x )

g_Z_bulk = symfun( aspect^2*sqrt(pi) / (4*a_m) * n(x) * m_i * (v_Ti(x))^2 * B_p * xi_p(x, Z), [x Z] )

g_Z_sum = symfun( g_Z_an(x) + g_Z_cx(x) + g_Z_bulk(x, Z), [x Z] )


%% ---------------- Orbit loss terms ------------------------------

g_OL = symfun( charge * n(x) * nu_eff(x) * sqrt(aspect) * rho_pi(x), x )

f_OL = symfun( g_OL(x) * exp(-sqrt(nu_ai(x) + Z^4)) / sqrt(nu_ai(x) + Z^4), [x Z] )

