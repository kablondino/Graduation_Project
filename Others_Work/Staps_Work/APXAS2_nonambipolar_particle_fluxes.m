%% A.2 Nonambipolar Particle Fluxes
% V0: 2016/05/24
% V1: 2016/11/09

clear all; close all; clc

%% Global parameters
% Constant degrees of freedom
a_cx        = 1;                                            % constant ~ 1
a_an        = 1;                                            % constant ~ 1
a_in0       = .01;                                          % constant ~ 1

% Plasma variables: n{e}, T{i}, n{0}
n           = 1e19;                                         % [1/m^3]
np          = 1e19;                                         % [1/m^3]
T           = 1e2:25:2e3;                                   % [eV]
Tp          = 1e3;                                          % [eV]
Z           = -3:.1:3;                                      % [-]
n_0         = @(T) 4e17*a_in0*(T/100)^(3/4);                % [1/m^3]

% Physical constants
e           = 1.602e-19;                                    % [J/eV]
k           = 8.617e-5;                                     % [eV/K]
m_e         = 9.109e-31;                                    % [kg]
m_i         = 1.673e-27;                                    % [kg]
c           = 3e8;                                          % [m/s]
e_0         = 8.854187817e-12;                              % [F/m]
mu_0        = 4*pi*1e-7;                                    % [N/A^2]
R_cx        = 1e-14;                                        % [m^3/s]

% JET tokamak
av          = 2.10;                                         % [m]
ah          = 1.25;                                         % [m]
am          = sqrt((av^2+ah^2)/2);                          % [m]
R           = 1.3;                                          % [m]
aspect      = am/R;                                         % [-]
Ip          = 3.2e6;                                        % [A]
Bt          = 3.45;                                         % [T]
Bp          = mu_0*Ip/(2*pi*am);                            % [T]
q           = aspect*Bt/Bp;                                 % [-]
B           = sqrt(Bt^2+Bp^2);                              % [T]

% Length scales
L_n         = 0.020;                                        % [m]
L_T         = 0.015;                                        % [m]
W_ped       = 0.040;                                        % [m]
L_nn        = L_n/W_ped;                                    % [-]
L_Tn        = L_T/W_ped;                                    % [-]

% Parameters
e_p         = 1+(m_i*n+m_e*n)/(e_0*B^2);                    % [-]
v_Ti        = @(T) (2*e*T/m_i)^(1/2);                       % [m/s]
rho_pi      = @(T) m_i*v_Ti(T)/(e*Bp);                      % [m]
v_Te        = @(T) (2*e*T/m_e)^(1/2);                       % [m/s]
rho_pe      = @(T) m_e*v_Te(T)/(e*Bp);                      % [m]
omega_bi    = @(T) (aspect)^(3/2)*v_Ti(T)/(q*R);            % [1/s]
omega_be    = @(T) (aspect)^(3/2)*v_Te(T)/(q*R);            % [1/s]

% Frequencies
nu_ei       = @(T) 1.33e5*(n*1e-20)/(T*1e-3)^(3/2);         % [1/s]
nu_ii       = @(T) 1.2*(m_e/m_i)^(1/2)*nu_ei(T);            % [1/s]
nu_in0      = @(T) a_in0*omega_bi(T);                       % [1/s]
nu_eff      = @(T) nu_ii(T)+nu_in0(T);                      % [1/s]
nu_efi      = @(T) nu_ii(T)/omega_bi(T);                    % [-]
nu_efe      = @(T) nu_ei(T)/omega_be(T);                    % [-]

%% Bulk viscosity
g_bv = [];
for j = 1:length(T)
    for i = 1:length(Z)
        up(i,j)     = (v_Ti(T(j))*Bp/B)*(Z(i)+rho_pi(T(j))/2*(1/L_n + 1/L_T));
        upm         = Z(i);
        % Energy integrals: Ip, It
        Fp_BV       = @(x,y) 1/pi*(1    ).*x.^2.*exp(-x).*(nu_efi(T(j))*aspect^(3/2)*(nu_ei(T(j))./nu_ii(T(j))*sqrt(x)))./((y+upm./sqrt(x)).^2+(nu_efi(T(j))*aspect^(3/2)*(nu_ei(T(j))./nu_ii(T(j))*sqrt(x))).^2);
        Ip_BV       = integral2(Fp_BV,0,nu_efi(T(j))^(1/2),-1,1);
        Ft_BV       = @(x,y) 1/pi*(5/2-x).*x.^2.*exp(-x).*(nu_efi(T(j))*aspect^(3/2)*(nu_ei(T(j))./nu_ii(T(j))*sqrt(x)))./((y+upm./sqrt(x)).^2+(nu_efi(T(j))*aspect^(3/2)*(nu_ei(T(j))./nu_ii(T(j))*sqrt(x))).^2);
        It_BV       = integral2(Ft_BV,0,nu_efi(T(j))^(1/2),-1,1);
        % Bulk viscosity flux
        g_tbv(i,j)  = (pi)^(1/2)*aspect^2*n*m_i*v_Ti(T(j))*B/(4*am)*(Ip_BV);
        g_pbv(i,j)  = (pi)^(1/2)*aspect^2*n*m_i*v_Ti(T(j))*B/(4*am)*(It_BV);
        up0(i,j)    = -rho_pi(T(j))*v_Ti(T(j))/(2*L_T);
        G_bv(i,j)   = g_tbv(i,j)*up(i,j)+g_pbv(i,j)*up0(i,j);
    end
    fprintf('%0.5g / %0.5g \n',j,length(T))
end

%% Plot mesh plot: (T,Z,G{bv})
[Tm,Zm] = meshgrid(T,Z);
figure; mesh(Tm*1e-3,Zm,-G_bv);
zlim([-1 2]); colorbar;
view([-45 45]); shading interp;

% Format
title('Bulk viscosity current','interpreter','latex')
xlabel('$T$ [keV]','interpreter','latex')
ylabel('$Z$','interpreter','latex')
zlabel('$- e \Gamma^{bv}$ [A/m$^{2}$]','interpreter','latex')

% Save as .pdf file
savepdf('APXAS2_bulk_viscosity_3d')

%% Plot G{bv}(T = 1e3 [eV])
Zp      = min(Z):.01:max(Z);
G_bvp   = interp1(Z,G_bv(:,T==Tp),Zp,'spline');
figure; plot(Zp,-G_bvp)

% Format
title(['Bulk viscosity current (T = ' num2str(Tp/1e3) ' [keV])'],'interpreter','latex')
xlabel('$Z$ [keV]','interpreter','latex')
ylabel('$- e \Gamma^{bv}$ [A/m$^{2}$]','interpreter','latex')

% Save as .pdf file
savepdf('APXAS2_bulk_viscosity_2d')

%% Anomalous loss

for j = 1:length(T)
    for i = 1:length(Z)
        g_an(i,j)   = -( (pi)^(1/2)*aspect^2 )/( 2*am )*( rho_pe(T(j))*e*T(j) )/( B )*n;
        G_an(i,j)   = g_an(i,j)*(1/L_n + a_an/L_T + Z(i)/rho_pi(T(j)));
    end
end

%% Plot mesh plot: (T,Z,G{an})
[Tm,Zm] = meshgrid(T,Z);
figure; mesh(Tm*1e-3,Zm,G_an);
colorbar; view([135 45]); shading interp;

% Format
title('Anomalous current','interpreter','latex')
xlabel('$T$ [keV]','interpreter','latex')
ylabel('$Z$','interpreter','latex')
zlabel('$e \Gamma^{an}$ [A/m$^{2}$]','interpreter','latex')

% Save as .pdf file
savepdf('APXAS2_anomalous_3d')

%% Plot G{an}(T = 1e3 [eV])
Zp      = min(Z):.01:max(Z);
G_anp   = interp1(Z,G_an(:,T==Tp),Zp,'spline');
figure; plot(Zp,-G_anp)

% Format
title(['Anomalous current (T = ' num2str(Tp/1e3) ' [keV])'],'interpreter','latex')
xlabel('$Z$ [keV]','interpreter','latex')
ylabel('$e \Gamma^{an}$ [A/m$^{2}$]','interpreter','latex')

% Save as .pdf file
savepdf('APXAS2_anomalous_2d')

%% Orbit loss flux

for j = 1:length(T)
    for i = 1:length(Z)
        g_ol(i,j)   = e*n*nu_eff(T(j))*sqrt(aspect)*rho_pi(T(j));
        G_ol(i,j)   = g_ol(i,j)*exp( - (nu_efi(T(j))+Z(i)^4)^(1/2) ) / ...
                                   (   (nu_efi(T(j))+Z(i)^4)^(1/2) );
    end
end

%% Plot mesh plot: (T,Z,G{ol})
[Tm,Zm] = meshgrid(T,Z);
figure; mesh(Tm*1e-3,Zm,G_ol);
colorbar; view([-45 45]); shading interp;

% Format
title('Orbit loss current','interpreter','latex')
xlabel('$T$ [keV]','interpreter','latex')
ylabel('$Z$','interpreter','latex')
zlabel('$e \Gamma^{ol}$ [A/m$^{2}$]','interpreter','latex')

% Save as .pdf file
savepdf('APXAS2_orbit_loss_3d')

%% Plot G{ol}(T = 1e3 [eV])
Zp      = min(Z):.01:max(Z);
G_olp   = interp1(Z,G_ol(:,T==Tp),Zp,'spline');
figure; plot(Zp,-G_olp)

% Format
title(['Orbit loss current (T = ' num2str(Tp/1e3) ' [keV])'],'interpreter','latex')
xlabel('$Z$ [keV]','interpreter','latex')
ylabel('$- e \Gamma^{ol}$ [A/m$^{2}$]','interpreter','latex')

% Save as .pdf file
savepdf('APXAS2_orbit_loss_2d')

%% Charge exchange flux

for j = 1:length(T)
    for i = 1:length(Z)
        g_cx(i,j)   = -( m_i*n_0(T(j))*R_cx*n*T(j) )/(Bp^2);
        G_cx(i,j)   = g_cx(i,j)*( 1/L_n + a_cx/L_T - Z(i)/rho_pi(T(j)) );
    end
end

%% Plot mesh plot: (T,Z,G{cx})
[Tm,Zm] = meshgrid(T,Z);
figure; mesh(Tm*1e-3,Zm,-G_cx);
colorbar; view([-45 45]); shading interp;

% Format
title('Charge exchange current','interpreter','latex')
xlabel('$T$ [keV]','interpreter','latex')
ylabel('$Z$','interpreter','latex')
zlabel('$- e \Gamma^{cx}$ [A/m$^{2}$]','interpreter','latex')

% Save as .pdf file
savepdf('APXAS2_charge_exchange_3d')

%% Plot G{cx}(T = 1e3 [eV])
Zp      = min(Z):.01:max(Z);
G_cxp   = interp1(Z,G_cx(:,T==Tp),Zp,'spline');
figure; plot(Zp,-G_cxp)

% Format
title(['Charge exchange current (T = ' num2str(Tp/1e3) ' [keV])'],'interpreter','latex')
xlabel('$Z$ [keV]','interpreter','latex')
ylabel('$- e \Gamma^{cx}$ [A/m$^{2}$]','interpreter','latex')

% Save as .pdf file
savepdf('APXAS2_charge_exchange_2d')
