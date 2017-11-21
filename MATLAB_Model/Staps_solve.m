function Staps_solve

global L epsilon mu zeta c_n c_T D_min D_max a1 a2 a3 g1 g2 g3 Z_S Gamma_c gamma q_c lambda_n lambda_T lambda_Z
L = 5.0;
epsilon = 1.0/25.0;
mu = 1.0/20.0;
zeta = 0.5;
c_n = -1.1;
c_T = -0.9;

D_min = 2.0/5.0;
D_max = 2.0;

% Coefficient for Stap's diffusivity
alpha_sup = 0.5;

% Coefficients for the Flow-Shear model of diffusivity
a1 = 1;
a2 = 0;
a3 = 1;

% Values for the polynomial G; g1 = a, g2 = b, g3 = c
g1 = 1.5;
g2 = 2.0;
g3 = -1.0;
Z_S = -1.5;

Gamma_c = -4.0/5.0;
gamma = 5.0/3.0;
q_c = -4.0;
lambda_n = 1.25;
lambda_T = 1.5;
lambda_Z = 1.25;

x = linspace(0, 5.0, 10);
t = linspace(0, 1.0, 10);

sol = pdepe(0, @the_pde, @the_ic, @the_bc, x, t);
n = sol(:,:,1);
Temp = sol(:,:,2);
Z = sol(:,:,3);

figure;
surf(x, t, n);
title('Density');
xlabel('x');
ylabel('t');

figure;
surf(x, t, Temp);
title('Temperature');
xlabel('x');
ylabel('t');

figure;
surf(x, t, Z);
title('Electric Field');
xlabel('x');
ylabel('t');

%figure;
%plot(x,Z(end,:));
%title('A Test');

	function [c, f, s] = the_pde(x, t, u, DuDx)
	% Zohm's model for diffusivity
%	Diffusivity = 0.5*(D_max + D_min) + 0.5*(D_max - D_min)*tanh(u(3));
	% Staps' model for diffusivity
	Diffusivity = D_min + (D_max - D_min) / (1 + alpha_sup*(DuDx(3))^2)
	% Flow-Shear model for diffusivity
%	Diffusivity = D_min + (D_max - D_min) / (1 + a1*u(3)^2 + a2*u(3)*DuDx(3) + a3*(DuDx(3))^2);
	G = g1 + g2*(u(3) - Z_S) + g3*(u(3) - Z_S)^3;

	c = [1.0; 1.0; 1.0];

	f = [DuDx(1); (1.0 / zeta)*DuDx(2); (mu / epsilon)*DuDx(3)] * Diffusivity;

	s = [0; (zeta + 1.0)/zeta * Diffusivity / u(1) * DuDx(1) * DuDx(2); c_n * u(2) / (epsilon*(DuDx(1))^2) + c_T / (epsilon*u(1)*DuDx(2)) + G/epsilon];
	end

	function u0 = the_ic(x)
	density_initial = -Gamma_c*lambda_n * (1 + x / lambda_n);
	temperature_initial = (gamma - 1) / Gamma_c * (1 - (lambda_n / (zeta*lambda_T + lambda_n))*(1 + x / lambda_n)^(-zeta));
	u0 = [density_initial; temperature_initial; Z_S*(1.0 - tanh((L / 2.0)*(x - 1.0)))];
	end

	function [pl, ql, pr, qr] = the_bc(xl, ul, xr, ur, t)
	%pl = -[ul(1)/lambda_n; (1.0/zeta)*ul(2)/lambda_T; mu/epsilon*ul(3)/lambda_Z] * (D_min + (D_max - D_min) / (1 + a1*(ul(3))^2 + a2*ul(3)*(ul(3) / lambda_Z) + a3*(ul(3) / lambda_Z)^2));
	pl = -[ul(1)/lambda_n; (1.0/zeta)*ul(2)/lambda_T; mu/epsilon*ul(3)/lambda_Z] * (0.5*(D_max + D_min) + 0.5*(D_max - D_min)*tanh(ul(3)));
	ql = [1; 1; 1];
	pr = [Gamma_c; ((gamma - 1)*q_c - ur(2)*Gamma_c) / ur(1); 0];
	qr = [1; 1; 1];
	end

end

