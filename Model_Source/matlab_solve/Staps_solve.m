function Staps_solve

global epsilon mu zeta c_n c_T D_min D_max alpha_sup g1 g2 g3 Z_S Gamma_c gamma q_c lambda_n lambda_T lambda_Z% D_shear
epsilon = 1.0/25.0;
mu = 1.0/20.0;
zeta = 0.5;
c_n = -1.1;
c_T = -0.9;

D_min = 2.0/5.0;
D_max = 2.0;
alpha_sup = 0.5;

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

x = linspace(0, 1, 100);
t = linspace(0, 1, 100);

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

	function [c, f, s] = the_pde(x, t, u, DuDx)
	D_shear = D_min + (D_max - D_min) / (1 + alpha_sup*(DuDx(3))^2);
	G = g1 + g2*(u(3) - Z_S) + g3*(u(3) - Z_S)^3;

	c = [1.0; 1.0; 1.0];

	f = [DuDx(1); (1.0 / zeta)*DuDx(2); (mu / epsilon)*DuDx(3)] * D_shear;

	s = [0; (zeta + 1.0)/zeta * D_shear / u(1) * DuDx(1) * DuDx(2); c_n * u(2) / (epsilon*(DuDx(1))^2) + c_T / (epsilon*u(1)*DuDx(2)) + G/epsilon];
	end

	function u0 = the_ic(x)
%	u0 = [-(Gamma_c*lambda_n)/(D_shear(0)); 2.0e3*x; Z_S*(1.0 - tanh(0.5*(x - 1.0)))];
	u0 = [1.0e19*x; 2.0e3*x; Z_S*(1.0 - tanh(0.5*(x - 1.0)))];
	end

	function [pl, ql, pr, qr] = the_bc(xl, ul, xr, ur, t)
	pl = -[ul(1)/lambda_n; (1.0/zeta)*ul(2)/lambda_T; mu/epsilon*ul(3)/lambda_Z];% * D_shear;
	ql = [1; 1; 1];
	pr = [Gamma_c; ((gamma - 1)*q_c - ur(2)*Gamma_c) / ur(1); 0];
	qr = [1; 1; 1];
	end

end

