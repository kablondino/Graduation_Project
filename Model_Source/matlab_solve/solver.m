function pde_model

m = 0;
x = linspace(0, 1, 100);
t = linspace(0, 2, 50);

sol = pdepe(m, @the_pde, @the_ic, @the_bc, x, t);
Z = sol(:,:,1);

surf(x, t, Z)
title('The first attempt')
xlabel('x')
ylabel('t')

figure
plot(x, Z(end,:))
title('Solution at t = 2')
xlabel('x')
ylabel('Z(x, 2)')

% --------------------------------------------------------------
% The model PDE itself
function [c, f, s] = the_pde(x, t, Z, DZDx)
c = symfun( m_i*B_p^2 / (charge*B^2) * n(x) * Temp(x), x );
f = symfun( m_i / (e*rho_pi(x)) * n(x) * Temp(x) * DZDx, [x Z] );
s = symfun( B_p^2 * (g_n_sum*L_n + g_T_sum*L_T + g_Z_sum*Z), [x Z] );

% --------------------------------------------------------------
% Initial condition
function Z0 = the_ic(x)
Z0 = 0;

% --------------------------------------------------------------
% Boundary conditions
function [pl, ql, pr, qr] = the_bc(xl, Zl, xr, Zr, t)
pl = 0;
ql = 0;
pr = Zr;
qr = 1;
