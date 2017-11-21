syms density(x)

global L Z_S D_min D_max alpha_sup Z density
L = 5.0;
Z_S = -1.5;
D_min = 2.0/5.0;
D_max = 2.0;
alpha_sup = 0.5;

%x = linspace(0, L, 100);
Z = symfun( Z_S*(1.0 - tanh((L/2.0)*(x - 1.0))), x );
%D_shear = symfun( , x );
conds = [density(0) == 0 Z(0) == 0]

[V] = odeToVectorField( diff(D_min + (D_max - D_min) / (1 + alpha_sup*(diff(Z))^2) * diff(density)) == 0, conds )

%M = matlabFunction(V, 'vars', {'x','Y'})

%sol = ode15i( M, [0 L], [0 0] );

%fplot( @(x)deval(sol,x,1), [0 0], [0 0])

