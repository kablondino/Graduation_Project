reset()
# Import file with parameters, etc.
#load("/home/kabv/Documents/Masters/Graduation_Project/Flux_Calculations/parameters.sage")

var('x')
Z = function('Z')(x)
D = function('D')(x)
n = function('n')(x)

L = 5.0
D_min = 2.0/5.0
D_max = 2.0
alpha_sup = 0.5
Z_S = 0.5

# TO BE FILLED IN
lambda_n = ??
Gamma_c = ??

#var('D_min,D_max,alpha_sup,Z_S,L')

Z = Z_S*(1 - tanh((L/2)*(x - 1)))

D = D_min + (D_max - D_min) / (1 + alpha_sup*(derivative(Z,x))^2)

# To print or not to print
#to_print = False

# Print out the function's results
def to_print():
	if to_print == True:
		print "Z = "+str(Z.full_simplify())
		print "dZ/dx = "+str(derivative(Z,x).full_simplify())
		print "d^2Z/dx^2 = "+str(derivative(Z,x,2).full_simplify())
		print "D = "+str(D.full_simplify())
		print "dD/dx = "+str(derivative(D,x).full_simplify())

density_initial_cond = derivative( D*derivative(n,x), x ) == 0
f = desolve( density_initial_cond, n, ivar=x, show_method=True )
pretty_print( f )
#pretty_print(density_initial_cond)

## ----------------- Plot everything! ---------------------
#plotting_colors = ['red','green','blue', 'magenta', 'cyan', 'yellow']
#
#Z_title = 'Initial Condition of $Z(x) = Z_S\\left[1 - \\tanh\\left(\\frac{Lx - L}{2}\\right)\\right]$ and its Derivatives'
#D_title = 'Diffusivity $D(\partial_x Z) = D_{\\tt{min}} + \\frac{D_{\\tt{max}} - D_{\\tt{min}}}{1 + \\alpha_{\\tt{sup}}\\left(\partial_x Z\\right)^2}$ with $Z(x,0) = Z_S\\left[1 - \\tanh\\left(\\frac{Lx - L}{2}\\right)\\right]$'
#
#Z_functions_plots = sum([plot(derivative(Z, x, i).subs(Z_S = 0.5, L = 5.0), (x , 0, 3.0), gridlines=True, legend_label='$\partial_x^'+str(i)+'Z$', color=plotting_colors[i], title=Z_title) for i in range(3)])
#Z_functions_plots.set_legend_options(font_size=14, loc='lower right')
#
#Z_parameter_box = text('$Z_S = 0.5$ \n $L = 5.0$', (3.0,-1.5), bounding_box={'boxstyle':'round', 'fc':'w'}, fontsize=12, color='black', horizontal_alignment='right')
#
##show(Z_functions_plots + Z_parameter_box)
#
#D_functions_plots = sum([plot(derivative(D, x, i).subs(D_min = 2/5, D_max = 2, alpha_sup = 0.5, Z_S = 0.5, L = 5.0), (x, 0, 3.0), gridlines=True, legend_label='$\partial_x^'+str(i)+'D$', color=plotting_colors[i+3], title=D_title) for i in range(2)])
#D_functions_plots.set_legend_options(font_size=14, loc='lower right')
#
#D_parameter_box = text("$Z_S = 0.5$ \n $L = 5.0$ \n $D_{\\tt{min}} = 0.4$ \n $D_{\\tt{max}} = 2.0$ \n $\\alpha_{\\tt{sup}} = 0.5$", (2.65,-1.7), bounding_box={'boxstyle':'round', 'fc':'w'}, fontsize=12, color='black', horizontal_alignment='right')

#show(Z_functions_plots + D_functions_plots + D_parameter_box)

