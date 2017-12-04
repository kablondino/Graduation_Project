## File to (attempt to) calculate the initial conditions for density and tempurature
# It also plots Z and the selected models of diffusivity

reset()
# Import file with parameters, etc.
#load("/home/kabv/Documents/Masters/Graduation_Project/Flux_Calculations/parameters.sage")

var('x')
Z = function('Z')(x)
D = function('D')(x)
n = function('n')(x)

# Parameters
L = 5.0
zeta = 0.5
Gamma_c = -4.0/5.0
q_c = -4.0
gamma = 5.0/3.0

lambda_n = 1.25
lambda_T = 3.0/2.0

D_min = 2.0/5.0
D_max = 2.0
alpha_sup = 0.5

Z_S = -1.5

# Coefficients for the Flow-Shear model's denominator
a1 = 0.7
a2 = 1.25
a3 = 0.5

# Initial condition for Z used by Staps
Z(x) = Z_S*(1.0 - tanh((L/2.0)*(x - 1.0)))

diff_array = ["Zohm", "Staps", "Flow-Shear"]

diff_selector = diff_array[2]

def diffusivity(aDiff_selector):
	# Zohm model of diffusivity; considered non-realistic
	if aDiff_selector.lower() == "zohm":
		return ((D_max + D_min) / 2.0 + ((D_max + D_min)*tanh(Z)) / 2.0)

	# Staps' model of diffusivity (not dependent on Z, but only on dZ/dx)
	elif aDiff_selector.lower() == "staps":
		return D_min + (D_max - D_min) / (1.0 + alpha_sup*(derivative(Z,x))**2)

	# Flow-Shear model (dependent on both Z and dZ/dx)
	elif aDiff_selector.lower() == "flow" or "flow-shear" or "flowshear" or "shear":
		return D_min + (D_max - D_min) / (1.0 + a1*Z**2.0 + a2*Z*derivative(Z,x) + a3*(derivative(Z,x))**2)

	else:
		print "Diffusivity form not properly selected."

def diffusivity_title(aDiff_selector):
	# Zohm model of diffusivity; considered non-realistic
	if aDiff_selector.lower() == "zohm":
		return r"Zohm's Model: $D \,=\, \frac{D_{max} + D_{min}}{2} + \frac{(D_{max} - D_{min})\tanh(Z)}{2}$"

	# Staps' model of diffusivity (not dependent on Z, but only on dZ/dx)
	elif aDiff_selector.lower() == "staps":
		return r"Staps' Model: $D \,=\, D_{min} + \frac{D_{max} - D_{min}}{1 + \alpha_{sup} \cdot (Z^{\prime})^2}$"

	# Flow-Shear model (dependent on both Z and dZ/dx)
	elif aDiff_selector.lower() == "flow" or "flow-shear" or "flowshear" or "shear":
		return r"Flow-Shear Model: $D \,=\, D_{min} + \frac{D_{max} - D_{min}}{1 + a_1 Z^2 + a_2 Z \cdot Z^{\prime} + a_3 (Z^{\prime})^2}$"

	else:
		print "Diffusivity form not properly selected."

# ----------------- Initial Conditions ---------------------
if diff_selector.lower() == "zohm":
	density_initial(x) = 1.0e-2*-(Gamma_c*lambda_n / diffusivity(diff_selector))* (1 + x/lambda_n)
else:
	density_initial(x) = -(Gamma_c*lambda_n / diffusivity(diff_selector))* (1 + x/lambda_n)

temp_initial(x) = q_c * ((gamma - 1.0) / Gamma_c) * (1.0 - (lambda_n / (zeta*lambda_T + lambda_n)) * (1.0 + x / lambda_n)**(-zeta))


# Print out the function's results
def to_print():
	if toggle_print == True:
		print "Z = "+str(Z.full_simplify())
		print "dZ/dx = "+str(derivative(Z,x).full_simplify())
		print "d^2Z/dx^2 = "+str(derivative(Z,x,2).full_simplify())
		print "D = "+str(D.full_simplify())
		print "dD/dx = "+str(derivative(D,x).full_simplify())


## ----------------- Plot everything! ---------------------
plotting_colors = ['red', 'green', 'blue', 'magenta', 'cyan', 'yellow']

the_title = r"Initial Conditions with " + diffusivity_title(diff_selector)
plot_min, plot_max = 0.0, 3.0

density_plot = plot(density_initial, (x, plot_min, plot_max), gridlines=True, legend_label='$n_i$', color=plotting_colors[2])

temp_plot = plot(temp_initial, (x, plot_min, plot_max), gridlines=True, legend_label='$T_i$', color=plotting_colors[1])

Z_plot = plot(Z, (x, plot_min, plot_max), gridlines=True, legend_label='$Z_i$', color=plotting_colors[0])
shear_plot = plot(Z.derivative(x), (x, plot_min, plot_max), gridlines=True, legend_label=r'$\frac{\partial Z_i}{\partial x}$', color=plotting_colors[5])
shear_plot.set_legend_options(loc='lower right')


D_plot = plot(diffusivity(diff_selector), (x, plot_min, plot_max), gridlines=True, legend_label='$D$', color=plotting_colors[3])

#show(density_plot + temp_plot + Z_plot + D_plot + shear_plot, title=the_title)

# Plot each diffusivity form
D_plots = sum([plot(diffusivity(diff_array[i]), (x, plot_min, plot_max), gridlines=True, legend_label=diff_array[i], color=plotting_colors[i]) for i in range(3)])
#show(D_plots, title= "Diffusivity Models")

# ----------------- ANIMATIONS ----------------------------
# Loop over diffusivity parameters
#D_shear_plots = plot(diffusivity("flow-shear"), (x, plot_min, plot_max), gridlines=True, legend_label= r'$a_1$')
D_shear_plots_a1 = []
import numpy
for j in numpy.arange(0.8, 2.1, 0.1):
	a1 = j
	D_shear_plots_a1.append(plot(diffusivity("flow-shear"), (x, plot_min, plot_max), legend_label="$a_1 = $"+str(a1), ymin=0.0, gridlines=True))#, color=plotting_colors[numpy.random.randint(0,0)]))
#show(sum(D_shear_plots), title="Flow-Shear Diffusivity", gridlines=True)
# Reset Coefficients
a1 = 1.0

D_shear_plots_a2 = []
for k in numpy.arange(0.4, 1.5, 0.1):
	a2 = k
	D_shear_plots_a2.append(plot(diffusivity("flow-shear"), (x, plot_min, plot_max), legend_label="$a_2 = $"+str(a2), ymin=0.0, gridlines=True))#, color=plotting_colors[numpy.random.randint(0,0)]))

a2 = 1.0

D_shear_plots_a3 = []
for l in numpy.arange(0.4, 1.9, 0.1):
	a3 = l
	D_shear_plots_a3.append(plot(diffusivity("flow-shear"), (x, plot_min, plot_max), legend_label="$a_3 = $"+str(a3), ymin=0.0, gridlines=True))#, color=plotting_colors[numpy.random.randint(0,0)]))

#animate(D_shear_plots_a2).show()

# ----------------- INTERACT ------------------------------
#@interact
#def _(a1=(0.8, 2.1, step_size=0.1)):
#	show(plot(diffusivity()))

## ---------------- OLD -----------------------------------
#Z_title = 'Initial Condition of $Z(x) = Z_S\\left[1 - \\tanh\\left(\\frac{Lx - L}{2}\\right)\\right]$ and its Derivatives'
#D_title = 'Diffusivity $D(\partial_x Z) = D_{\\tt{min}} + \\frac{D_{\\tt{max}} - D_{\\tt{min}}}{1 + \\alpha_{\\tt{sup}}\\left(\partial_x Z\\right)^2}$ with $Z(x,0) = Z_S\\left[1 - \\tanh\\left(\\frac{Lx - L}{2}\\right)\\right]$'

#Z_functions_plots = sum([plot(derivative(Z, x, i).subs(Z_S = 0.5, L = 5.0), (x , 0, 3.0), gridlines=True, legend_label='$\partial_x^'+str(i)+'Z$', color=plotting_colors[i], title=Z_title) for i in range(3)])
#Z_functions_plots.set_legend_options(font_size=14, loc='lower right')

#Z_parameter_box = text('$Z_S = 0.5$ \n $L = 5.0$', (3.0,-1.5), bounding_box={'boxstyle':'round', 'fc':'w'}, fontsize=12, color='black', horizontal_alignment='right')

#show(Z_functions_plots + Z_parameter_box)

#D_functions_plots = sum([plot(derivative(D, x, i).subs(D_min = 2/5, D_max = 2, alpha_sup = 0.5, Z_S = 0.5, L = 5.0), (x, 0, 3.0), gridlines=True, legend_label='$\partial_x^'+str(i)+'D$', color=plotting_colors[i+3], title=D_title) for i in range(2)])
#D_functions_plots.set_legend_options(font_size=14, loc='lower right')

#D_parameter_box = text("$Z_S = 0.5$ \n $L = 5.0$ \n $D_{\\tt{min}} = 0.4$ \n $D_{\\tt{max}} = 2.0$ \n $\\alpha_{\\tt{sup}} = 0.5$", (2.65,-1.7), bounding_box={'boxstyle':'round', 'fc':'w'}, fontsize=12, color='black', horizontal_alignment='right')

#show(D_functions_plots + D_parameter_box)

