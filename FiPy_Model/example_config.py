"""
	This is an example configuration file. Read the other
	source files to see what values get defaulted.
"""

# Particle and heat flux from the core
Gamma_c = -4.0 / 5.0
q_c = -4.0

# Number of cells
nx = 1000

# Length of domain
L = 4.0

# Total time steps; should be ~ L^2 / D
total_timeSteps = 200

# Size of time step denominator (delta t)
# Either mu or epsilon is in the numerator
timeStep_denom = 20.0

# Choose the Diffusivity model, as a string (case does not matter)
# D_Zohm, D_Staps, and D_Shear are the possibilities
D_choice = "D_Staps"
# Coefficient of (Z')**beta in Stap's diffusivity
alpha_sup = 0.5
# Exponent of Z' in Stap's diffusivity
beta = 1.5

# Boolean, to choose what mode to have as initial conditions
initial_H_mode = True

# Boolean, to choose either the original numerical model, or
# the full flux model
original_model = True

# Choose numerical values for non-gradient Z-equation
# Currently, choices are Staps, Paquay, and some variant of g_grad
numerical_choice = "Staps"

# Plot details
plot_title = "H--Mode Start; $D \sim (Z^\prime)^{{{:01.2f}}}$".format(beta)\
		+"\n" + r"$\Gamma_c = {:.2f},\, T \,=\, {:d},\,$"\
		.format(Gamma_c, total_timeSteps)\
		+ r"$\Delta t \,=\, \frac{{\epsilon}}{{{:02.0f}}}$"\
		.format(timeStep_denom)

# Maximum x on the plots
plotx_max = L
ploty_max = None
# Maximum and minimum y values on the plots
aux1y_min = None
aux1y_max = None
aux2y_min = None
aux2y_max = None

# Should the results be saved to file?
save_plots = False
save_TSVs = False

save_directory = ""

