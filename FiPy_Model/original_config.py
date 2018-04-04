"""
	This is an example configuration file. Read the other
	source files to see what values get defaulted.
"""

# Particle and heat flux from the core
Gamma_c = -0.8
q_c = 5.0*Gamma_c

# Number of cells
nx = 500

# Boolean, to choose either the original numerical model, or
# the full flux model. Note that it sets the length of the domain
# L to be 4.0 AU in the original, and 0.03 m in the flux model.
original_model = True

res_tol = 1.0e-6

# Total time steps; should be ~ L^2 / D
total_timeSteps = 2000

# Size of time step denominator (delta t)
# Either mu or epsilon is in the numerator
timeStep_denom = 15.0

# Choose the Diffusivity model, as a string (case does not matter)
# D_Zohm, D_Staps, and D_Shear are the possibilities
D_choice = "D_Shear"
# Coefficient of (Z')**beta in Stap's diffusivity
alpha_sup = 0.5
# Exponent of Z' in Stap's diffusivity
beta = 2.0

# Boolean, to choose what mode to have as initial conditions
initial_H_mode = False

# Boolean, to show the initial conditions
show_initial = False

# Choose numerical values for non-gradient Z-equation
# Currently, choices are Staps, Paquay, and some variant of g_grad
numerical_choice = "Staps"

# Plot details
plot_title = "H--Mode Start; $D \sim (Z^\prime)^{{-{:01.2f}}}$".format(beta)\
		+"\n" + r"$\Gamma_c = {:.3g},\, T \,=\, {:d},\,$"\
		.format(Gamma_c, total_timeSteps)\
		+ r"$\Delta t \,=\, \frac{{\epsilon}}{{{:2.0f}}}$"\
		.format(timeStep_denom)

# Maximum x on the plots
ploty_max = None
aux_plots = False
# Maximum and minimum y values on the plots
aux1y_min = None
aux1y_max = None
aux2y_min = None
aux2y_max = None

# Should the results be saved to file?
save_plots = False
save_TSVs = False

save_directory = ""

