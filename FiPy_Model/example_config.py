"""
	This is an example configuration file. Read the other
	source files to see what values get defaulted.
"""

# Particle and heat flux from the core
Gamma_c = -4.0 / 5.0
q_c = -4.0
print type(Gamma_c)

# Number of cells
nx = 200

# Length of domain
L = 4.0

# Total time steps; should be ~ L^2 / D
total_timeSteps = 200

# Size of time step denominator (delta t)
# Either mu or epsilon is in the numerator
timeStep_denom = 12.0

# Choose the Diffusivity model, as a string (case does not matter)
# D_Zohm, D_Staps, and D_Shear are the possibilities
D_choice = "D_Staps"

# Boolean, to choose what mode to have as initial conditions
initial_H_mode = False

# Boolean, to choose what units are used
SI_units = False

# Choose numerical values for non-gradient Z-equation
# Currently, choices are Staps, Paquay, and some variant of g_grad
numerical_choice = "Staps"

# Choose the largest acceptable residual when sweeping
res_tol = 1.0e-6

# Plot details
plot_title = "L--Mode Start; Staps' parameters; \n" +\
		r"$\Gamma_c =$" +str(Gamma_c)+ \
		r", $T = $"+str(total_timeSteps)+ r", $\Delta t = \mu / $"\
		+str(timeStep_denom)

# Maximum x on the plots
plotx_max = L
# Maximum and minimum y values on the plots
aux1y_min = None
aux1y_max = None
aux2y_min = None
aux2y_max = None

# Should the results be saved to file?
save_plots = True
save_TSVs = True

save_directory = "Standard_4_5"

