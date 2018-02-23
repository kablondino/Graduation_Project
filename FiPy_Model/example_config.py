"""
	This is an example configuration file.
	Every configuration file MUST include these variables,
	exactly.
"""

# Number of cells
nx = 500

# Length of domain
L = 5.0

# Total time steps; should be ~ L^2 / D
total_time = 300

# Size of time step denominator (delta t)
# Either mu or epsilon is in the numerator
timeStep_denom = 40.0

# Choose the Diffusivity model, as a string (case does not matter)
# D_Zohm, D_Staps, and D_Shear are the possibilities
D_choice = "D_Staps"

# Boolean, to choose what mode to have as initial conditions
initial_H_mode = True

# Choose numerical values for non-gradient Z-equation
# Currently, choices are Staps, Paquay, and some variant of g_grad
numerical_choice = "Staps"

# Plot details
plot_title = "GMRES H--Mode Start; $t = $"+str(total_time)+\
			r", $\Delta t = \epsilon / $"+str(timeStep_denom)

# Should the results be saved to file?
save_plots = False
save_TSVs = False

save_directory = "Test"

