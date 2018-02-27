"""
	This is an example configuration file. Read the other
	source files to see what values get defaulted.
"""

# Number of cells
nx = 200

# Length of domain
L = 4.0

# Total time steps; should be ~ L^2 / D
total_timeSteps = 2000

# Size of time step denominator (delta t)
# Either mu or epsilon is in the numerator
timeStep_denom = 10.0

# Choose the Diffusivity model, as a string (case does not matter)
# D_Zohm, D_Staps, and D_Shear are the possibilities
D_choice = "D_Staps"

# Boolean, to choose what mode to have as initial conditions
initial_H_mode = True

# Choose numerical values for non-gradient Z-equation
# Currently, choices are Staps, Paquay, and some variant of g_grad
numerical_choice = "Staps"

# Choose the largest acceptable residual when sweeping
res_tol = 1.0e-6

# Plot details
plot_title = "GMRES H--Mode Start; Staps' parameters \n $T = $"\
		+str(total_timeSteps)+ r", $\Delta t = \mu / $"+str(timeStep_denom)

# Maximum x on the plot
plot_max = L

# Should the results be saved to file?
save_plots = False
save_TSVs = False

save_directory = "Test"

