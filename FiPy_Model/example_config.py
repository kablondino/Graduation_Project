"""
	This is an example configuration file.
"""

# Boolean, to choose what mode to have as initial conditions
Initial_H_mode = True

# Number of cells
nx = 100

# Choose the Diffusivity model, as a string (case does not matter)
# D_Zohm, D_Staps, and D_Shear are the possibilities
D_choice = "D_Staps"

# Choose numerical values for non-gradient Z-equation
# Currently, choices are Staps, Paquay, and some variant of g_grad
numerical_choice = "Staps"

