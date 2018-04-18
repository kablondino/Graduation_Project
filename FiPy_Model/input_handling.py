"""
	This file deals with all of the inputs for the system. The current set
	of inputs are the following:
	nx:				int		The number of grid points
	total_timeSteps:int		The total number of time steps
	timeStep_denom:	float	The denominator of the time step; DEPRICATED
	timeStep:		float	The overall dt in solving
	res_tol			float	The tolerance of the residual
	Gamma_c			float	The particle flux from the core
	q_c				float	The heat flux from the core
	numerical_parameter:	string		The set of predetermined parameters
	D_choice:		string	The model of the diffusivity
	initial_H_mode:	bool	Start in L-- or H--mode?
	original_model:	bool	Which Z-equation model? 'False' is flux model.
	view_matplotlib	bool	Should use base (non-FiPy) matplotlib?
	show_initial:	bool	Show the initial conditions
	plot_title:		string	The title of the plot; can be formatted
	ploty_max:		float	The maximum y-value on the plot
	aux_plots:		bool	Turns on specified auxiliary plots
	aux_vars:		list	List of strings for auxiliary plots
	aux_titles:		list	Title of the aux plots
	aux_ymin:		float	Minimum y value of the aux plots
	aux_ymax:		float	Maximum y value of the aux plots

	save_directory: string	The name of the saving directory, from current
							directory being run.
	save_plots:		bool	Should the plots be saved?
	save_TSVs:		bool	Should TSV files be generated and saves?

	Each possible input also has a default value, if nothing is set.

	NOTE that the auxiliary variables 'aux_vars' is only checked for
	variable type (list of strings) in this file. The
	validity of the contents is checked in the solving file.
"""
import sys

# Import variables from job configuration file, as cmd line argument
config = __import__ (sys.argv[1].replace('.py',''))

parameter_sets = ["staps", "paquay", "g_grad", "gradient_model"]
diffusivity_models = ["d_zohm", "d_staps", "d_shear", "d_flow_shear",\
		"d_weymiens_l"]


# -------------- Check Configuration Variables ------------
# Z-equation model choice
if type(getattr(config, 'original_model', None)) != bool:
	config.original_Z_model = True
	print "Defaulted to using the original numerical model for Z."


# Particle and heat fluxes from the core
if (type(getattr(config, 'Gamma_c', None)) != float\
		and type(getattr(config, 'Gamma_c', None)) != int):
	try:
		config.Gamma_c = int(input("The particle flux from the core Gamma_c is not chosen properly. Choose a floating-point value: "))
	except (EOFError, NameError, SyntaxError) as e:
		config.Gamma_c = -4.0 / 5.0
		print "Gamma_c defaulted to -0.8"

if (type(getattr(config, 'q_c', None)) != float\
		and type(getattr(config, 'Gamma_c', None)) != int):
	try:
		config.Gamma_c = int(input("The heat flux from the core q_c is not chosen properly. Choose a floating-point value: "))
	except (EOFError, NameError, SyntaxError) as e:
		config.q_c = config.Gamma_c * 5.0
		print "q_c defaulted to 5.0 * Gamma_c"


# Check the choice for numerical parameters
if (getattr(config, 'numerical_choice', "").lower() not in parameter_sets\
		or type(getattr(config, 'numerical_choice', None)) != str):
	try:
		config.numerical_choice =\
			raw_input("The set of parameters not properly chosen. Choose from the following: Staps, Paquay -> ")

		if config.numerical_choice.lower() not in parameter_sets:
			raise TypeError

	except (EOFError, TypeError):
		config.numerical_choice = "Paquay"
		print "Numerical choice defaulted to Paquay's set."


# Grid points
if ((type(getattr(config, 'nx', None)) != int and\
		type(getattr(config, 'nx', None)) != float) or\
		getattr(config, 'nx', None) <= 0):
	try:
		config.nx = int(input("nx (Grid number) not properly defined. Enter a positive integer value: "))

		if config.nx <= 0:
			raise ValueError

		print "nx set to " + str(config.nx)

	except (EOFError, NameError, SyntaxError, ValueError) as e:
		config.nx = 100
		print "nx defaulted to 100."

if type(config.nx) == float:
	config.nx = int(config.nx)

# Domain size		NOW DEPRICATED!
#if ((type(getattr(config, 'L', None)) != float and\
#		type(getattr(config, 'L', None)) != int) or\
#		getattr(config, 'L', None) <= 0.0):
#	try:
#		config.L = float(input("Length of domain not properly defined. Enter floating-point value: "))
#
#		if config.L <= 0:
#			raise NameError
#
#		print "L set to " + str(config.L)
#
#	except (NameError, SyntaxError, EOFError, ValueError):
#		config.L = 4.0
#		print "L defaulted to 4.0"
#
#if type(config.L) == int:
#	config.L = float(config.L)


# Total number of time steps
if ((type(getattr(config, 'total_timeSteps', None)) != int and\
		type(getattr(config, 'total_timeSteps', None)) != float) or\
		getattr(config, 'total_timeSteps', None) <= 0.0):
	try:
		config.total_timeSteps = int(input("Total number of time steps not properly defined. Enter integer value: "))

		if config.total_timeSteps <= 0:
			raise ValueError

		print "Total # of time steps set to " + str(config.total_timeSteps)

	except (NameError, SyntaxError, EOFError, ValueError):
		config.total_timeSteps = 100
		print "Total time steps defaulted to 100."

if type(config.total_timeSteps) == float:
	config.total_timeSteps = int(config.total_timeSteps)


# Denominator of the delta t
# DEPRICATED!
#if ((type(getattr(config, 'timeStep_denom', None)) != float and\
#		type(getattr(config, 'timeStep_denom', None)) != int) or\
#		getattr(config, 'timeStep_denom', None) <= 0.0):
#	try:
#		config.timeStep_denom = float(input("The denomintor of the time step size is not properly defined. Enter floating-point value: "))
#
#		if config.timeStep_denom <= 0.0:
#			raise ValueError
#
#		print "The denominator of the time step size is set to "\
#				+ str(config.timeStep_denom)
#
#	except (NameError, SyntaxError, EOFError, ValueError):
#		config.timeStep_denom = 15.0
#		print "The denominator of the time step size is defaulted to 15.0"
#
#if type(config.timeStep_denom) == int:
#	config.timeStep_denom = float(config.timeStep_denom)


# Time step
if ((type(getattr(config, 'timeStep', None)) != float and\
		type(getattr(config, 'timeStep', None))) or\
		getattr(config, 'timeStep', None) <= 0.0):
	try:
		config.timeStep_denom = float(input("The time step size is not properly defined. Enter floating-point value: "))

		if config.timeStep <= 0.0:
			raise ValueError

		print "The time step size is set to " + str(config.timeStep)

	except (NameError, SyntaxError, EOFError, ValueError):
		if config.original_model == True:
			config.timeStep = 1.0 / 375.0
			print "The time step is defaulted to 1.0 / 375.0."
		elif config.original_model == False:
			config.timeStep = 1.0e-9
			print "The time step is defaulted to 1.0e-9."

if type(config.timeStep) == int:
	config.timeStep = float(config.timeStep)


# Residual tolerance
if ((type(getattr(config, 'res_tol', None)) != float and\
		type(getattr(config, 'res_tol', Non)) != int) or\
		getattr(config, 'res_tol', None) <= 0.0):
	try:
		config.res_tol = float(input("The residual tolerance is not properly set. Enter a positive floating-point value: "))

		if config.res_tol <= 0.0:
			raise ValueError

		print "The residual tolerance is defaulted to " +str(config.res_tol)

	except (NameError, SyntaxError, EOFError, ValueError):
		if config.original_model == True:
			config.res_tol = 1.0e-6
			print "The residual tolerance is defaulted to 1.0e-6."
		elif config.original_model == False:
			config.res_tol = 1.0e14
			print "The residual tolerance is defaulted to 1.0e14."

if type(config.res_tol) == int:
	config.res_tol = float(config.res_tol)


# Choice of the diffusivity model
if (getattr(config, 'D_choice', "").lower() not in diffusivity_models or\
		type(getattr(config, 'D_choice', None)) != str):
	try:
		config.D_choice = raw_input("The diffusivity model is not properly chosen. Choose from the following: Zohm, Weymiens_L, Staps, Shear -> ")

		if config.D_choice.lower() not in diffusivity_models:
			raise IndexError()

	except (IndexError, EOFError):
		config.D_choice = "d_staps"
		print "Diffusivity model defaulted to Staps'."


# Initial starting mode
if type(getattr(config, 'initial_H_mode', None)) != bool:
	config.initial_H_mode = False
	print "Defaulted to starting in L--mode."


# Show initial conditions
if type(getattr(config, 'show_initial', None)) != bool:
	config.show_initial = False
	print "The initial conditions will not be shown."


# Matplotlib option
if type(getattr(config, 'view_matplotlib', None)) != bool:
	config.view_matplotlib = False
	print "The default FiPy viewers will be used."


# Plot title
if not hasattr(config, 'plot_title'):
	config.plot_title = ""


## Auxiliary plots
# Check for type
if (type(getattr(config, 'aux_plots', None)) != bool or\
		type(getattr(config, 'aux_vars', None)) != list):
	config.aux_plots = False

if config.aux_plots == True:
	# Create aux_titles, _ymin, and _ymax lists if they don't exist
	if not hasattr(config, 'aux_titles'):
		config.aux_titles = []
	if not hasattr(config, 'aux_ymin'):
		config.aux_ymin = []
	if not hasattr(config, 'aux_ymax'):
		config.aux_ymax = []
	
	# Forces all variable calls and titles to strings
	if all(isinstance(i, str) for i in config.aux_vars) == False:
		for i in range(len(config.aux_vars)):
			config.aux_vars[i] = str(config.aux_vars[i])
	if all(isinstance(i, str) for i in config.aux_titles) == False:
		for i in range(len(config.aux_titles)):
			config.aux_titles[i] = str(config.aux_titles[i])

	# Make the aux_titles, _ymin, and _ymax lists long enough
	while len(config.aux_vars) > len(config.aux_titles):
		config.aux_titles.append(None)
	while len(config.aux_vars) > len(config.aux_ymin):
		config.aux_ymin.append(None)
	while len(config.aux_vars) > len(config.aux_ymax):
		config.aux_ymax.append(None)

	for j in range(len(config.aux_vars)):
		if type(config.aux_titles[j]) != str:
			config.aux_titles[j] = None

		# If the datamins/maxes for aux plots are bad, set them to None
		if (type(config.aux_ymin[j]) != float and\
				type(config.aux_ymin[j]) != int):
			config.aux_ymin[j] = None
		if (type(config.aux_ymax[j]) != float and\
				type(config.aux_ymax[j]) != int):
			config.aux_ymax[j] = None


# Makes sure that the saved directory is a string
if hasattr(config, 'save_directory'):
	if config.save_directory != str:
		config.save_directory = str(config.save_directory)


# If saving data is enabled, but not a directory, exit the run.
if (getattr(config, 'save_directory', None) == None and\
		(getattr(config, 'save_plots', False) == True or\
		getattr(config, 'save_TSVs', False) == True)):
	sys.exit("No directory specified for saving specified files. Exiting...")


# Assumes save_directory exists, but not written correctly as a string
if hasattr(config, 'save_directory'):
	if type(config.save_directory) != str:
		config.save_directory = str(config.save_directory)

# If save_plots and/or TSVs does not exist or not booleans, set to False
if (not hasattr(config, 'save_plots') or\
		type(getattr(config, 'save_plots', None)) != bool):
	config.save_plots = False
if (not hasattr(config, 'save_TSVs') or\
		type(getattr(config, 'save_TSVs', None)) != bool):
	config.save_TSVs = False

