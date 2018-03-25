"""
	This file deals with all of the inputs for the system. The current set
	of inputs are the following:
	nx:				int		The number of grid points
	L:				float	The length of the domain
	total_timeSteps:int		The total number of time steps
	timeStep_denom:	float	The denominator of the time step
	Gamma_c			float	The particle flux from the core
	q_c				float	The heat flux from the core
	numerical_parameter:	string		The set of predetermined parameters
	D_choice:		string	The model of the diffusivity
	initial_H_mode:	bool	Start in L-- or H--mode?
	original_model:	bool	Which Z-equation model? 'False' is flux model.
	SI_units		bool	Determines the units of the profiles
	plot_title:		string	The title of the plot; can be formatted
	plotx_max:		float	The maximum x-value on the plot (min is always 0)
	ploty_max:		float	The maximum y-value on the plot
	aux_plots:		bool	Turns on specified auxiliary plots
	aux_title1:		str		Title of the first auxiliary plot
	aux_title2:		str		............ second ...
	aux_plot_var:	tuple	Specified variables for the auxiliary plots

	save_directory: string	The name of the saving directory, from current
							directory being run.
	save_plots:		bool	Should the plots be saved?
	save_TSVs:		bool	Should TSV files be generated and saves?

	Each possible input also has a default value, if nothing is set.
"""
import sys

# Import variables from job configuration file, as cmd line argument
config = __import__ (sys.argv[1].replace('.py',''))

parameter_sets = ["staps", "paquay", "g_grad", "gradient_model"]
diffusivity_models = ["d_zohm", "d_staps", "d_shear", "d_flow_shear"]


# -------------- Check Configuration Variables ------------
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

# Domain size
if ((type(getattr(config, 'L', None)) != float and\
		type(getattr(config, 'L', None)) != int) or\
		getattr(config, 'L', None) <= 0.0):
	try:
		config.L = float(input("Length of domain not properly defined. Enter floating-point value: "))

		if config.L <= 0:
			raise NameError

		print "L set to " + str(config.L)

	except (NameError, SyntaxError, EOFError, ValueError):
		config.L = 4.0
		print "L defaulted to 4.0"

if type(config.L) == int:
	config.L = float(config.L)

# Choice of the diffusivity model
if (getattr(config, 'D_choice', "").lower() not in diffusivity_models or\
		type(getattr(config, 'D_choice', None)) != str):
	try:
		config.D_choice = raw_input("The diffusivity model is not properly chosen. Choose from the following: Zohm, Staps, Shear -> ")

		if config.D_choice.lower() not in diffusivity_models:
			raise IndexError()

	except (IndexError, EOFError):
		config.D_choice = "d_staps"
		print "Diffusivity model defaulted to Staps'."


# Initial starting mode
if type(getattr(config, 'initial_H_mode', None)) != bool:
	config.initial_H_mode = False
	print "Defaulted to starting in L--mode."


# Z-equation model choice
if type(getattr(config, 'original_model', None)) != bool:
	config.original_Z_model = True
	print "Defaulted to using the original numerical model for Z."


# Use SI units, or arbitrary?
if type(getattr(config, 'SI_units', False)) != bool:
	config.SI_units = False
	print "The units of the profiles are arbitrary."


# Total number of time steps
if ((type(getattr(config, 'total_timeSteps', None)) != int and\
		type(getattr(config, 'total_timeSteps', None)) != float) or\
		getattr(config, 'total_timeSteps', None) <= 0):
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
if ((type(getattr(config, 'timeStep_denom', None)) != float and\
		type(getattr(config, 'timeStep_denom', None)) != int) or\
		getattr(config, 'timeStep_denom', None) <= 0.0):
	try:
		config.timeStep_denom = float(input("The denomintor of the time step size is not properly defined. Enter floating-point value: "))

		if config.timeStep_denom <= 0:
			raise ValueError

		print "The denominator of the time step size is set to "\
				+ str(config.timeStep_denom)

	except (NameError, SyntaxError, EOFError, ValueError):
		config.timeStep_denom = 15.0
		print "The denominator of the time step size is defaulted to 15.0"

if type(config.timeStep_denom) == int:
	config.timeStep_denom = float(config.timeStep_denom)


# Set auxiliary plots
if type(getattr(config, 'aux_plots', None)) != bool:
	config.aux_plots = False

# Plot titles
if not hasattr(config, 'plot_title'):
	config.plot_title = ""
if not hasattr(config, 'aux_title1'):
	config.aux_title1 = ""
if not hasattr(config, 'aux_title2'):
	config.aux_title2 = ""

# Maximum x on the plots
if (type(getattr(config, 'plotx_max', None)) != int and\
		type(getattr(config, 'plotx_max', None)) != float):
	config.plotx_max = config.L

# MINIMUM y values on the plots
if (type(getattr(config, 'ploty_min', None)) != int and\
		type(getattr(config, 'ploty_min', None)) != float):
	config.ploty_min = None
if (type(getattr(config, 'aux1y_min', None)) != int and\
		type(getattr(config, 'aux1y_min', None)) != float):
	config.aux1y_min = None
if (type(getattr(config, 'aux2y_min', None)) != int and\
		type(getattr(config, 'aux2y_min', None)) != float):
	config.aux2y_min = None
# MAXIMUM y values on the plots
if (type(getattr(config, 'ploty_max', None)) != int and\
		type(getattr(config, 'ploty_max', None)) != float):
	config.ploty_max = None
if (type(getattr(config, 'aux1y_max', None)) != int and\
		type(getattr(config, 'aux1y_max', None)) != float):
	config.aux1y_max = None
if (type(getattr(config, 'aux2y_max', None)) != int and\
		type(getattr(config, 'aux2y_max', None)) != float):
	config.aux2y_max = None


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

