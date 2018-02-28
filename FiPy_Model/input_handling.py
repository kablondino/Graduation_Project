"""
	This file deals with all of the inputs for the system. The current set
	of inputs are the following:
	nx:				int		The number of grid points
	L:				float	The length of the domain
	total_timeSteps:int		The total number of time steps
	timeStep_denom:	float	The denominator of the time step
	numerical_parameter:	string		The set of parameters
	D_choice:		string	The model of the diffusivity
	initial_H_mode:	bool	Start in L-- or H--mode?
	res_tol:		float	The residual tolerance while solving
	plot_title:		string	The title of the plot; can be formatted
	plot_max:		float	The maximum x-value on the plot

	save_directory: string	The name of the saving directory, from current\
							directory being run.
	save_plots:		bool	Should the plots be saved?
	save_TSVs:		bool	Should TSV files be generated and saves?

	Each possible input also has a default value, if nothing is set.
"""
import sys

# Import variables from job configuration file
config = __import__ (sys.argv[1].replace('.py',''))

parameter_sets = ["staps", "paquay", "g_grad", "gradient_model"]
diffusivity_models = ["d_zohm", "d_staps", "d_shear", "d_flow_shear"]


# -------------- Check Configuration Variables ------------
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

print config.nx

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


# Maximum allowable residual when sweeping
if ((type(getattr(config, 'res_tol', None)) != float and\
		type(getattr(config, 'res_tol', None)) != int) or\
		getattr(config, 'res_tol', None) <= 0):
	try:
		config.res_tol = float(input("The maximum allowable tolerance for the sweep residual is improperly set. Enter a floating-point value: "))

		if config.res_tol <= 0:
			raise ValueError

		print "The maximum residual tolerance is set to " +str(config.res_tol)

	except (NameError, SyntaxError, EOFError, ValueError):
		config.res_tol = 1.0e-6
		print "The maximum residual tolerance is defaulted to 1.0e-6."


# Plot title
if not hasattr(config, 'plot_title'):
	config.plot_title = ""

# Maximum x on the plot
if (type(getattr(config, 'plot_max', None)) != int and\
		type(getattr(config, 'plot_max', None)) != float):
	config.plot_max = config.L


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
		type(getattr(config, 'save_plots', None)) != bool) :
	config.save_plots = False
if (not hasattr(config, 'save_TSVs') or\
		type(getattr(config, 'save_TSVs', None)) != bool) :
	config.save_TSVs = False

