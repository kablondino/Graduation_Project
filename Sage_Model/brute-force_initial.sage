"""
	Calculates the fits for the initial density, temperature,
	and Z profiles. It includes both L-- and H--mode fits.
"""

reset()

# Density in L--Mode
AU_density(x) = 1.5*x / 4.0 + 0.5

# Density in H--mode
AU_density1 = (0.5/1.5)*x + 0.5
AU_density2 = 3.0*x - 3.5
AU_density3 = (0.5/1.5)*x + (11.0/6.0)*1.0

AU_density_piece = piecewise([([0,1.5], AU_density1),\
		((1.5,2.0), AU_density2),\
		([2.0,4.0], AU_density3)])

# Plot density profiles
show(plot(AU_density, (x,0,4), color='blue', legend_label="L--Mode")\
		+ plot(AU_density_piece, (x,0,4), color='magenta',\
		legend_label="H--mode: Linear Piecewise") ,axes_labels=["$x$", ""],\
		title="Initial Density Profiles in AU")

# Temperature in L--mode
var('x, p1,p2,p3')
L_mode_parabola(x) = p1 + p2*x + p3*x^2

AU_temperature_L_data = [(0.0,1.2),(2.0,1.65),(4.0,2.2)]
AU_temperature_L_parabola = L_mode_parabola.subs(find_fit(AU_temperature_L_data,\
		L_mode_parabola, solution_dict=True))
AU_L_plot = plot(AU_temperature_L_parabola, (x,0,4.0), color='green',\
		legend_label="L--Mode")

SI_temperature_L_data = [(0.0,100),(2.0,220),(4.0,400)]
SI_L_plot = plot(L_mode_parabola.subs(find_fit(SI_temperature_L_data,\
		L_mode_parabola, solution_dict=True)), (x,0,4.0), color='green',\
		legend_label="L--Mode")


# Temperature in H--mode
var('a,b')
line_model(x) = a + b*x

SI_Tpiece1(x) = line_model.subs(find_fit([(0,100),(1.5,180)], line_model, solution_dict=True))
SI_Tpiece2(x) = line_model.subs(find_fit([(1.5,180),(2,480)], line_model, solution_dict=True))
SI_Tpiece3(x) = line_model.subs(find_fit([(2,480),(4,600)], line_model, solution_dict=True))

SI_H_mode_piece = piecewise([([0,1.5], SI_Tpiece1),\
		((1.5,2),SI_Tpiece2),\
		([2,4],SI_Tpiece3)])

SI_H_piece_plot = plot(SI_H_mode_piece, (x,0,4), color='blue',\
		legend_label="H--Mode: Linear Piecewise")


AU_Tpiece1(x) = line_model.subs(find_fit([(0,1.2),(1.5,1.5)], line_model, solution_dict=True))
AU_Tpiece2(x) = line_model.subs(find_fit([(1.5,1.5),(2,2.4)], line_model, solution_dict=True))
AU_Tpiece3(x) = AU_Tpiece1(x) + AU_Tpiece2(2.0) - AU_Tpiece1(2.0)
AU_H_mode_piece = piecewise([([0,1.5], AU_Tpiece1),\
		((1.5,2),AU_Tpiece2),\
		([2,4],AU_Tpiece3)])
AU_H_piece_plot = plot(AU_H_mode_piece, (x,0,4), color='blue',\
		legend_label="H--Mode: Linear Piecewise")

# Plot temperature profiles in SI
show(list_plot(SI_temperature_L_data, color='green')\
#		+ list_plot(SI_temperature_H_data, color='red') + fitted_H_plot\
		+ SI_L_plot + SI_H_piece_plot,\
		axes_labels=["$x$", "eV"], title="Initial Temperature Profiles in SI")

# Plot temperature profiles in AU
show(list_plot(AU_temperature_L_data, color='green')\
		+ AU_L_plot + AU_H_piece_plot,\
		axes_labels=["$x$", ""], title="Initial Temperature Profiles in AU")

# Logistic fit for Z profile
var('core, k, x0, a,b,c,d')
H_mode_logistic(x) = core / (1 + exp(-k*(x - x0)))

Z_H_logistic(x) = H_mode_logistic.subs(core=3.0, k=-12.0, x0=1.75)
Z_H_logistic_plot = plot(Z_H_logistic, (x,0,4.0),\
		color='red', legend_label="H--Mode: Logistic Function")

## Z profile
# Linear, piecewise. It does not work in FiPy
Z1(x) = 3.0
Z2(x) = line_model.subs(find_fit([(1.5,3.0),(2,0.0)], line_model, solution_dict=True))
Z3(x) = 0.0

Z_piece = piecewise([([0,1.5], Z1),\
		((1.5,2),Z2),\
		([2,4],Z3)])
Z_piece_plot = plot(Z_piece, (x,0,4), color='orange',\
		legend_label="Z H--Mode: Linear Piecewise")

show(Z_H_logistic_plot + Z_piece_plot, title="Initial $Z$ Profile in H--Mode")

