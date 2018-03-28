"""
	Calculates the fits for the initial density, temperature,
	and Z profiles. It includes both L-- and H--mode fits.
"""

reset()

# Length of domain, in meters; for use in SI
L = 0.03

var('core, k, x0, a,b,c,d')
line_model(x) = a + b*x
H_mode_logistic(x) = core / (1 + exp(-k*(x - x0)))

## Density in L--Mode
AU_density(x) = 1.5*x / 4.0 + 0.5

## Density in H--mode
AU_density1 = (0.5/1.5)*x + 0.5
AU_density2 = 3.0*x - 3.5
AU_density3 = (0.5/1.5)*x + (11.0/6.0)

AU_density_piece = piecewise([([0,1.5], AU_density1),\
		((1.5,2.0), AU_density2),\
		([2.0,4.0], AU_density3)])

# Density in L--Mode for SI
SI_L_density(x) = ((2.0e19-0.5e19) / L)*x + 0.5e19

#show(plot(SI_density, (x,0,L), legend_label="L--Mode",\
#		axes_labels=["$x$ (m)", "m$^{-3}$"], title="Initial Density Profile"))

# Plot AU density profiles
show(plot(AU_density, (x,0,4.0), color='blue', legend_label="L--Mode")\
		+ plot(AU_density_piece, (x,0,4.0), color='magenta',\
		legend_label="H--mode: Linear Piecewise") ,axes_labels=["$x$", ""],\
		title="Initial Density Profiles in AU")

# Density in SI
SI_density1(x) = 0.375e21*x + 0.5e19
SI_density2(x) = 3.375e21*(x - 0.01) + SI_density1(0.01)
SI_density3(x) = 0.375e21*(x - 0.015) + SI_density2(0.015)
SI_H_density_piece = piecewise([([0.0,0.01], SI_density1),\
		((0.01,0.015), SI_density2), ([0.015,L], SI_density3)])

show(plot(SI_L_density, (x,0,L), color='blue', legend_label="L--Mode")\
		+ plot(SI_H_density_piece, (x,0,L), color='magenta',\
		legend_label="H--Mode", axes_labels=["$x$", "m$^{-3}$"],\
		title="Initial Density Profiles in SI"))

# Temperature in L--mode
var('x, p1,p2,p3,p4')
L_mode_parabola(x) = p1 + p2*x + p3*x^2
L_mode_sqroot(x) = p1 + sqrt(p2*x + p3)
L_mode_log(x) = p1*log(p2*x + p3)
L_mode_inverse(x) = p1 + (p2*x)^p3

#AU_temperature_L_data = [(0.0,1.2),(L/2,1.65),(L,2.2)]
#AU_temperature_L_parabola = L_mode_parabola.subs(find_fit(AU_temperature_L_data,\
#		L_mode_parabola, solution_dict=True))
#AU_L_plot = plot(AU_temperature_L_parabola, (x,0,L), color='green',\
#		legend_label="L--Mode")

SI_temperature_L_data = [(0.0,100),(L/2,300),(L,400)]
SI_temp_L_line_data = [(0.0,100),(L/2, 400-100),(L,400)]

SI_L_plot = plot(L_mode_log.subs(find_fit(SI_temperature_L_data,\
		L_mode_log, solution_dict=True)), (x,0,L), color='red',\
		legend_label="L--Mode")

SI_L_line_plot = plot(line_model.subs(find_fit(SI_temp_L_line_data,\
		line_model, solution_dict=True)), (x,0,L), color='green',\
		legend_label="L--Mode")


# Temperature in H--mode
#SI_Tpiece1(x) = line_model.subs(find_fit([(0,100),(0.01,180)], line_model, solution_dict=True))
#SI_Tpiece2(x) = line_model.subs(find_fit([(0.01,180),(0.015,480)], line_model, solution_dict=True))
#SI_Tpiece3(x) = line_model.subs(find_fit([(0.015,480),(L,600)], line_model, solution_dict=True))

SI_Tpiece1(x) = 8.0e3*x + 100
SI_Tpiece2(x) = 6.0e4*x - 420
SI_Tpiece3(x) = 8.0e3*x + 360

SI_H_mode_piece = piecewise([([0,0.01], SI_Tpiece1),\
		((0.01,0.015),SI_Tpiece2),\
		([0.015,L],SI_Tpiece3)])

SI_H_piece_plot = plot(SI_H_mode_piece, (x,0,L), color='blue',\
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
		+ list_plot(SI_temp_L_line_data, color='red')\
		+ SI_L_plot\
		+ SI_H_piece_plot
		+ SI_L_line_plot\
		,\
		axes_labels=["$x$", "eV"], title="Initial Temperature Profiles")

# Plot temperature profiles in AU
#show(list_plot(AU_temperature_L_data, color='green')\
#		+ AU_L_plot + AU_H_piece_plot,\
#		axes_labels=["$x$", ""], title="Initial Temperature Profiles in AU")

# Logistic fit for Z profile
Z_H_logistic(x) = H_mode_logistic.subs(core=3.0, k=-12.0, x0=1.75)
Z_H_logistic_plot = plot(Z_H_logistic, (x,0,L),\
		color='red', legend_label="H--Mode: Logistic Function")

## Z profile
# Linear, piecewise. It does not work in FiPy
#Z1(x) = 3.0
#Z2(x) = line_model.subs(find_fit([(1.5,3.0),(2,0.0)], line_model, solution_dict=True))
#Z3(x) = 0.0
#
#Z_piece = piecewise([([0,1.5], Z1),\
#		((1.5,2),Z2),\
#		([2,4],Z3)])
#Z_piece_plot = plot(Z_piece, (x,0,4), color='orange',\
#		legend_label="Z H--Mode: Linear Piecewise")
#
#show(Z_H_logistic_plot + Z_piece_plot, title="Initial $Z$ Profile in H--Mode")

