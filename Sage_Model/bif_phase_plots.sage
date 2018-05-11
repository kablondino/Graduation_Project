reset()

import numpy

var('x, a, b, x_dot,t')
x_min, x_max = -2.0, 2.0

"""
	ROOT FINDING function, for all of the roots.
	Returns a sorted list of UNIQUE roots. This does NOT work when there
	is multiplicity in the roots. Use `w.roots()` to symbolically check
	for multiplicity.
"""
def find_root_recursive(func, lower_bound, upper_bound, tol=1.0e-12):
	L = []
	try:
		x0 = find_root(func, lower_bound, upper_bound)
		L.append(x0)
		L += find_root_recursive(func, lower_bound, x0-tol, tol)
		L += find_root_recursive(func, x0+tol, upper_bound, tol)
	except:
		pass
	return L

# ----------------- \dot{x} vs x plots --------------------
x_dot = -(a + b*x - x^3)

a_list = [-2.0, 0.0, 2.0]
b_list = [-2.5, -0.5, 1.0, 3.0]


plot_list = []
point_list = []

a_selection = -0.5
b_selection = -0.8

# Create list of randomized integers from 0 to the length of a_list or b_list
random_a_index = [k for k in range(len(a_list))]
shuffle(random_a_index)
#random_b_index = [j for j in range(len(b_list))]
#shuffle(random_b_index)

# Options for plotting
the_title = '$\dot{x} = -(a + bx - x^3)$'
ax_labels = ['$x$', r'$\dot{x}$']
the_font_size = 54

for i in range(len(a_list)):
	this_loops_color = rainbow(len(a_list)+2)[i+2]

	the_label = '$a = ' + str(a_list[i].n(prec=10)) + '$'

	plot_list.append(plot(x_dot.subs(a=a_list[i], b=b_selection),\
			(x, x_min, x_max), color=this_loops_color, gridlines='major',\
			thickness = 4.0, title=the_title, frame=True, axes=False,\
			legend_label=the_label, axes_labels=ax_labels, figsize=16,\
			fontsize=the_font_size, typeset='type1'))

	# Create list of points of zeros
	particular_roots = find_root_recursive(x_dot.subs(a=a_list[i], b=b_selection), x_min, x_max)
#	particular_roots = [find_root(x_dot.subs(a=a_list[i], b=b_selection), x_min, x_max)]
	for j in particular_roots:
		point_list.append(ellipse((j,0), 0.03, 0.2, color=this_loops_color,\
				thickness=4.0, aspect_ratio='automatic'))


# Create the parameter box for b
b_parameter_box1 = text('$b = ' + str(b_selection.n(prec=10)) + '$',\
		(x_max-0.1, -10.0), bounding_box={'boxstyle':'round', 'fc':'w'},\
		fontsize=the_font_size, color='black', horizontal_alignment='right')


# Tangent root function
#tangent_root_fun = x_dot.subs(a=sqrt(4*b_selection^3/27), b=b_selection)

# Equation with exactly 2 ROOTS
#plot_list.append(plot(tangent_root_fun, (x, x_min, x_max), color='black',\
#		legend_label=r'$a = \sqrt{\frac{4\,b^3}{27}}$', thickness=5.0))

## Find lower, tangent root
#point_list.append(ellipse((find_root(tangent_root_fun, x_min, 0),0), 0.03,\
#		0.2, color='black', thickness=5.0, aspect_ratio='automatic'))

# Find upper, non-tangent root, and make point
#point_list.append(ellipse((find_root_recursive(tangent_root_fun, 0, x_max)[0],\
#		0), 0.03, 0.2, color='black', thickness=5.0, aspect_ratio='automatic'))

# Set legend options
combined_plots1 = sum(plot_list) + sum(point_list)
combined_plots1.set_legend_options(font_size=the_font_size-14)
combined_with_box1 = combined_plots1 + b_parameter_box1

show(combined_with_box1)

# ----------------- Co-dimension 2 Surface ----------------
#f3 = diff(x,t) == -(a + b*x - x^3)
#a_min, a_max = -2, 2
#co_2_surface = implicit_plot3d(f3, (a,a_min,a_max), (b,-1.0,2.0),\
#		(x,x_min,x_max), color=(x,colormaps.gist_rainbow))
#
#co_2_surface.show()

