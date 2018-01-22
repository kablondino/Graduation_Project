reset()

import numpy

var('x,t,a,b')

# ----------------- Co-dimension 1 ------------------------
#a_min, a_max = -1.5, 0.5
#x_min, x_max = -1.25, 1.25

## Co-dimension 1 equation
#f1 = diff(x,t) == a + x^2
#
#co_1 = plot_vector_field((0, f1 / sqrt(a^2 + x^2)), (a, a_min, a_max), (x, x_min, x_max), color='blue', axes=True, frame=True)
#
#
#co_1_root = solve(a + x^2 == 0, x)
#co_1_lower = implicit_plot(co_1_root[0], (a,a_min,a_max), (x,x_min,x_max), linewidth=2.5, axes_labels=['$a$','$x$'], color='red', plot_points=2000, gridlines=False, fontsize=20, title="Co-dimension 1 Fold")
#co_1_upper = implicit_plot(co_1_root[1], (a,a_min,a_max), (x,x_min,x_max), linewidth=2.5, color='green', plot_points=2000, figsize=10, typeset='type1')
#
## Zero point
#co_1_zero = point((0,0), size=30, color='black')
#
#co_1_total = co_1_lower + co_1_upper + co_1 + co_1_zero
#
##co_1_total.show()
#
#
## ----------------- Co-dimension 2 ------------------------
#a_min, a_max = -1.2, 1.2
#x_min, x_max = -1.5, 1.5
#var('b')
#
## Co-dimension 2 equation
#f2 = diff(x,t) == -(a + b*x - x^3).subs(b=1.1)
#
#co_2 = plot_vector_field((0, f2), (a, a_min, a_max), (x, x_min, x_max), color='blue', axes=True, frame=True)
#
## Points
#turning_point_a = 1/3*sqrt(3)
#turning_point_x = solve(f2.subs(a=turning_point_a), x)[0].right().real()
#
#co_2_low = implicit_plot(f2.right() == 0, (a, a_min, a_max), (x, x_min, -turning_point_a), linewidth=2.5, axes_labels=['$a$','$x$'], axes=True, fontsize=20, title="Co-dimension 2 Fold")
#co_2_mid = implicit_plot(f2.right() == 0, (a, a_min, a_max), (x, -turning_point_a, turning_point_a), linewidth=2.5, color='red', axes=True)
#co_2_high = implicit_plot(f2.right() == 0, (a, a_min, a_max), (x, turning_point_a, x_max), linewidth=2.5, color='green', axes=True, figsize=10, typeset='type1')
#
##co_2_zero1 = point((turning_point_a, turning_point_x), size=40, color='black')
##co_2_zero2 = point((-turning_point_a, -turning_point_x), size=40, color='black')
#
#co_2_total = co_2 + co_2_low + co_2_mid + co_2_high

#co_2_total.show()

# ----------------- \dot{x} vs x plots --------------------
var('x_dot')
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

x_dot = -(a + b*x - x^3)

a_list = [-2.0, 0.0, 2.0]
b_list = [-2.5, -0.5, 1.0, 3.0]


plot_list = []
point_list = []

a_selection = -0.5
b_selection = 1.1

# Create list of randomized integers from 0 to the length of a_list or b_list
random_a_index = [k for k in range(len(a_list))]
shuffle(random_a_index)
#random_b_index = [j for j in range(len(b_list))]
#shuffle(random_b_index)

# Options for plotting
the_title = '$\dot{x} = -(a + bx - x^3)$'
ax_labels = ['$x$', r'$\frac{{d}x}{{d}t}$']
the_font_size = 54

for i in range(len(a_list)):
	this_loops_color = rainbow(len(a_list)+2)[i+2]

	the_label = '$a = ' + str(a_list[i].n(prec=10)) + '$'

	plot_list.append(plot(x_dot.subs(a=a_list[i], b=b_selection),\
			(x, x_min, x_max), color=this_loops_color, gridlines='major',\
			thickness = 5.0, title=the_title, frame=True, axes=False,\
			legend_label=the_label, axes_labels=ax_labels, figsize=16,\
			fontsize=the_font_size, typeset='type1'))

	# Create list of points of zeros
	particular_roots = find_root_recursive(x_dot.subs(a=a_list[i], b=b_selection), x_min, x_max)
#	particular_roots = [find_root(x_dot.subs(a=a_list[i], b=b_selection), x_min, x_max)]
	for j in particular_roots:
		point_list.append(ellipse((j,0), 0.03, 0.2, color=this_loops_color,\
				thickness=5.0, aspect_ratio='automatic'))


# Create the parameter box for b
b_parameter_box = text('$b = ' + str(b_selection.n(prec=10)) + '$',\
		(x_max-0.1, -6.4), bounding_box={'boxstyle':'round', 'fc':'w'},\
		fontsize=the_font_size, color='black', horizontal_alignment='right')


# Tangent root function
tangent_root_fun = x_dot.subs(a=sqrt(4*b_selection^3/27), b=b_selection)

# Equation with exactly 2 ROOTS
plot_list.append(plot(tangent_root_fun, (x, x_min, x_max), color='black',\
		legend_label=r'$a = \sqrt{\frac{4\,b^3}{27}}$', thickness=5.0))

## Find lower, tangent root
point_list.append(ellipse((find_root(tangent_root_fun, x_min, 0),0), 0.03,\
		0.2, color='black', thickness=5.0, aspect_ratio='automatic'))

# Find upper, non-tangent root, and make point
point_list.append(ellipse((find_root_recursive(tangent_root_fun, 0, x_max)[0],\
		0), 0.03, 0.2, color='black', thickness=5.0, aspect_ratio='automatic'))

# Set legend options
combined_plots = sum(plot_list) + sum(point_list)
combined_plots.set_legend_options(font_size=the_font_size-14)
combined_with_box = combined_plots + b_parameter_box

show(combined_with_box)

# ----------------- Co-dimension 2 Surface ----------------
#f3 = diff(x,t) == -(a + b*x - x^3)
#co_2_surface = implicit_plot3d(f3, (a,a_min,a_max), (b,-1.0,2.0),\
#		(x,x_min,x_max), color=(x,colormaps.gist_rainbow))

#co_2_surface.show()

