reset()

import numpy

"""
	This file attempts to calculate the critical b value in
	the co-dimension 2 bifurcation such that a jump occurs,
	but only barely.
"""

var('x,t,a,b')
# ----------------- Co-dimension 2 ------------------------
a_min, a_max = -1.0, 1.0
x_min, x_max = -1.4, 1.4
var('b')

# Co-dimension 2 equation
f2 = diff(x,t) == -(a + b*x - x^3)

# Points
turning_point_a = 1/3*sqrt(3)
turning_point_x = solve(f2.subs(a=turning_point_a), x)[0].right().real()

the_axes_labels = ['$a$','$x$']
the_figsize = 10
the_title = "Co-dimension 2 Fold"
the_linewidth = 1.0
the_fontsize = 10

b_values = numpy.arange(-0.4, 0.6, 0.2)

the_colors = rainbow(len(b_values))

co_2_loop = []
b_text = []

for i in range(len(b_values)):
	co_2_loop.append(implicit_plot(f2.right().subs(b=b_values[i]) == 0,\
			(a, a_min, a_max), (x, x_min, x_max), linewidth=the_linewidth,\
			color=the_colors[i], axes_labels=the_axes_labels,\
			axes=True,\
			fontsize=the_fontsize, title=the_title))

	b_text.append(text('$b = ' + str(b_values[i]) + '$\n', (0.5, -0.1 - i/10.0),\
			color=the_colors[i]))


the_plots = sum(co_2_loop) + sum(b_text)

show(the_plots)

## OLD splitting of function
#co_2_low = implicit_plot(f2.right() == 0, (a, a_min, a_max),\
#		(x, x_min, -turning_point_a), linewidth=the_linewidth,\
#		axes_labels=the_axes_labels, axes=True, fontsize=the_fontsize,\
#		title=the_title)
#co_2_mid = implicit_plot(f2.right() == 0, (a, a_min, a_max),\
#		(x, -turning_point_a, turning_point_a), linewidth=the_linewidth,\
#		color='red', axes=True)
#co_2_high = implicit_plot(f2.right() == 0, (a, a_min, a_max),\
#		(x, turning_point_a, x_max), linewidth=the_linewidth,\
#		color='green', axes=True, typeset='type1')

#co_2_zero1 = point((turning_point_a, turning_point_x), size=40, color='black')
#co_2_zero2 = point((-turning_point_a, -turning_point_x), size=40, color='black')


