reset()

import numpy

var('x,t,a,b')
# ----------------- Co-dimension 2 ------------------------
a_min, a_max = -1.1, 1.1
x_min, x_max = -1.4, 1.4
var('b')

# Co-dimension 2 equation
f2 = diff(x,t) == -(a + b*x - x^3)

the_axes_labels = ['$a$','$x$']
the_figsize = 10
the_title = "Co-dimension 2 Fold"
the_linewidth = 2.0
the_fontsize = 20

b_values = [-0.5, 0.0, 1.0]
b_plots = []
the_label = []
the_colors = rainbow(len(b_values))

for i in range(len(b_values)):

	b_plots.append(implicit_plot(f2.right().real().subs(b=b_values[i]) == 0,\
			(a, a_min, a_max), (x, x_min, x_max), linewidth=the_linewidth,\
			color=the_colors[i], axes_labels=the_axes_labels,\
			axes=True, gridlines=True, figsize=10,\
			fontsize=the_fontsize, title=the_title, typeset='type1'))

boxes = [text('$b = -0.5$', (0.9, 0.7),\
		bounding_box={'boxstyle':'round', 'fc':'w'},\
		fontsize=the_fontsize-4, color='black')]
boxes.append(text('$b = 0.0$', (0.95, 1.0),\
		bounding_box={'boxstyle':'round', 'fc':'w'},\
		fontsize=the_fontsize-4, color='black'))
boxes.append(text('$b = 1.0$',(0.2, 1.2),\
		bounding_box={'boxstyle':'round', 'fc':'w'},\
		fontsize=the_fontsize-4, color='black'))
boxes.append(text('$a - bx + x^3 = 0$', (0.5, -1.2),\
		bounding_box={'boxstyle':'round', 'fc':'w'},\
		fontsize=the_fontsize, color='black'))

combined_plots = sum(b_plots) + sum(boxes)

show(combined_plots)

