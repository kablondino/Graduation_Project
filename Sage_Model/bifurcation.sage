reset()

var('x,t,a')

# ----------------- Co-dimension 1 ------------------------
a_min, a_max = -1.5, 0.5
x_min, x_max = -1.25, 1.25

# Co-dimension 1 equation
f1 = diff(x,t) == a + x^2

co_1 = plot_vector_field((0, f1 / sqrt(a^2 + x^2)), (a, a_min, a_max), (x, x_min, x_max), color='blue', axes=True, frame=True)


co_1_root = solve(a + x^2 == 0, x)
co_1_lower = implicit_plot(co_1_root[0], (a,a_min,a_max), (x,x_min,x_max), linewidth=2.5, axes_labels=['$a$','$x$'], color='red', plot_points=2000, gridlines=False, fontsize=20, title="Co-dimension 1 Fold")
co_1_upper = implicit_plot(co_1_root[1], (a,a_min,a_max), (x,x_min,x_max), linewidth=2.5, color='green', plot_points=2000, figsize=10, typeset='type1')

# Zero point
co_1_zero = point((0,0), size=30, color='black')

co_1_total = co_1_lower + co_1_upper + co_1 + co_1_zero

#co_1_total.show()


# ----------------- Co-dimension 2 ------------------------
a_min, a_max = -1.2, 1.2
x_min, x_max = -1.5, 1.5
var('b')

# Co-dimension 2 equation
f2 = diff(x,t) == -(a + b*x - x^3).subs(b=1.1)

co_2 = plot_vector_field((0, f2), (a, a_min, a_max), (x, x_min, x_max), color='blue', axes=True, frame=True)

# Points
turning_point_a = 1/3*sqrt(3)
turning_point_x = solve(f2.subs(a=turning_point_a), x)[0].right().real()

co_2_low = implicit_plot(f2.right() == 0, (a, a_min, a_max), (x, x_min, -turning_point_a), linewidth=2.5, axes_labels=['$a$','$x$'], axes=True, fontsize=20, title="Co-dimension 2 Fold")
co_2_mid = implicit_plot(f2.right() == 0, (a, a_min, a_max), (x, -turning_point_a, turning_point_a), linewidth=2.5, color='red', axes=True)
co_2_high = implicit_plot(f2.right() == 0, (a, a_min, a_max), (x, turning_point_a, x_max), linewidth=2.5, color='green', axes=True, figsize=10, typeset='type1')

#co_2_zero1 = point((turning_point_a, turning_point_x), size=40, color='black')
#co_2_zero2 = point((-turning_point_a, -turning_point_x), size=40, color='black')

co_2_total = co_2 + co_2_low + co_2_mid + co_2_high

#co_2_total.show()

import numpy
# \dot{x} vs x phase plots
var('x_dot')
co_2_dot = implicit_plot(x_dot == -(a + b*x - x^3).subs(a=-2.5, b=-2.5), (x, x_min, x_max), (x_dot, -4.0, 4.0), axes=True, aspect_ratio='automatic')
for i in numpy.arange(-2.0, 2.5, 0.5):
	for j in numpy.arange(-2.0, 2.5, 0.5):
		co_2_dot = co_2_dot + implicit_plot(x_dot == -(a + b*x - x^3).subs(a=i, b=j), (x, x_min, x_max), (x_dot, -4.0, 4.0), axes=True, aspect_ratio='automatic')

co_2_dot.show()

# ----------------- Co-dimension 2 Surface ----------------
#f3 = diff(x,t) == -(a + b*x - x^3)
#co_2_surface = implicit_plot3d(f3, (a,a_min,a_max), (b,-1.0,2.0), (x,x_min,x_max), color=(x,colormaps.gist_rainbow))

#co_2_surface.show()

