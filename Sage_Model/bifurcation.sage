reset()

var('x,t,a')

# ----------------- Co-dimension 1 ------------------------
a_min, a_max = -1.5, 0.5
x_min, x_max = -1.25, 1.25

# Co-dimension 1 equation
f1 = diff(x,t) == a + x^2

co_1 = plot_vector_field((0, f1 / sqrt(a^2 + x^2)), (a, a_min, a_max), (x, x_min, x_max), color='blue', axes=True, frame=True)


co_1_root = solve(a + x^2 == 0, x)
co_1_lower = implicit_plot(co_1_root[0], (a,a_min,a_max), (x,x_min,x_max), linewidth=2.5, axes_labels=['$a$','$x$'], color='red', plot_points=2000, gridlines=False, fontsize=18, title="Co-dimension 1 Fold")
co_1_upper = implicit_plot(co_1_root[1], (a,a_min,a_max), (x,x_min,x_max), linewidth=2.5, color='green', plot_points=2000)

co_1_zero = point((0,0), size=30, color='black')

co_1_total = co_1_lower + co_1_upper + co_1 + co_1_zero

#show(co_1_total, figsize=[8,10], typeset='type1')


# ----------------- Co-dimension 2 ------------------------
a_min, a_max = -1.0, 1.0
x_min, x_max = -1.5, 1.5
var('b')

# Co-dimension 2 equation
f2 = diff(x,t) == (a + b*x - x^3).subs(b=1)

co_2 = plot_vector_field((0, f2 / sqrt(a^2 + (x - x^3)^2)), (a, a_min, a_max), (x, x_min, x_max), color='blue', axes=True, frame=True, fontsize=16, title="Co-dimension 2 Fold")

turning_point_a = 1/3*sqrt(3)
turning_point_x = solve(f2.subs(a=turning_point_a), x)[0].right().real()
co_2_low = implicit_plot(f2.right() == 0, (a, a_min, a_max), (x, x_min, -turning_point_a), linewidth=2.5, axes_labels=['$a$','$x$'], aspect_ratio='automatic')
co_2_mid = implicit_plot(f2.right() == 0, (a, a_min, a_max), (x, -turning_point_a, turning_point_a), linewidth=2.5, color='red', aspect_ratio='automatic')
co_2_high = implicit_plot(f2.right() == 0, (a, a_min, a_max), (x, turning_point_a, x_max), linewidth=2.5, color='green', aspect_ratio='automatic')

#co_2_zero1 = point((turning_point_a, turning_point_x), size=40, color='black')
#co_2_zero2 = point((-turning_point_a, -turning_point_x), size=40, color='black')

co_2_total = co_2 + co_2_low + co_2_mid + co_2_high

show(co_2_total, figsize=[8,10], typeset='type1')

