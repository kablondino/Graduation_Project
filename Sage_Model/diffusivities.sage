reset()

var('x')

L = 5.0
Z_S = 3/2
D_min = 2/5
D_max = 5

Z = Z_S*(1 - tanh((L*x - L) / 2))

D_Zohm = (D_min + D_max) / 2 + ((D_max - D_min)*tanh(Z)) / 2
D_Staps = D_min + (D_max - D_min) / (1 + (diff(Z,x))^2)
D_Biglari = D_min / (1 + (diff(Z,x))^(2/3))

a1, a2, a3 = 1.0, 0.0, 1.5	# ASSUMES a2 = 0
D_Shear = D_min + (D_max - D_min) / (1.0 + a1*(Z)**2 + a2*Z*(diff(Z,x)) + a3*(diff(Z,x)^2))

the_plots = []
the_colors = rainbow(5)

the_plots.append(plot(D_Zohm, (x,0,L-1.5), color=the_colors[0], legend_label=r"$\sim \tanh(Z)$"))
the_plots.append(plot(D_Staps, (x,0,L-1.5), color=the_colors[1], legend_label=r"$\sim \left[1 + (Z^{\prime})^2\right]^{-1}$"))
the_plots.append(plot(D_Wilson1, (x,0,L-1.5), color=the_colors[2], legend_label=r"$\sim \left[1 + (Z^{\prime})^{2/3}}\right]^{-1}"))
the_plots.append(plot(D_Wilson2, (x,0,L-1.5), color=the_colors[3], legend_label=r"$\sim \left[1 + (Z^{\prime})^2\right]^{-3/2}$"))
the_plots.append(plot(D_Shear, (x,0,L-1.5), color=the_colors[4], legend_label=r"$\sim \left[1 + a_1 Z^2 + a_2 (Z^{\prime})^2\right]^{-1}$", gridlines=True, axes_labels=['$x$','$D$']))

combined_plots = sum(the_plots)

show(sum(the_plots))

