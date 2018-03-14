reset()

var('z')

plasma_disp(z) = I*sqrt(pi)*exp(-z^2) * erfc(-I*z)

the_title = r"$X(z) = i\sqrt{\pi} \, e^{-z^2} \, $erfc$(-i z) = 2i \, e^{-z^2} \int_{-\infty}^{ix} e^{-t^2} \, dt$"

real_plot = plot(plasma_disp(z).real(), (z,-5,5), color='red',\
		legend_label=r"Re$(X)$", axes_labels=[r"$z$", ""],\
		title=the_title)
imag_plot = plot(plasma_disp(z).imag(), (z,-5,5), color='green',\
		legend_label=r"Im$(X)$")

combined_plot = real_plot + imag_plot
combined_plot.show()

