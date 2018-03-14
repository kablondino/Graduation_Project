reset()

plasma_disp(x) = I*sqrt(pi)*exp(-x^2) * erfc(-I*x)

the_title = r"$X(z) = i\sqrt{\pi} \, e^{-z^2} \, $erfc$(-i z) = 2i \, e^{-z^2} \int_-\infty^{ix} e^{-t^2} \, dt$"

real_plot = plot(plasma_disp(x).real(), (x,-5,5), color='red',\
		legend_label=r"Re$(X)$", axes_labels=[r"$z$", ""], typeset='type1',\
		title=the_title)
imag_plot = plot(plasma_disp(x).imag(), (x,-5,5), color='green',\
		legend_label=r"Im$(X)$")

combined_plot = real_plot + imag_plot
combined_plot.show()

