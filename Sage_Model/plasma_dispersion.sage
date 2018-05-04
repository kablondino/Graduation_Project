reset()

var('x,z')

plasma_disp(z) = I*sqrt(pi)*exp(-z^2) * erfc(-I*z)

complex_one(x) = imag(plasma_disp.subs(z=x + 1*i))
complex_onetenth(x) = imag(plasma_disp.subs(z=x + 0.1*i))

complex_one_plot = plot(complex_one, (x,-4,4), color='olive', legend_label=r"Im$[X(Z + 1i)]$")
complex_onetenth_plot = plot(complex_onetenth, (x,-4,4), color='turquoise', legend_label=r"Im$[X(Z + 0.1i)]$")

non_complex_plot = plot(sqrt(pi)*exp(-(x)^2), (x,-4,4), color='blue', legend_label=r"$\sqrt{\pi} \, \exp(-Z^2)$", axes_labels=['$Z$', ''], gridlines=True)

combined_plots = complex_one_plot + complex_onetenth_plot + non_complex_plot
combined_plots.show()

