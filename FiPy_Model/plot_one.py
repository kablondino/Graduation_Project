import matplotlib.pyplot as plt
import numpy

plt.rcParams.update({'font.size': 24})

data = numpy.loadtxt('./Original_Flow_Shear/0135.tsv', unpack=True, skiprows=1)

x, density, temperature, Z, Diffusivity = data[:]

density_plot = plt.plot(x, density, label="$n$", linewidth=2)
temperature_plot = plt.plot(x, temperature, label="$T$", linewidth=2)
Z_plot = plt.plot(x, -Z, label="$-Z$", linewidth=2)
Diffusivity_plot = plt.plot(x, Diffusivity, label="$D$", linewidth=2)

plt.ylim(-0.2, 5.2)

#the_title = r"$\Gamma_c = -0.8, \, D \sim 1 / [1 + 0.1 (Z)^2 + 0.5 (Z^\prime)^{2}],$"\
#		+"\n" + r"$t \,=\, 135, \, \Delta t \,=\, 1 / 375$"
the_title = r"$t \,=\, 135$"

plt.title(the_title, fontsize='xx-large')
plt.xlabel(r"$x$", fontsize='xx-large')

plt.legend(loc='best')

plt.grid(True)

plt.show()

