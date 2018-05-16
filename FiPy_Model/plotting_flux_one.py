import matplotlib.pyplot as plt
import numpy
from matplotlib import ticker

import os
import sys

# Import DIRECTORY NAME in which to do the bidding
data_file = sys.argv[1]


# Function to remove first and last y-axis tick labels
def remove_yticks(*the_axes):
	for axis in the_axes:
		plt.setp(axis.get_yticklabels()[0], visible=False)
		plt.setp(axis.get_yticklabels()[-1], visible=False)


data = numpy.loadtxt(data_file, unpack=True, skiprows=1)

# The data
x, density, temperature, Z, Diffusivity, D_an, Gamma_an, n_0, cx_rate,\
		Gamma_cx, D_bulk, Gamma_bulk, Gamma_ol = data[:]

# Generate the figure
fig_flux = plt.figure()
fig_flux = plt.subplots(4, 1, sharex=True, squeeze=True, figsize=(12,20))[0]
fig_flux.subplots_adjust(hspace=0)

ax_list = []

## TOP FLUX PLOT (Gamma_an)
ax_list.append(plt.subplot(4,1,1))
Gamma_an_plot = ax_list[0].plot(x, Gamma_an, label=r"$\Gamma_e^{an}$",\
		color='darkcyan', linewidth=2)
ax_list[0].set_ylabel(r"$\Gamma_e^{an}$", fontsize='xx-large',\
		rotation=0, labelpad=20)
ax_list[0].tick_params(axis='x', labelbottom='off')
ax_list[0].yaxis.set_major_locator(ticker.MaxNLocator(4))
ax_list[0].tick_params(axis='y', labelsize=24)

ax_list[0].grid(True)
#ax_list[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))

#plt.xlim(0.0,0.006)


## SECOND FLUX PLOT (Gamma_cx)
ax_list.append(plt.subplot(4,1,2))
Gamma_cx_plot = ax_list[1].plot(x, Gamma_cx, label=r"$\Gamma_i^{cx}$",\
		color='darkcyan', linewidth=2)
ax_list[1].set_ylabel(r"$\Gamma_i^{cx}$", fontsize='xx-large',\
		rotation=0, labelpad=20)
ax_list[1].tick_params(axis='x', labelbottom='off')
ax_list[1].tick_params(axis='y', labelsize=24)

ax_list[1].grid(True)
#ax_list[1].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax_list[1].yaxis.set_major_locator(ticker.MaxNLocator(4))

#plt.xlim(0.0,0.006)


## THIRD FLUX PLOT (Gamma_bulk and D_bulk)
ax_list.append(plt.subplot(4,1,3))
Gamma_bulk_plot = ax_list[2].plot(x, Gamma_bulk,\
		label=r"$\Gamma_i^{\pi\parallel}$", color='darkcyan', linewidth=2)
ax_list[2].set_ylabel(r"$\Gamma_i^{\pi\parallel}$", fontsize='xx-large',\
		rotation=0, labelpad=20)
ax_list[2].tick_params(axis='x', labelbottom='off')
ax_list[2].tick_params(axis='y', labelsize=24)

ax_list[2].grid(True)
#ax_list[2].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax_list[2].yaxis.set_major_locator(ticker.MaxNLocator(4))

#plt.xlim(0.0,0.006)


## FOURTH FLUX PLOT (Gamma_ol)
ax_list.append(plt.subplot(4,1,4))
Gamma_ol_plot = ax_list[3].plot(x, Gamma_ol, label=r"$\Gamma_i^{ol}$",\
		color='darkcyan', linewidth=2)
ax_list[3].set_ylabel(r"$\Gamma_i^{ol}$", fontsize='xx-large',\
		rotation=0, labelpad=20)
ax_list[3].tick_params(axis='y', labelsize=24)

ax_list[3].grid(True)
#ax_list[3].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax_list[3].yaxis.set_major_locator(ticker.MaxNLocator(4))

# The x-axis label and ticks
ax_list[3].grid(True)
ax_list[3].set_xlabel(r"$x$", fontsize='xx-large')

remove_yticks(ax_list[0], ax_list[1], ax_list[2], ax_list[3])

for k in range(len(ax_list)):
	ax_list[k].yaxis.tick_right()
	ax_list[k].tick_params(axis='y', labelsize='medium')

#plt.xlim(0.0,0.006)

fig_flux.suptitle(r"Fluxes, $t \,=\, 1600 \,\cdot\, 0.5 \, \mu$s", fontsize=24)
fig_flux.tight_layout(pad=0.05, w_pad=0.0)
#fig_flux.subplots_adjust(top=0.95)

#fig_flux.savefig(data_directory +'/Gamma'+ filename_sans_ext +'.png')

plt.show()

