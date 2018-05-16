import matplotlib.pyplot as plt
import numpy
from matplotlib.ticker import FormatStrFormatter

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

# Generate the figures
fig_state = plt.figure()
fig_state = plt.subplots(2, 1, sharex=True, squeeze=True, figsize=(12,12))[0]
fig_state.subplots_adjust(hspace=0)

ax_list = []

# TOP STATE PLOT
ax_list.append(plt.subplot(2,1,1))
density_plot = ax_list[0].plot(x, density, label=r"$n$",\
		color='blue', linewidth=2)
ax_list[0].set_ylabel(r"$n$", fontsize='xx-large', rotation=0, labelpad=20)
ax_list[0].tick_params(axis='x', labelbottom='off') # Remove top plot's x-axis labels
ax_list[0].tick_params(axis='y', labelsize=24)

ax_list.append(ax_list[0].twinx())
#ax_list[0].set_ylim(0.0, 1.2e18)
#ax_list[1].set_ylim(140, 250)
temp_plot = ax_list[1].plot(x, temperature, label=r"$T$",\
		color='red', linewidth=2)
ax_list[1].set_ylabel(r"$T$", fontsize='xx-large', rotation=0, labelpad=20)
ax_list[1].tick_params(axis='y', labelsize=24)

ax_list[1].set_yticks(numpy.linspace(ax_list[1].get_yticks()[0],\
		ax_list[1].get_yticks()[-1], len(ax_list[0].get_yticks())))

ax_list[0].grid(True)
#ax_list[0].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
ax_list[1].grid(True)
ax_list[1].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))


#plt.xlim(0.0,0.006)


# BOTTOM STATE PLOT
ax_list.append(plt.subplot(2,1,2))
Z_plot = ax_list[2].plot(x, Z, label=r"$Z$", color='green', linewidth=2)
ax_list[2].set_ylabel(r"$Z$", fontsize='xx-large', rotation=0, labelpad=20)
ax_list[2].tick_params(axis='y', labelsize=24)
plt.ylim(-0.1, 5.1)

ax_list.append(ax_list[2].twinx())
D_plot = ax_list[3].plot(x, Diffusivity, label=r"$D$",\
		color='orange', linewidth=2)
ax_list[3].set_ylabel(r"$D$", fontsize='xx-large', rotation=0, labelpad=20)
ax_list[3].tick_params(axis='y', labelsize=24)
plt.ylim(-0.1, 5.1)

ax_list[2].set_yticks(numpy.linspace(ax_list[2].get_yticks()[0],\
		ax_list[2].get_yticks()[-1], len(ax_list[3].get_yticks())))

ax_list[2].set_ylim(-0.1, 5.1)
ax_list[3].set_ylim(-0.1, 5.1)

ax_list[2].grid(True)
ax_list[2].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_list[3].grid(True)
ax_list[3].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax_list[2].set_xlabel(r"$x$", fontsize='xx-large')

remove_yticks(ax_list[0], ax_list[1], ax_list[2], ax_list[3])

top_plot = density_plot + temp_plot
top_labels = [l.get_label() for l in top_plot]
ax_list[1].legend(top_plot, top_labels, loc='lower right', fontsize='xx-large')
bottom_plot = Z_plot + D_plot
bottom_labels = [l.get_label() for l in bottom_plot]
ax_list[3].legend(bottom_plot, bottom_labels, loc='lower right', fontsize='xx-large')

#plt.xlim(0.0,0.006)

#fig_state.suptitle(r"$D \sim 1 / [1 + 0.001 (Z)^2 + 0.0005 (Z\,^\prime)^2]$"+ "\n"\
#		+ r"$\Gamma_c = -1.0\times 10^{20}$, $t \,=\, 1060 \,\cdot\, 0.5 \, \mu$s",\
#		fontsize=20)
fig_state.suptitle(r"$\Gamma_c \,=\, -10^{21}$, $t \,=\, 190 \,\cdot\, 5 \, \mu$s",\
		fontsize=24)
fig_state.tight_layout(pad=0.2, w_pad=0.0)
plt.subplots_adjust(top=0.94)

plt.show()

