import matplotlib.pyplot as plt
import numpy
from matplotlib.ticker import FormatStrFormatter

import os
import sys

# Import DIRECTORY NAME in which to do the bidding
data_directory = sys.argv[1]


# Create list of files
file_list = []

# Function to remove first and last y-axis tick labels
def remove_yticks(*the_axes):
	for axis in the_axes:
		plt.setp(axis.get_yticklabels()[0], visible=False)
		plt.setp(axis.get_yticklabels()[-1], visible=False)


# Make the file_list and ALL the data
for filename in os.listdir('./'+str(data_directory)):
	if filename.endswith(".tsv"):
		file_list.append(filename)

# Sort the file_list
file_list.sort(); big_data_list = []; ax_list = []

for filename in file_list:
	big_data_list.append(numpy.genfromtxt(data_directory+'/'+filename,\
			delimiter='\t', unpack=True, skip_header=1))

# Stack the arrays
stacked_data = numpy.dstack(tuple(big_data_list))
# Get the minima and maxima
the_mins = numpy.amin(numpy.amin(stacked_data, axis=2), axis=1)
the_maxs = numpy.amax(numpy.amax(stacked_data, axis=2), axis=1)

big_data_list = []
del big_data_list

i = 0 # Looping counter
# State plots
for filename in file_list:
	filename_sans_ext = os.path.splitext(filename)[0]
	data = numpy.loadtxt(data_directory+'/'+filename, unpack=True, skiprows=1)

	# The data
	x, density, temperature, Z, Diffusivity, D_an, Gamma_an, n_0, cx_rate,\
			Gamma_cx, D_bulk, Gamma_bulk, Gamma_ol = data[:]

	# Generate the figures
	fig_state = plt.figure()
	fig_state = plt.subplots(2, 1, sharex=True, squeeze=True, figsize=(12,12))[0]
	fig_state.subplots_adjust(hspace=0)


	# TOP STATE PLOT
	ax_list.append(plt.subplot(2,1,1))
	density_plot = ax_list[0].plot(x, density, label=r"$n$",\
			color='blue', linewidth=2)
	ax_list[0].set_ylabel(r"$n$", fontsize='large', rotation=0, labelpad=20)
	ax_list[0].tick_params(axis='x', labelbottom='off') # Remove top plot's x-axis labels
	plt.ylim((the_mins[1], the_maxs[1]))

	ax_list.append(ax_list[0].twinx())
	temp_plot = ax_list[1].plot(x, temperature, label=r"$T$",\
			color='red', linewidth=2)
	ax_list[1].set_ylabel(r"$T$", fontsize='large', rotation=0, labelpad=20)
	plt.ylim((the_mins[2], the_maxs[2]))

	ax_list[1].set_yticks(numpy.linspace(ax_list[1].get_yticks()[0],\
			ax_list[1].get_yticks()[-1], len(ax_list[0].get_yticks())))

	ax_list[0].grid(True)
	ax_list[0].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
	ax_list[1].grid(True)
	ax_list[1].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))


	# BOTTOM STATE PLOT
	ax_list.append(plt.subplot(2,1,2))
	Z_plot = ax_list[2].plot(x, Z, label=r"$Z$", color='green', linewidth=2)
	ax_list[2].set_ylabel(r"$Z$", fontsize='large', rotation=0, labelpad=20)
	plt.ylim((0.0, the_maxs[4]+0.1))

	ax_list.append(ax_list[2].twinx())
	D_plot = ax_list[3].plot(x, Diffusivity, label=r"$D$",\
			color='orange', linewidth=2)
	ax_list[3].set_ylabel(r"$D$", fontsize='large', rotation=0, labelpad=20)
	plt.ylim((0.0, the_maxs[4]+0.1))

#	ax_list[2].set_yticks(numpy.linspace(ax_list[2].get_yticks()[0],\
#			ax_list[2].get_yticks()[-1], len(ax_list[3].get_yticks())))

	ax_list[2].grid(True)
	ax_list[2].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
	ax_list[3].grid(True)
	ax_list[3].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

	ax_list[2].set_xlabel(r"$x$", fontsize='large')

	remove_yticks(ax_list[0], ax_list[1], ax_list[2], ax_list[3])

	top_plot = density_plot + temp_plot
	top_labels = [l.get_label() for l in top_plot]
	ax_list[1].legend(top_plot, top_labels, loc='best', fontsize='medium')
	bottom_plot = Z_plot + D_plot
	bottom_labels = [l.get_label() for l in bottom_plot]
	ax_list[3].legend(bottom_plot, bottom_labels, loc='best')


	fig_state.suptitle(r"$\Gamma_c = -6.0\times 10^{20}$, $D \sim 1 / [1 + 0.01 (Z)^2 + 0.001 (Z_x)^2]$"+ "\n" +"$\mu = 0.05$, $t = " +str(int(filename_sans_ext))+ "$",\
			fontsize=22)
	fig_state.tight_layout(pad=0.2, w_pad=0.0)
	plt.subplots_adjust(top=0.9)

	fig_state.savefig(data_directory +'/'+ filename_sans_ext +'.png')

#	print str(i) + "\t| Saved " + str(filename_sans_ext) + ".png"
#	fig_state.show()
#	raw_input("BREAK"); break
	i = i + 1
	# Clear things
	plt.clf(); fig_state.clf(); ax_list = []


