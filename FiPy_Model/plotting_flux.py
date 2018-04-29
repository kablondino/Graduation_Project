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
			delimiter='\t', unpack=True, skiprows=1))

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
	x, density, temperature, Z, Diffusivity, D_an, Gamma_an, Gamma_cx,\
			D_bulk, Gamma_bulk, Gamma_ol = data[:]

	# Generate the figure
	fig_flux = plt.figure()
	fig_flux = plt.subplots(4, 1, sharex=True, squeeze=True, figsize=(12,20))[0]
	fig_flux.subplots_adjust(hspace=0)


	## TOP FLUX PLOT (Gamma_an)
	ax_list.append(plt.subplot(4,1,1))
	Gamma_an_plot = ax_list[0].plot(x, Gamma_an, label=r"$\Gamma_e^{an}$",\
			color='darkcyan', linewidth=2)
	ax_list[0].set_ylabel(r"$\Gamma_e^{an}$", fontsize='large',\
			rotation=0, labelpad=20)
	ax_list[0].tick_params(axis='x', labelbottom='off')
	plt.ylim((the_mins[6], the_maxs[6]))
	remove_yticks(ax_list[0])

	ax_list[0].grid(True)
	ax_list[0].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))


	## SECOND FLUX PLOT (Gamma_cx)
	ax_list.append(plt.subplot(4,1,2))
	Gamma_cx_plot = ax_list[1].plot(x, Gamma_cx, label=r"$\Gamma_i^{cx}$",\
			color='darkcyan', linewidth=2)
	ax_list[1].set_ylabel(r"$\Gamma_i^{cx}$", fontsize='large',\
			rotation=0, labelpad=20)
	ax_list[1].tick_params(axis='x', labelbottom='off')
	plt.ylim((the_mins[7], the_maxs[7]))
	remove_yticks(ax_list[1])

	ax_list[1].grid(True)
	ax_list[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))


	## THIRD FLUX PLOT (Gamma_bulk and D_bulk)
	ax_list.append(plt.subplot(4,1,3))
	Gamma_bulk_plot = ax_list[2].plot(x, Gamma_bulk,\
			label=r"$\Gamma_i^{\pi\parallel}$", color='darkcyan', linewidth=2)
	ax_list[2].set_ylabel(r"$\Gamma_i^{\pi\parallel}$", fontsize='large',\
			rotation=0, labelpad=20)
	ax_list[2].tick_params(axis='x', labelbottom='off')
	plt.ylim((the_mins[9], the_maxs[9]))
	remove_yticks(ax_list[2])

	ax_list[2].grid(True)
	ax_list[2].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))


	## FOURTH FLUX PLOT (Gamma_ol)
	ax_list.append(plt.subplot(4,1,4))
	Gamma_ol_plot = ax_list[3].plot(x, Gamma_ol, label=r"$\Gamma_i^{ol}$",\
			color='darkcyan', linewidth=2)
	ax_list[3].set_ylabel(r"$\Gamma_i^{ol}$", fontsize='large',\
			rotation=0, labelpad=20)
	ax_list[3].tick_params(axis='x', labelbottom='off')
	plt.ylim((the_mins[10], the_maxs[10]))
	remove_yticks(ax_list[3])

	ax_list[3].grid(True)
	ax_list[3].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

	# The x-axis label and ticks
	ax_list[3].grid(True)
	ax_list[3].set_xlabel(r"$x$", fontsize='large')



	for k in range(len(ax_list)):
		ax_list[k].yaxis.tick_right()

	fig_flux.suptitle(r"$\Gamma_c = -1.0\times 10^{18}~, ~~ D = 1 / (1 + 0.01\,Z^2 + 0.001\,(Z^\prime)^{2})$" +"\n"+ "$t = " +str(int(filename_sans_ext))+ "$",\
			fontsize=22)
	fig_flux.tight_layout(pad=0.1, w_pad=0.0)
	plt.subplots_adjust(top=0.92)

	fig_flux.savefig(data_directory +'/Gamma'+ filename_sans_ext +'.png')

#	print str(i) + "\t| Saved " + str(filename_sans_ext) + ".png"
	i = i + 1
	fig_flux.show()
	raw_input("ADSF")
	# Clear things
	plt.clf(); fig_flux.clf(); ax_list = []
	break


