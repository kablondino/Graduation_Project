import matplotlib.pyplot as plt
import numpy
from matplotlib.ticker import FormatStrFormatter
import csv

import os
import sys

# Import DIRECTORY NAME in which to do the bidding
data_directory = sys.argv[1]


# Create list of files
file_list = []
big_data = []

# Function to remove first and last y-axis tick labels
def remove_yticks(*the_axes):
	for axis in the_axes:
		plt.setp(axis.get_yticklabels()[0], visible=False)
		plt.setp(axis.get_yticklabels()[-1], visible=False)

def min_max(Afile_list, Abig_data, Acolumn):
	current_column.append(Abig_data[0][:, Acolumn])
	return current_column
#	the_axis.ylim(min(), max())

# Make the file_list and ALL the data
for filename in os.listdir('./'+str(data_directory)):
	file_list.append(filename)
#	big_data.append(numpy.loadtxt(data_directory+'/'+filename, skiprows=1))


# Sort the file_list
file_list.sort()
x, density, temperature, Z, Diffusivity = [], [], [], [], []


i = 0 # Looping counter
for filename in file_list:
	filename_sans_ext = int(os.path.splitext(filename)[0])
#	data = numpy.loadtxt(data_directory+'/'+filename, skiprows=1)

	# The data
#	x = big_data[0][:,0]
	x.append(numpy.genfromtxt(data_directory+'/'+filename, skip_header=1, unpack=True))
#	density = big_data[0][:,1]
#	temperature = big_data[0][:,2]
#	Z = big_data[0][:,3]
#	Diffusivity = big_data[0][:,4]
	print x

	# Generate the figure
	fig = plt.subplots(2, 1, sharex=True, squeeze=True, figsize=(12,12))[0]
	fig.subplots_adjust(hspace=0)
	
	i = i+1
	if i == 2:
		break

	# TOP PLOT
#	ax_list.append(plt.subplot(2,1,1))
#	density_plot = ax_list[0].plot(x, density, label=r"$n$", color='blue')
#	ax_list[0].set_ylabel(r"$n$", fontsize='large', rotation=0, labelpad=20)
#	ax_list[0].tick_params(axis='x', labelbottom='off') # Remove top plot's x-axis labels
#
#	ax_list.append(ax_list[0].twinx())
#	temp_plot = ax_list[1].plot(x, temperature, label=r"$T$", color='red')
#	ax_list[1].set_ylabel(r"$T$", fontsize='large', rotation=0, labelpad=20)
#
#	ax_list[1].set_yticks(numpy.linspace(ax_list[1].get_yticks()[0], ax_list[1].get_yticks()[-1],\
#			len(ax_list[0].get_yticks())))
#
#	ax_list[0].grid(True)
#	ax_list[1].grid(True)
#	ax_list[1].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
#
##	# BOTTOM PLOT
##	ax_list[2] = plt.subplot(2,1,2)
##	Z_plot = ax_list[2].plot(x, Z, label=r"$Z$", color='green')
##	ax_list[2].set_ylabel(r"$Z$", fontsize='large', rotation=0, labelpad=20)
##
##	ax_list[4] = ax_list[2].twinx()
##	D_plot = ax_list[4].plot(x, Diffusivity, label=r"$D$", color='orange')
##	ax_list[4].set_ylabel(r"$D$", fontsize='large', rotation=0, labelpad=20)
##
##	ax_list[2].set_yticks(numpy.linspace(ax_Z.get_yticks()[0], ax_Z.get_yticks()[-1],\
##			len(ax_list[4].get_yticks())))
##
##	ax_list[2].grid(True)
##	ax_list[2].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
##	ax_list[4].grid(True)
##	ax_list[4].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
##
##	ax_list[2].set_xlabel(r"$x$", fontsize='large')
##
##	remove_yticks(ax_list[0], ax_list[1], ax_list[2], ax_list[4])
##
##	top_plot = density_plot + temp_plot
##	top_labels = [l.get_label() for l in top_plot]
##	ax_list[1].legend(top_plot, top_labels, loc='best', fontsize='medium')
##	bottom_plot = Z_plot + D_plot
##	bottom_labels = [l.get_label() for l in bottom_plot]
##	ax_list[4].legend(bottom_plot, bottom_labels, loc='best')
##
##
##	fig.suptitle(r"$\Gamma_c = -1.0\times 10^{18}~, ~~ D = 1 / (1 + 0.01\,Z^2 + 0.001\,(Z^\prime)^{-2})$" +"\n"+ "$t = " +str(filename_sans_ext)+ "$")
##	fig.tight_layout(pad=0.2, w_pad=0.0)
##	plt.subplots_adjust(top=0.9)
##
###	plt.savefig(data_directory +'/'+ filename_sans_ext +'.png')
##
##	print str(i) + "\t| Saved " + str(filename_sans_ext) + ".png"
##	i = i + 1
##	plt.clf()
##
###plt.show()
##
