import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt
import numpy

import os
import sys

# Import DIRECTORY NAME in which to do the bidding
data_directory = sys.argv[1]


# Create list of files
file_list = []

for filename in os.listdir('./'+str(data_directory)):
	if filename.endswith(".tsv"):
		file_list.append(filename)

i = 0 # Looping counter
for filename in file_list:
	data = numpy.loadtxt(data_directory+'/'+filename, skiprows=1)
	x = data[:,0]
	density = data[:,1]
	temperature = data[:,2]
	Z = data[:,3]
	Diffusivity = data[:,4]

	filename_sans_ext = os.path.splitext(filename)[0]

	plt.plot(x, density)
	plt.plot(x, temperature)
	plt.plot(x, -Z)
	plt.plot(x, Diffusivity)
	plt.legend([r"$n$", r"$T$", r"$Z$", r"$D$"], loc='upper right')
	plt.title(r"$\Gamma_c = -0.8~, ~~ D \sim 1 / (1 + (Z^\prime)^{-2})$" +"\n"+ "$t = "+ str(filename_sans_ext)+ "$")
	plt.savefig(data_directory +'/'+ filename_sans_ext +'.png')

#	print str(i) + "\t| Saved " + str(filename_sans_ext) + ".png"
	i = i + 1
	plt.clf()

