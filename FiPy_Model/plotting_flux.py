import matplotlib.pyplot as plt
import numpy

import os
import sys

# Import DIRECTORY NAME in which to do the bidding
#data_directory = sys.argv[1]


# Create list of files
file_list = []

#for filename in os.listdir('./'+str(data_directory)):
#	file_list.append(filename)

i = 0 # Looping counter
#for filename in file_list:
data = numpy.loadtxt('0021.tsv', skiprows=1)
x = data[:,0]
density = data[:,1]
temperature = data[:,2]
Z = data[:,3]
Diffusivity = data[:,4]

#	filename_sans_ext = os.path.splitext(filename)[0]

fig = plt.subplots(2, 1, sharex=True)[0]
fig.subplots_adjust(hspace=0)

ax_n = plt.subplot(2,1,1)
ax_n.plot(x, density, color='blue')
ax_n.set_ylabel(r"$n$", fontsize='x-large')

ax_T = ax_n.twinx()
ax_T.plot(x, temperature, color='red')
ax_T.set_ylabel(r"$T$", fontsize='x-large')
print ax_n.get_yticks()
print ax_T.get_yticks()

ax_T.set_yticks(numpy.linspace(ax_T.get_yticks()[0], ax_T.get_yticks()[-1],\
		len(ax_n.get_yticks())))

print ax_T.get_yticks()

ax_n.grid(True)
ax_T.grid(True)

ax_Z = plt.subplot(2,1,2)
ax_Z.plot(x, Z, color='green')

ax_D = ax_Z.twinx()
ax_D.plot(x, Diffusivity, color='orange')

ax_Z.set_yticks(numpy.linspace(ax_Z.get_yticks()[0], ax_Z.get_yticks()[-1],\
		len(ax_D.get_yticks())))

ax_Z.grid(True)
ax_D.grid(True)

fig.suptitle(r"$\Gamma_c = -1.0~, ~~ D \sim 1 / (1 + (Z^\prime)^{-2})$" +"\n"+ "$t = " +"$")

#	plt.savefig(data_directory +'/'+ filename_sans_ext +'.png')

#	print str(i) + "\t| Saved " + str(filename_sans_ext) + ".png"
i = i + 1
#plt.clf()

plt.show()

