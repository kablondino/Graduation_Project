from fipy import Grid1D, CellVariable, Viewer
from fipy.tools import numerix, dump

import numpy
from scipy.special import wofz


# ----------------- Mesh Generation -----------------------
L = 4.0
nx = 100
mesh = Grid1D(nx=nx, Lx=L)
x = mesh.cellCenters[0] # Cell position

plasma_disp_real = CellVariable(name=r"Re$(X(z))$", mesh=mesh)
plasma_disp_imag = CellVariable(name=r"Im$(X(z))$", mesh=mesh)

full_plasma_disp = numpy.array(1j*numerix.sqrt(numerix.pi)\
		* numerix.exp(-x**2)*wofz(x))
print full_plasma_disp

plasma_disp_real.setValue(full_plasma_disp.real)
plasma_disp_imag.setValue(full_plasma_disp.imag)

initial_viewer = Viewer((plasma_disp_real, plasma_disp_imag),\
		xmin=0.0, legend='best')
raw_input("Pause for Initial Conditions")

