"""
	This file generates the 1D mesh and the 4 cell
	variables needed for the model. It can be argued that
	Diffusivity should not be one of them.
"""
from fipy import Grid1D, CellVariable

# ----------------- Mesh Generation -----------------------
nx = 100
L = 5.0
mesh = Grid1D(nx=nx, Lx=L)

x = mesh.cellCenters[0]

# ----------------- Variable Declarations -----------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True)

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True)

Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True)

Diffusivity = CellVariable(name=r"$D$", mesh=mesh, hasOld=True)

