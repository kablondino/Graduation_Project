from fipy import Grid1D, CellVariable

nx = 100
L = 5.0
mesh = Grid1D(nx=nx, Lx=L)

x = mesh.cellCenters[0]

# ----------------- Variable Declarations -----------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True)

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True)

Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True)

Diffusivity = CellVariable(name=r"$D$", mesh=mesh, hasOld=True)

