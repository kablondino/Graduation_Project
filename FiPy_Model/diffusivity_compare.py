from fipy import Grid1D, CellVariable, Viewer, MatplotlibViewer, TSVViewer
from matplotlib import pyplot

from fipy.tools import numerix

nx = 1000
L = 5.0
mesh = Grid1D(nx=nx, Lx=L)

x = mesh.cellCenters[0]

# Parameters
zeta = 0.5
Gamma_c = -4.0/5.0
q_c = -4.0
gamma = 5.0/3.0

lambda_n = 5.0/4.0
lambda_T = 3.0/2.0
lambda_Z = 5.0/4.0

D_max = 2.0
D_min = 2.0/5.0

epsilon = 1.0/25.0
mu = 1.0/20.0

c_n = -1.1
c_T = -0.9

a = 3.0/2.0
b = 2.0
c = -1.0
Z_S = -3.0/2.0


# ----------------- Variable Declarations -----------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True)

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True)

Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True)

## Diffusivity
D_Zohm = CellVariable(name="Zohm", mesh=mesh, hasOld=True)
D_Staps = CellVariable(name="Staps", mesh=mesh, hasOld=True)
D_Shear = CellVariable(name="Flow-Shear", mesh=mesh, hasOld=True)


# ------------- Initial Conditions and Diffusivity---------
Z0 = Z_S*(1 - numerix.tanh((L*x - L) / 2))
Z.setValue(Z0)

# Zohm's model
D_Zohm.setValue((D_max + D_min) / 2.0 + ((D_max - D_min)*numerix.tanh(Z)) / 2.0)
# Stap's Model
alpha_sup = 0.5
D_Staps.setValue(D_min + (D_max - D_min) / (1.0 + alpha_sup*(Z.grad.mag)**2))
# Flow-Shear Model
a1, a2, a3 = 0.7, 1.25, 0.5
D_Shear.setValue(D_min + (D_max - D_min) / (1.0 + a1*Z**2 + a2*Z*(Z.grad) + a3*(Z.grad)**2))


# ----------------- Boundary Conditions -------------------
# Z Boundary Conditions:
#	d/dx(Z(0,t)) == Z / lambda_Z
#	mu*D/epsilon * d/dx(Z(L,t)) == 0
#	d^2/dx^2(Z(0,t)) == 0
Z.faceGrad.constrain(Z.faceValue / lambda_Z, mesh.facesLeft)
Z.constrain(0.0, mesh.facesRight)

Z.grad.faceGrad.constrain(0.0, mesh.facesLeft)


# ----------------- PDE Declarations ----------------------

#initial_viewer = MatplotlibViewer((D_Zohm, D_Staps, D_Shear), title=r"Diffusivity Models: FiPy", xmin=0.0, xmax=3.0)
#pyplot.grid(True)
#pyplot.axes().set_aspect('equal')

viewer = None
if __name__ == '__main__':
	try:
		viewer = Viewer(vars=D_Zohm, datamin=0.0, datamax=3.0)
		viewer.plotMesh()
	except: print "Unable to create viewer"

#TSVViewer(vars=(D_Zohm, D_Staps, D_Shear)).plot(filename="Diffusivity_data.tsv")

raw_input("Pause for Initial")


