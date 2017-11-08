from fenics import *
import sys
sys.path.append("/home/kabv/Documents/Masters/Graduation_Project/Model_Source/")
from parameters import *
import coefficients as coef

# Create a mesh
nx = 100
L = 1.0
mesh = IntervalMesh(nx, 0, L)
V = VectorFunctionSpace(mesh, 'P', 1)

# Set up independent and dependent variables
x = SpatialCoordinate(mesh)
#z = TrialFunction(V)
v = TestFunction(V)
Z = Function(V)


## ---------------- Staps Values --------------------------
epsilon = 1.0/25.0
mu = 1.0/20.0
zeta = 0.5
c_n = -1.1
c_T = -0.9

D_min = 2.0/5.0
D_max = 2.0
alpha_sup = 0.5

# Values for the polynomial G; g1 = a, g2 = b, g3 = c
g1 = 1.5
g2 = 2.0
g3 = -1.0
Z_S = -1.5

Gamma_c = -4.0/5.0
gamma = 5.0/3.0
q_c = -4.0
lambda_n = 1.25
lambda_T = 1.5
lambda_Z = 1.25

# Diffusivity
D_shear = D_min + (D_max - D_min) / (1 + alpha_sup*(Z)**2)
#D_shear = Expression('D_min + (D_max - D_min) / (1 + alpha_sup*(Z)**2)', degree=1, D_min=D_min, D_max=D_max, alpha_sup=alpha_sup,Z=Z)
print D_shear
def boundary(x, on_boundary):
	return on_boundary

# WRTIE A PROPER BOUNDARY CONDITION!
#bc = 



