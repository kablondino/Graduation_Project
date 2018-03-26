"""
	This file generates the 1D mesh and the 4 cell
	state variables needed for the model (including Diffusivity).
"""

from fipy import Grid1D, CellVariable
from fipy.tools import numerix, dump

from parameters import *


# ----------------- Mesh Generation -----------------------
mesh = Grid1D(nx=config.nx, Lx=config.L)

x = mesh.cellCenters[0] # Cell position
X = mesh.faceCenters[0] # Face position, if needed


# -------------- State Variable Declarations --------------
density = CellVariable(name=r"$n$", mesh=mesh, hasOld=True,\
		unit="1.0e19*m**-3")

temperature = CellVariable(name=r"$T$", mesh=mesh, hasOld=True, unit="eV")

U = CellVariable(name=r"$U$", mesh=mesh, hasOld=True, unit="eV/m**3")

Z = CellVariable(name=r"$Z$", mesh=mesh, hasOld=True, unit="")

Diffusivity = CellVariable(name=r"$D$", mesh=mesh, hasOld=True, unit="m**2/s")

U.setValue(density*temperature / (gamma - 1.0))

# Variables that properly add the units to the respective terms
density_grad_unit = PhysicalField(value=1.0, unit="m**-4")
temp_grad_unit = PhysicalField(value=1/charge.numericValue, unit="eV/m")
Z_transient_unit = PhysicalField(value=1.0, unit="1/s")
Z_diffusion_unit = PhysicalField(value=1.0, unit="m**-2")


# -------------- Other Variable Declarations --------------
# Neutrals density in use for CX friction
n_0 = CellVariable(name=r"$n_0$", mesh=mesh)

# Thermal velocities (most probable)
v_Ti = CellVariable(name=r"$v_{th,i}$", mesh=mesh, unit="m/s")
v_Te = CellVariable(name=r"$v_{th,e}$", mesh=mesh, unit="m/s")

# Poloidal gyro-(Larmor) radii
rho_pi = CellVariable(name=r"$\rho_{\theta i}$", mesh=mesh, unit="m")
rho_pe = CellVariable(name=r"$\rho_{\theta e}$", mesh=mesh, unit="m")

# Banana orbit bounce frequencies
omega_bi = CellVariable(name=r"$\omega_{bi}$", mesh=mesh, unit="1 / s")
omega_be = CellVariable(name=r"$\omega_{be}$", mesh=mesh, unit="1 / s")

# Collision Frequencies
nu_ei = CellVariable(name=r"$\nu_{ei}$", mesh=mesh, unit="1 / s")
nu_ii = CellVariable(name=r"$\nu_{ii}$", mesh=mesh, unit="1 / s")
nu_in0 = CellVariable(name=r"$\nu_{i0}$", mesh=mesh, unit="1 / s")
nu_eff = CellVariable(name=r"$\nu_{eff}$", mesh=mesh, unit="1 / s")
nu_ai = CellVariable(name=r"$\nu_{*i}$", mesh=mesh, unit="")
nu_ae = CellVariable(name=r"$\nu_{*e}$", mesh=mesh, unit="")

## Electron Anomalous Diffusion
D_an = CellVariable(name=r"$D_{an}$", mesh=mesh, unit="m**2 / s")
g_n_an = CellVariable(name=r"$g_n^{an}$", mesh=mesh, unit="A/m")
g_T_an = CellVariable(name=r"$g_T^{an}$", mesh=mesh, unit="A/m")
g_Z_an = CellVariable(name=r"$g_Z^{an}$", mesh=mesh, unit="A/m**2")
Gamma_an = CellVariable(name=r"$\Gamma_e^{an}$", mesh=mesh,\
		unit="1 / (m**2 * s)")

## Charge Exchange Friction
g_n_cx = CellVariable(name=r"$g_n^{cx}$", mesh=mesh, unit="A/m")
g_T_cx = CellVariable(name=r"$g_T^{cx}$", mesh=mesh, unit="A/m")
g_Z_cx = CellVariable(name=r"$g_Z^{cx}$", mesh=mesh, unit="A/m**2")
Gamma_cx = CellVariable(name=r"$\Gamma_i^{cx}$", mesh=mesh,\
		unit="1 / (m**2 * s)")

## Ion Bulk (Parallel) Viscosity
# Consolidated constant to reduce clutter, listed as N on reference
#N = CellVariable(name="Consolidated Bulk Viscosity value", mesh=mesh)
## xi_p integral
#xi_p = CellVariable(name=r"$\xi_\theta$", mesh=mesh)
## xi_t integral
#xi_t = CellVariable(name=r"$\xi_\phi$", mesh=mesh)

g_n_bulk = CellVariable(name=r"$g_n^{\pi\parallel}$", mesh=mesh, unit="A/m")
g_T_bulk = CellVariable(name=r"$g_T^{\pi\parallel}$", mesh=mesh, unit="A/m")
g_Z_bulk = CellVariable(name=r"$g_Z^{\pi\parallel}$", mesh=mesh, unit="A/m**2")
Gamma_bulk = CellVariable(name=r"$\Gamma_i^{\pi\parallel}$", mesh=mesh,\
		unit="1 / (m**2 * s)")

## Ion Orbit Loss
g_OL = CellVariable(name=r"$g^{OL}$", mesh=mesh)
Gamma_OL = CellVariable(name=r"$\Gamma_i^{OL}$", mesh=mesh,\
		unit="1 / (m**2 * s)")

variable_dictionary = {\
		'x': x, 'density': density, 'temperature': temperature, 'Z': Z,\
		'Diffusivity': Diffusivity, 'n_0': n_0, 'v_Ti': v_Ti, 'v_Te': v_Te,\
		'rho_pi': rho_pi, 'rho_pe': rho_pe, 'omega_bi': omega_bi,\
		'omega_be': omega_be, 'nu_ei': nu_ei, 'nu_ii': nu_ii,\
		'nu_in0': nu_in0, 'nu_eff': nu_eff, 'nu_ai': nu_ai, 'nu_ae': nu_ae,\
		'D_an': D_an, 'g_n_an': g_n_an, 'g_T_an': g_T_an, 'g_Z_an': g_Z_an,\
		'Gamma_an': Gamma_an,\
		'g_n_cx': g_n_cx, 'g_T_cx': g_T_cx, 'g_Z_cx': g_Z_cx,\
		'Gamma_cx': Gamma_cx, 'g_OL': g_OL, 'Gamma_OL': Gamma_OL,\
		'g_n_bulk': g_n_bulk, 'g_T_bulk': g_T_bulk, 'g_Z_bulk': g_Z_bulk,\
		'Gamma_bulk': Gamma_bulk
		}

