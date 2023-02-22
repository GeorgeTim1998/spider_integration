# All that is needed here is probably dp/dpsi and df/dpsi
import numpy as np
from helpers import eqdsk_equlx_helper as eq
from math import pi, sqrt, floor
import matplotlib.pyplot as pyplot
import create_gmsh_mesh_from_points as gmsh
from support import lao_hash

number_format = "%16.9E" 
formating = "%s\n" % (number_format * 5)

#%% Used for comparing with Spider
Re = lao_hash()['Re'] # ellipse center which sometimes can be R0 in your writings
ell_a = lao_hash()['a']
E = 1
p0 = 10000
psi0 = 0.5 # на это я обезразмеривал
psi_max = 0.002613455467247852 # это максимальная psi/ Поток на 1 рад  
F2_0 = 0.25
Fpl_vs_Fvac_ratio = 0.1

#%%
ACASE48, SPIDER, IDUM, MESHR, MESHZ = 'KIAM', 'SPIDER', 3, 128, 257
RDIM, ZDIM, RCENTR, RLEFT, ZMID = 0, 0, Re, 0, 0
RMAXIS, ZMAXIS, UM, UP, BCENTR =  Re, 0, (psi_max)*2*pi, 0, F2_0**0.5 / Re
CURRENT, RX1, ZX1, RX2, ZX2 = 0.15E+07, 0, 0, 0, 0
ZMAXIS, UXN, UX1, UX2, XDUM = 0, 0, 0, 0, 0

pprime = eq.pprime_linear_profile(p0, (psi0)*2*pi, MESHR)
ffprim = eq.ffprim_linear_profile(F2_0, (psi0)*2*pi, MESHR)

pres, fpol = eq.restore_pres_n_fpol(UM, UP, MESHR, pprime, ffprim, BCENTR, RCENTR)

u = np.zeros(MESHR * MESHZ)
q = np.zeros(len(pprime))

NXB,NBLM = 89, 89

boundary = eq.ellipse_boundary(Re, ell_a, E, NXB)
#%% plot boundary
# pyplot.scatter(boundary.reshape((-1, 2))[:, 0], boundary.reshape((-1, 2))[:, 1])
# pyplot.gca().set_aspect("equal")
# pyplot.show()

#%% 
folder = 'Files/fenics_vs_spider'
filename = 'eqdsk_equilx'
with open("%s/%s" % (folder, filename), 'w') as file:
  file.write("%8s%8s%36d%4d%4d\n" % (ACASE48, SPIDER, IDUM, MESHR, MESHZ))
  file.write(formating % (RDIM, ZDIM, RCENTR, RLEFT, ZMID))
  file.write(formating % (RMAXIS, ZMAXIS, UM, UP, BCENTR))
  file.write(formating % (CURRENT, RX1, ZX1, RX2, ZX2))
  file.write(formating % (ZMAXIS, UXN, UX1, UX2, XDUM))
  
  for array in [fpol, pres, ffprim, pprime, u, q]:
    eq.write_np_array_to_file(file, array, number_format)
  
  file.write("%5d%5d\n" % (NXB, NBLM))
  
  eq.write_np_array_to_file(file, boundary, number_format)
  eq.write_np_array_to_file(file, boundary, number_format)

eq.plot_pres_fpol_n_ders(pres, fpol, pprime, ffprim)

gmsh.create_gmsh_mesh_from_points(boundary)