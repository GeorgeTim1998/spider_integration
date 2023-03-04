# All that is needed here is probably dp/dpsi and df/dpsi
import numpy as np
from helpers import eqdsk_equlx_helper as eq
from math import pi, sqrt, floor
import matplotlib.pyplot as pyplot
import create_gmsh_mesh_from_points as gmsh
from support import lao_hash, print_colored
from pathlib import Path

number_format = "%16.9E" 
formating = "%s\n" % (number_format * 5)

#%% Saving to folders

spider_folder = "WK_linear_profs_no2pi_vs_fenics_slimmest"
files_dir = "Files"
my_dir = "Files/%s" % spider_folder

#%% Used for comparing with Spider
Re = lao_hash()['Re'] # ellipse center which sometimes can be R0 in your writings
ell_a = lao_hash()['a']
R_plasma_boundary = Re + ell_a
elongations = lao_hash()['E']

F2_0 = lao_hash()['F2_0']
I = lao_hash()['I']

betta_pol = lao_hash()['betta_pol']
betta_tor = lao_hash()['betta_tor']

alpha = 0.1
kappa = 1e-2

for elongation in elongations:
#%% Define all data
  ACASE48, SPIDER, IDUM, MESHR, MESHZ = 'KIAM', 'SPIDER', 3, 128, 257
  RDIM, ZDIM, RCENTR, RLEFT, ZMID = 0, 0, Re, 0, 0
  RMAXIS, ZMAXIS, UM, UP, BCENTR =  Re, 0, 1, 0, F2_0**0.5 / R_plasma_boundary
  CURRENT, RX1, ZX1, RX2, ZX2 = I/2/pi, 0, 0, 0, 0 # spider uses current for calculating actual psi...
  ZMAXIS, UXN, UX1, UX2, XDUM = 0, 0, 0, 0, 0

  pprime, ffprim = eq.exponential_derivarives_profiles(Re, betta_pol, betta_tor, alpha, kappa, MESHR)

  pres, fpol = eq.restore_pres_n_fpol(UM, UP, MESHR, pprime, ffprim, BCENTR, RCENTR)

  u = np.zeros(MESHR * MESHZ)
  q = np.zeros(len(pprime))

  NXB,NBLM = 89, 89 

  boundary = eq.ellipse_boundary(Re, ell_a, elongation, NXB)
  
#%% Write to eqdsk_equilx
  folder = "%s/%.1f" % (files_dir, elongation)
  Path(folder).mkdir(parents=True, exist_ok=True)
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

  print_colored("eqdsk_equilx saved to:", color='green', white_str="%s/%s" % (folder, filename))

  # eq.plot_pres_fpol_n_ders(pres, fpol, pprime, ffprim)

  gmsh.create_gmsh_mesh_from_points(boundary, ell_a, elongation)