import numpy as np
from helpers import eqdsk_equlx_helper as eq
from math import pi, sqrt, floor
import matplotlib.pyplot as pyplot

number_format = "%16.9E" 
formating = "%s\n" % (number_format * 5)

ACASE48, SPIDER, IDUM, MESHR, MESHZ = 'KIAM', 'SPIDER', 3, 128, 257
RDIM, ZDIM, RCENTR, RLEFT, ZMID = 2.5, 4, 1.48, 0.5, 0
RMAXIS, ZMAXIS, UM, UP, BCENTR =  0.162692725E+01, 0, 0.602393949E+00, 0.157520013E+00, 0.200000000E+01
CURRENT, RX1, ZX1, RX2, ZX2 = 0.15E+07, 0, 0, 0, 0
ZMAXIS, UXN, UX1, UX2, XDUM = 0, 0, 0, 0, 0

pprime = eq.default_pprime()
ffprim = eq.default_ffprim()

pres, fpol = eq.restore_pres_n_fpol(UM, UP, MESHR, pprime, ffprim, BCENTR, RCENTR)

u = np.zeros(MESHR * MESHZ)
q = np.zeros(len(pprime))

NXB,NBLM = 89, 89

boundary = eq.default_plasma_boundary()

folder = 'Files'
filename = 'eqdsk_equilx'
with open("%s/%s" % (folder, filename), 'w') as file:
  file.write("%8s%8s%36d%4d%4d\n" % (ACASE48, SPIDER, IDUM, MESHR, MESHZ))
  file.write(formating % (RDIM, ZDIM, RCENTR, RLEFT, ZMID))
  file.write(formating % (RMAXIS, ZMAXIS, UM, UP, BCENTR))
  file.write(formating % (CURRENT, RX1, ZX1, RX2, ZX2))
  file.write(formating % (ZMAXIS, UXN, UX1, UX2, XDUM))
  
  for array in [fpol, pres, pprime, ffprim, u, q]:
    eq.write_np_array_to_file(file, array, number_format)
  
  file.write("%5d%5d\n" % (NXB, NBLM))
  
  eq.write_np_array_to_file(file, boundary, number_format)
  eq.write_np_array_to_file(file, boundary, number_format)

exit()
psi = np.linspace(0, 1, len(fpol))

ax1= pyplot.subplot(231)
pyplot.scatter(psi, pres)
pyplot.ylabel('pres')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(232)
pyplot.scatter(psi, pprime)
pyplot.ylabel('pprim')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(233)
pyplot.scatter(psi, pres)
pyplot.ylabel('my_pres')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(234)
pyplot.scatter(psi, fpol)
pyplot.ylabel('fpol')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(235)
pyplot.scatter(psi, ffprim)
pyplot.ylabel('ffprim')
pyplot.xlabel('psi')
pyplot.grid(True)

ax = pyplot.subplot(236)
pyplot.scatter(psi, fpol)
pyplot.ylabel('my_fpol')
pyplot.xlabel('psi')
pyplot.grid(True)

pyplot.show()

pyplot.figure()
pyplot.scatter(psi, pres/pres)
pyplot.ylabel('pres multuplicator')
pyplot.xlabel('psi')
pyplot.grid(True)
pyplot.savefig("Pics/pres_vs_my_pres.png", dpi=240)

pyplot.figure()
pyplot.scatter(psi, fpol/fpol)
pyplot.ylabel('fpol multuplicator')
pyplot.xlabel('psi')
pyplot.grid(True)
pyplot.savefig("Pics/fpol_vs_my_fpol.png", dpi=240)

print(1)
