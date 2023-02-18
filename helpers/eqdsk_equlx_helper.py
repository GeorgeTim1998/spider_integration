import numpy as np
import math as mt
from math import pi

def write_np_array_to_file(file, array, number_format):
  columns = 5
  rows = mt.floor(len(array)/5)
  sep_index = columns*rows

  arr5 = array[0: sep_index]
  arrleftover = array[sep_index: :]

  arr5 = np.resize(arr5, (rows, columns))
  arrleftover = np.resize(arrleftover, (1, len(arrleftover)))
  np.savetxt(file, arr5, fmt=number_format, delimiter='')
  np.savetxt(file, arrleftover, fmt=number_format, delimiter='')
  
def restore_pres_n_fpol(um, up, meshr, pprime, ffprim, bcentr, rcentr):
  um = -um / (2*pi)
  up = -up / (2*pi)
  dpsi = (um - up)/(meshr-1)

  pres = np.zeros(len(pprime))
  fpol = np.zeros(len(pprime))
  fpol[-1] = bcentr * rcentr

  my_steps = np.flip(np.array(range(0, meshr-1)))
  for i in my_steps:
    pprim_c = -0.5 * (pprime[i+1] + pprime[i])
    fprim_c = -0.5 * (ffprim[i+1] + ffprim[i])
    
    pres[i]=(pres[i+1] + pprim_c * 2*pi*dpsi) # for some reason pres needs 2pi multiplicator to be actual pressure
    fpol[i]=mt.sqrt(fpol[i+1]**2 + 2*fprim_c * 2*pi*dpsi)

  return pres, fpol  
