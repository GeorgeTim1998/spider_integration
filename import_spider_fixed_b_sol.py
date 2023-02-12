import numpy as np
from math import pi
import support as sup
import matplotlib.pyplot as pyplot

path = '/media/george/part/Spider'
working_folder = 'WK_T15MD_fixed_boundary'
wr_file = 'spik.wr'

path_to_file = "%s/%s/%s" % (path, working_folder, wr_file)

def first_line(file):
  psi_size, spacial_size, size, psi_max = next(file).split()
  
  return int(psi_size), int(spacial_size), int(size), float(psi_max)

def second_line(file):
  psi_min, end_entry = next(file).split()
  
  return float(psi_min), int(end_entry)

def return_and_delete_range(data, data_range):
  return_array = data[0:data_range]
  data = np.delete(data, np.s_[0:data_range])
  
  return return_array, data

def restore_ppsi(psi, dpdpsi, pmin=0):
  psi = np.flip(psi)
  dpdpsi = np.flip(dpdpsi)

  ppsi = np.zeros(len(psi))
  
  ppsi[0] = pmin
  for i in range(1, len(psi)):
    ppsi[i] = dpdpsi[i] * (psi[i]-psi[i-1]) + ppsi[i-1]
    
  return np.flip(ppsi)
  
#%% Get data from file
with open(path_to_file, 'r') as file:
  psi_size, spacial_size, size, psi_max = first_line(file)
  psi_min, end_entry = second_line(file)

  sqrt_psi_norm = np.zeros(psi_size)
  dpdpsi = np.zeros(psi_size)
  dfdpsi = np.zeros(psi_size)
  
  rc = np.zeros(spacial_size)
  zc = np.zeros(spacial_size)
  
  rb = np.zeros(spacial_size)
  zb = np.zeros(spacial_size)
  
  ro = np.zeros(psi_size*spacial_size)
  
  q = np.zeros(psi_size)
  fvac = 0
  
  data = []
  for line in file: # read rest of lines
    data.append([float(num) for num in line.split()])

data = np.array(data).reshape(np.size(data))

sqrt_psi_norm, data = return_and_delete_range(data, psi_size)

dpdpsi, data = return_and_delete_range(data, psi_size) 
dfdpsi, data = return_and_delete_range(data, psi_size) 

rc, data = return_and_delete_range(data, spacial_size) 
zc, data = return_and_delete_range(data, spacial_size) 
rb, data = return_and_delete_range(data, spacial_size) 
zb, data = return_and_delete_range(data, spacial_size) 

ro, data = return_and_delete_range(data, psi_size*spacial_size) 

q, data = return_and_delete_range(data, psi_size)
fvac, data = return_and_delete_range(data, 1)

ro = ro.reshape(psi_size, spacial_size)

#%% Restore data needed fo fenics
psi = psi_max * (1 - sqrt_psi_norm**2) # This is magnetic flux/2pi. In spider flux is used/ Multiply by 2pi

ppsi = restore_ppsi(2*pi*psi, dpdpsi)
fpsi = restore_ppsi(2*pi*psi, dfdpsi, fvac)

I = np.ones((psi_size, spacial_size))
r_mesh = ro*(rb - rc) + I*rc
z_mesh = ro*(zb - zc) + I*zc

sqrt_psi_norm = (I.transpose() * psi).transpose()

figure = pyplot.contour(r_mesh, z_mesh, sqrt_psi_norm)
pyplot.colorbar(figure).set_label("\u03C8(r, z), Вб")
pyplot.gca().set_aspect("equal")
pyplot.xlim(0.8, 2.4)
pyplot.ylim(-1, 1)
pyplot.xlabel("r, м")
pyplot.ylabel("z, м")
pyplot.savefig('Pics/1.png', dpi=240, bbox_inches="tight")
print(1)
# pyplot.scatter(r_start_array, z_start_array)
# pyplot.scatter(r_end_array, z_end_array)
# pyplot.show()