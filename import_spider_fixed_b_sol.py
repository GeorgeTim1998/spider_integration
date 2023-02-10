import numpy as np
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

# for i, row in enumerate(ro)
I = np.ones((psi_size, spacial_size))
r_mesh = ro*(rb - rc) + I*rc
z_mesh = ro*(zb - zc) + I*zc

print(1)
# pyplot.scatter(r_start_array, z_start_array)
# pyplot.scatter(r_end_array, z_end_array)
# pyplot.show()