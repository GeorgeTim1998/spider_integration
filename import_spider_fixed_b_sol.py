import numpy as np
import support as sup

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
  
  r_start_array = np.zeros(spacial_size)
  z_start_array = np.zeros(spacial_size)
  
  r_end_array = np.zeros(spacial_size)
  z_end_array = np.zeros(spacial_size)
  
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
  
  r_start_array, data = return_and_delete_range(data, spacial_size) 
  z_start_array, data = return_and_delete_range(data, spacial_size) 
  r_end_array, data = return_and_delete_range(data, spacial_size) 
  z_end_array, data = return_and_delete_range(data, spacial_size) 
  
  ro, data = return_and_delete_range(data, psi_size*spacial_size) 
  
  q, data = return_and_delete_range(data, psi_size)
  fvac, data = return_and_delete_range(data, 1)
