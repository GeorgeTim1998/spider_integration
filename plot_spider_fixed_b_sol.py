import numpy as np
from math import pi
import support as sup
import matplotlib.pyplot as pyplot
from fenics import *
from scipy.interpolate import interp2d, LinearNDInterpolator
from fenics_support import countour_plot_via_mesh, assign_const_to_nan_in_expression, plot_1D, interpolate_spider_data_on_function_space
from helpers import eqdsk_equlx_helper as eq

path = '/media/george/part/Spider'
working_folder = 'WK_PRIME_restore_q'
pic_path = "Pics/%s.png" % working_folder
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

def restore_funcpsi(psi, dfuncdpsi, funcb=0):
  psi = np.flip(psi)
  dfuncdpsi = np.flip(dfuncdpsi)

  funcpsi = np.zeros(len(psi))
  
  funcpsi[0] = funcb
  for i in range(1, len(psi)):
    funcpsi[i] = dfuncdpsi[i] * (psi[i]-psi[i-1]) + funcpsi[i-1]
    
  return np.flip(funcpsi)
  
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

#%% Plot data needed from Spider
psi = psi_max * (1 - sqrt_psi_norm**2) # This is magnetic flux/2pi. In spider flux is used/ Multiply by 2pi

ppsi = restore_funcpsi(2*pi*psi, dpdpsi)
fpsi = restore_funcpsi(2*pi*psi, dfdpsi, fvac)

I = np.ones((psi_size, spacial_size))
r_mesh = ro*(rb - rc) + I*rc
z_mesh = ro*(zb - zc) + I*zc

psi_mesh = (I.transpose() * psi).transpose()
p_mesh = (I.transpose() * ppsi).transpose()
f_mesh = (I.transpose() * fpsi).transpose()
dfdpsi_mesh = (I.transpose() * dfdpsi).transpose()
dpdpsi_mesh = (I.transpose() * dpdpsi).transpose()
q_mesh = (I.transpose() * q).transpose()

figure = pyplot.contour(r_mesh, z_mesh, psi_mesh, 10)
pyplot.colorbar(figure).set_label("\u03C8(r, z), Вб")
pyplot.gca().set_aspect("equal")
pyplot.grid(True)
pyplot.xlim(0.8, 2.4)
pyplot.ylim(-1, 1)
pyplot.xlabel("r, м")
pyplot.ylabel("z, м")
pyplot.savefig(pic_path, dpi=240, bbox_inches="tight")
pyplot.close()
print("\n", 1)

#%% try to find what to do...
pprim = eq.default_pprime()
ffprim = eq.default_ffprim()

multiplicator_p = pprim[0]/dpdpsi_mesh[0,0]
multiplicator_f = ffprim[0]/f_mesh[0,0]*dpdpsi_mesh[0,0]


# plot_1D(r_mesh[:, 0], p_mesh[:, 0], xlabel='coord', ylabel='p_mesh')
# plot_1D(r_mesh[:, 0], multiplicator_p*dpdpsi_mesh[:, 0], xlabel='coord', ylabel='dpdpsi')
# plot_1D(r_mesh[:, 0], f_mesh[:, 0], xlabel='coord', ylabel='f_mesh')
# plot_1D(r_mesh[:, 0], 0.25*multiplicator_f*f_mesh[:, 0]*dfdpsi_mesh[:, 0], xlabel='coord', ylabel='dfdpsi')


#%% Import Spider solution to fenics
folder = sup.xml_files_folder()
filename = 'test'

xml_file = "%s/%s.xml" % (folder, filename)

gmsh = Mesh(xml_file)
V = FunctionSpace(gmsh, 'Lagrange', 1)

psi = interpolate_spider_data_on_function_space(r_mesh, z_mesh, psi_mesh, V)

countour_plot_via_mesh(gmsh, psi, levels=50, colorbar=True, grid=True, xlim=[0.8, 2.4], ylim=[-1, 1])