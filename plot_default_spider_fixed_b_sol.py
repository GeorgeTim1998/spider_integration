import numpy as np
from math import pi
import support as sup
import matplotlib.pyplot as pyplot
from fenics import *
from scipy.interpolate import interp2d, LinearNDInterpolator
from fenics_support import countour_plot_via_mesh, assign_const_to_nan_in_expression, plot_1D, interpolate_spider_data_on_function_space
from helpers import eqdsk_equlx_helper as eq

#%% Prescript things
path = '/media/george/part/Spider'
working_folder = 'WK_PRIME'
pic_path = "Pics/%s" % working_folder
wr_file = 'spik.wr'

path_to_file = "%s/%s/%s" % (path, working_folder, wr_file)

#%% Get data from file
with open(path_to_file, 'r') as file:
  psi_size, spacial_size, size, psi_max = sup.first_line(file)
  psi_min, end_entry = sup.second_line(file)

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

sqrt_psi_norm, data = sup.return_and_delete_range(data, psi_size)

dpdpsi, data = sup.return_and_delete_range(data, psi_size) 
dfdpsi, data = sup.return_and_delete_range(data, psi_size) 

rc, data = sup.return_and_delete_range(data, spacial_size) 
zc, data = sup.return_and_delete_range(data, spacial_size) 
rb, data = sup.return_and_delete_range(data, spacial_size) 
zb, data = sup.return_and_delete_range(data, spacial_size) 

ro, data = sup.return_and_delete_range(data, psi_size*spacial_size) 

q, data = sup.return_and_delete_range(data, psi_size)
fvac, data = sup.return_and_delete_range(data, 1)
fvac = fvac[0]

ro = ro.reshape(psi_size, spacial_size)

#%% Plot data needed from Spider
psi = psi_max * (1 - sqrt_psi_norm**2) # This is magnetic flux/2pi. In spider flux is used/ Multiply by 2pi?

# ppsi = sup.restore_funcpsi(2*pi*psi, dpdpsi)
# fpsi = sup.restore_funcpsi(2*pi*psi, dfdpsi, fvac)

dpdpsi = sup.restore_dpdpsi(dpdpsi)
dfdpsi = sup.restore_dfdpsi(dfdpsi)


ppsi, fpsi = sup.restore_pres_n_fpol(psi.max(), 0, len(psi), dpdpsi, dfdpsi, fvac)

I = np.ones((psi_size, spacial_size))
r_mesh = ro*(rb - rc) + I*rc
z_mesh = ro*(zb - zc) + I*zc

psi_mesh = (I.transpose() * psi).transpose()
p_mesh = (I.transpose() * ppsi).transpose()
f_mesh = (I.transpose() * fpsi).transpose()
dfdpsi_mesh = (I.transpose() * dfdpsi).transpose()
dpdpsi_mesh = (I.transpose() * dpdpsi).transpose()
q_mesh = (I.transpose() * q).transpose()

# sup.countour_plot_maxtrix(r_mesh, z_mesh, psi_mesh, 20, grid=True, colorbar=True, note='Spider', PATH=pic_path, plot_title='Spider')

#%% try to find what to do...
pprime = eq.default_pprime()
ffprim = eq.default_ffprim()

UM = eq.default_UM()
UP = eq.default_UP()
FPOLB = eq.default_fpolb()

pres, fpol = eq.restore_pres_n_fpol(UM, UP, len(pprime), pprime, ffprim, FPOLB**0.5, FPOLB**0.5)

# plot_1D(np.linspace(0, 1, len(fpol)), pres, xlabel='coord', ylabel='pres', PATH=pic_path)
# plot_1D(r_mesh[:, 0], p_mesh[:, 0], xlabel='coord', ylabel='p_mesh', PATH=pic_path)
plot_1D(np.linspace(0, 1, len(fpol)), fpol, xlabel='coord', ylabel='fpol', PATH=pic_path)
plot_1D(r_mesh[:, 0], f_mesh[:, 0], xlabel='coord', ylabel='f_mesh', PATH=pic_path)
# plot_1D(r_mesh[:, 0], dpdpsi_mesh[:, 0], xlabel='coord', ylabel='dpdpsi', PATH=pic_path)
# plot_1D(r_mesh[:, 0], dfdpsi_mesh[:, 0], xlabel='coord', ylabel='f*dfdpsi', PATH=pic_path)
