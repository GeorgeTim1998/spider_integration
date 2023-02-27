import numpy as np
from math import pi
import support as sup
import matplotlib.pyplot as pyplot
from fenics import *
from scipy.interpolate import interp2d, LinearNDInterpolator
from fenics_support import countour_plot_via_mesh, assign_const_to_nan_in_expression, plot_1D, interpolate_spider_data_on_function_space, M0
from helpers import eqdsk_equlx_helper as eq

#%% Prescript things
path = '/media/george/part/Spider'
working_folder = 'WK_linear_profs_no2pi_vs_fenics_slimmest'
pic_path = "Pics/%s" % working_folder
wr_file = 'spik.wr'

path_to_file = "%s/%s/%s" % (path, working_folder, wr_file)

sup.print_colored("Launch Spider program", 'red', "\n", ["bold"])

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

dpdpsi = sup.restore_dpdpsi(dpdpsi)
dfdpsi = sup.restore_dfdpsi(dfdpsi)

# ppsi = sup.restore_funcpsi(2*pi*psi, dpdpsi)
# fpsi = sup.restore_funcpsi(2*pi*psi, dfdpsi, fvac)

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

Jpsi = r_mesh[:, 0]*dpdpsi + dfdpsi/(M0*r_mesh[:, 0])

sup.print_colored("\u03C8 max =", color='green', white_str=psi_mesh.max())
sup.print_colored("p' =", color='green', white_str=dpdpsi.max())
sup.print_colored("F2' =", color='green', white_str=dfdpsi.max())
print("\n")

sup.countour_plot_maxtrix(r_mesh, z_mesh, psi_mesh, levels=[0.001, 0.002, 0.003, 0.01], grid=True, colorbar=True, note='3D countour plot from Spider saved to PATH:', PATH=pic_path, plot_title='Spider')

exit()
plot_1D(r_mesh[:, 0], dpdpsi, xlabel='coord', ylabel='dpdpsi. Spider', PATH=pic_path)
# plot_1D(r_mesh[:, 0], p_mesh[:, 0], xlabel='coord', ylabel='p_mesh. Spider', PATH=pic_path)
# plot_1D(r_mesh[:, 0], f_mesh[:, 0], xlabel='coord', ylabel='f_mesh. Spider', PATH=pic_path)
# plot_1D(r_mesh[:, 0], Jpsi, xlabel='coord', ylabel='Jpsi. Spider', PATH=pic_path)
#%% try to find what to do...
plot_1D(r_mesh[:, 0], p_mesh[:, 0], xlabel='coord', ylabel='p_mesh', PATH=pic_path)
plot_1D(r_mesh[:, 0], multiplicator_p*dpdpsi_mesh[:, 0], xlabel='coord', ylabel='dpdpsi', PATH=pic_path)
plot_1D(r_mesh[:, 0], f_mesh[:, 0], xlabel='coord', ylabel='f_mesh', PATH=pic_path)
plot_1D(r_mesh[:, 0], multiplicator_f*f_mesh[:, 0]*dfdpsi_mesh[:, 0], xlabel='coord', ylabel='f*dfdpsi', PATH=pic_path)

exit()
#%% Import Spider solution to fenics
filename = 'test'
folder = sup.xml_files_folder()
xml_file = "%s/%s.xml" % (folder, filename)

gmsh = Mesh(xml_file)
V = FunctionSpace(gmsh, 'Lagrange', 1)

psi = interpolate_spider_data_on_function_space(r_mesh, z_mesh, psi_mesh, V)
ppsi = interpolate_spider_data_on_function_space(r_mesh, z_mesh, p_mesh, V)
countour_plot_via_mesh(gmsh, psi, levels=50, colorbar=True, grid=True, xlim=[0.8, 2.4], ylim=[-1, 1], PATH=pic_path)
countour_plot_via_mesh(gmsh, ppsi, levels=50, colorbar=True, grid=True, xlim=[0.8, 2.4], ylim=[-1, 1], PATH=pic_path)