from load_data_from_matlab_scripts import LoadDataFromMatlabScripts
import matplotlib.pyplot as matplt
import numpy as np
import math as mh
import datetime
import time
import fenics as f
from scipy.interpolate import interp2d

matlab_data = LoadDataFromMatlabScripts()

#%% create mesh
min_point = f.Point(matlab_data.r.min(), matlab_data.z.min())
max_point = f.Point(matlab_data.r.max(), matlab_data.z.max())
mesh = f.RectangleMesh(min_point, max_point, len(matlab_data.z) - 1, len(matlab_data.r) - 1)

#%% create psi
V   = f.FunctionSpace(mesh, "Lagrange", 1)
psi = f.Function(V)

r = matlab_data.r
z = matlab_data.z

dof_indexes = f.vertex_to_dof_map(V) # this is probably what you need. it is an array of indexes of vertices. you need to do it like in your textbook!

interpolant = interp2d(r, z, matlab_data.psi, kind='linear', copy=False, bounds_error=True)

values = []
# for dof_index in dof_indexes:
#   coord = mesh.coordinates()[dof_index]
#   values.append(interpolant(coord[0], coord[1]))
for coord in mesh.coordinates():
  values.append(interpolant(coord[0], coord[1]))
  
psi.vector().set_local(values)
f.plot(psi)
matplt.show()
