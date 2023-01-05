from load_data_from_matlab_scripts import LoadDataFromMatlabScripts
import matplotlib.pyplot as matplt
import numpy as np
import math as mh
import datetime
import time
import fenics as f

def make_arrays_ready_for_tricontour(matlab_data):
  matlab_data.r = matlab_data.r.repeat(matlab_data.psi.shape[0])
  matlab_data.z = np.tile(matlab_data.z, matlab_data.psi.shape[1])
  matlab_data.psi = matlab_data.psi.transpose().reshape(mh.prod(matlab_data.psi.shape))

def Time_name():
  ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
  time_title = str(ttime)  # get current time to make figure name unique
  return time_title
        
def save_contour_plot(PATH):
  time_title = Time_name()

  path_my_file = 'Pics/%s/%s' % (PATH, time_title)
  file_path = "%s.png" % path_my_file
  
  matplt.savefig(file_path, dpi=200, bbox_inches="tight")
  matplt.close()

  print(file_path)
  time.sleep(1)

matlab_data = LoadDataFromMatlabScripts()

matplt.show()

# create mesh
min_point = f.Point(matlab_data.r.min(), matlab_data.z.min())
max_point = f.Point(matlab_data.r.max(), matlab_data.z.max())
mesh = f.RectangleMesh(min_point, max_point, len(matlab_data.z) - 1, len(matlab_data.r) - 1)
# mesh = f.RectangleMesh(min_point, max_point, len(matlab_data.r) - 1, len(matlab_data.z) - 1)
#

# create psi
V   = f.FunctionSpace(mesh, "Lagrange", 1)
psi = f.Function(V)

# psi.vector()[:] = matlab_data.psi.reshape(mh.prod(matlab_data.psi.shape))[f.vertex_to_dof_map(V)]
# psi.vector().set_local(matlab_data.psi.reshape(mh.prod(matlab_data.psi.shape))[f.vertex_to_dof_map(V)])
psi.vector().set_local(matlab_data.psi.reshape(mh.prod(matlab_data.psi.shape)))
# f.vertex_to_dof_map(V) this is probably what you need. it is an array of indexes of vertices. you need to do it like in your textbook!
f.plot(psi)
matplt.show()
#

# plot matlab psi data
make_arrays_ready_for_tricontour(matlab_data)

# matplt.gca().set_aspect("equal")
# matplt.tricontour(matlab_data.r, matlab_data.z, matlab_data.psi, 100)
# save_contour_plot('matlab_data_plots')
#