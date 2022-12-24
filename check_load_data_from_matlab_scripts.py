from load_data_from_matlab_scripts import LoadDataFromMatlabScripts
import matplotlib.pyplot as matplt
import numpy as np
import math as mh

def make_arrays_ready_for_tricontour(matlab_data):
  matlab_data.r = matlab_data.r.repeat(matlab_data.psi.shape[0])
  matlab_data.z = np.tile(matlab_data.z, matlab_data.psi.shape[1])
  matlab_data.psi = matlab_data.psi.transpose().reshape(mh.prod(matlab_data.psi.shape))
      

matlab_data = LoadDataFromMatlabScripts()
make_arrays_ready_for_tricontour(matlab_data)

matplt.gca().set_aspect("equal")
matplt.tricontour(matlab_data.r, matlab_data.z, matlab_data.psi, 100)
matplt.show()