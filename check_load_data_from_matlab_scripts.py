from load_data_from_matlab_scripts import LoadDataFromMatlabScripts
import matplotlib.pyplot as matplt
import numpy as np
import math as mh
import datetime
import time

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
make_arrays_ready_for_tricontour(matlab_data)

matplt.gca().set_aspect("equal")
matplt.tricontour(matlab_data.r, matlab_data.z, matlab_data.psi, 100)
save_contour_plot('matlab_data_plots')