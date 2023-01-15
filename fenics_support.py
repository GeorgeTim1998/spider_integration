from fenics import *
import time
import datetime
import numpy
import matplotlib.tri as tri
import matplotlib.pyplot as matplt

DPI = 240
PICS_FOLDER = 'Pics'

def operator_weights(V):
  r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
  r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such
  
  return r_2, r

def linear_profiles():
  p0 = 1
  f0_2 = 1
  psi0 = 1
  
  return p0, f0_2, psi0

def Time_name():
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime)  # get current time to make figure name unique
    return time_title
  
def countour_plot_via_mesh(gmsh, u, levels,
                           plot_title='',
                           xticks_array = [],
                           yticks_array = [],
                           grid = False,
                           colorbar = False):

  u_min = u.vector()[:].min()
  u_max = u.vector()[:].max()
  
  if u_min == u_max:
    print("Trivial solution. u = %s" % u_max)
  else:
    triang = tri.Triangulation(*gmsh.coordinates().reshape((-1, 2)).T, triangles=gmsh.cells())
    u_array = u.compute_vertex_values(gmsh)

    fig = matplt.tricontour(triang, u_array, levels, zorder=2)
    matplt.gca().set_aspect("equal")
    
    # matplt.xlim(gmsh.plot_domain[0], gmsh.plot_domain[1])
    # matplt.ylim(gmsh.plot_domain[2], gmsh.plot_domain[3])
    matplt.xlabel("r, см")
    matplt.ylabel("z, см")
    
    if xticks_array != []:
        matplt.xticks(numpy.array(xticks_array))
    if yticks_array != []:
        matplt.xticks(numpy.array(yticks_array))

    if grid == True:
        matplt.grid(True)
    
    if colorbar == True:
        matplt.colorbar(fig).set_label("\u03C8(r, z), Вб (СГС)")

    save_contour_plot(plot_title)

  return u_max

def save_contour_plot(plot_title):
  time_title = Time_name()

  path_my_file = '%s/%s' % (PICS_FOLDER, time_title)
  file_path = "%s.png" % path_my_file
  
  matplt.title(plot_title)
  matplt.savefig(file_path, dpi=DPI, bbox_inches="tight")
  matplt.close()

  print("3D countour plot saved to PATH: %s" % file_path)
  time.sleep(1)
  
