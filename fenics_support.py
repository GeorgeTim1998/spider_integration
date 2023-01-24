from fenics import *
import time
import datetime
import numpy
import matplotlib.tri as tri
import matplotlib.pyplot as matplt
import sympy

DPI = 240
PICS_FOLDER = 'Pics'

M0 = 1.25e-6

def operator_weights(V):
  r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
  r = interpolate(Expression('x[0]', degree = 1), V) # interpolation is needed so that 'a' could evaluate deriviations and such
  z = interpolate(Expression('x[1]', degree = 1), V) # interpolation is needed so that 'a' could evaluate deriviations and such
  
  return r_2, r, z

def linear_profiles():
  bp = 0.9
  f0_2 = 1
  psi0 = 1
  
  return bp, f0_2, psi0

def Time_name():
    ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
    time_title = str(ttime)  # get current time to make figure name unique
    return time_title

def point_source(I, sigma, R0):
  x = sympy.symbols('x[0]')
  z = sympy.symbols('x[1]')

  j0 = I/(pi*sigma**2)

  f_expr = j0 * sympy.exp(-((x-R0)**2 + z**2) / sigma**2)

  f_text = sympy.printing.ccode(f_expr)
  f_text = f_text.replace('exp', 'std::exp')

  return Expression(f_text, degree=2)
  
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
    gmsh_coordinates = gmsh.coordinates().reshape((-1, 2)).T
    triang = tri.Triangulation(*gmsh_coordinates, triangles=gmsh.cells())
    u_array = u.compute_vertex_values(gmsh)

    fig = matplt.tricontour(triang, u_array, levels, zorder=2)
    matplt.gca().set_aspect("equal")
    
    matplt.xlim(gmsh_coordinates[0].min(), gmsh_coordinates[0].max())
    matplt.ylim(gmsh_coordinates[1].min(), gmsh_coordinates[1].max())
    matplt.xlabel("r, м")
    matplt.ylabel("z, м")
    
    if xticks_array != []:
        matplt.xticks(numpy.array(xticks_array))
    if yticks_array != []:
        matplt.xticks(numpy.array(yticks_array))

    if grid == True:
        matplt.grid(True)
    
    if colorbar == True:
        matplt.colorbar(fig).set_label("\u03C8(r, z), Вб")

    save_contour_plot(plot_title)

  return u_max

def save_contour_plot(plot_title = ''):
  time_title = Time_name()

  path_my_file = '%s/%s' % (PICS_FOLDER, time_title)
  file_path = "%s.png" % path_my_file
  
  matplt.title(plot_title)
  matplt.savefig(file_path, dpi=DPI, bbox_inches="tight")
  matplt.close()

  print("3D countour plot saved to PATH: %s" % file_path)
  time.sleep(1)

def boundary_length(ds): 
  # checked against circle tor
  return assemble(Constant(1) * ds)

def form_vector_space(u):
  V = u.function_space()
  gmsh = V.mesh()
  degree = V.ufl_element().degree()
  W = VectorFunctionSpace(gmsh, 'P', degree)
  
  return W

def calculate_Bp(u, r, W):
  # this works correctly for u = r**2
  # I get [0, 2]
  grad_u = project(grad(u), W)

  return project(as_vector( (-1/r * grad_u[1], 1/r * grad_u[0]) ), W)

def circulation(vector, n, ds):
  # checked for er. Circ = 0. OK
  # checked for z*er. Circ = pi*a**2. OK
  # checked for point source. Problem statement was false. Missing function. True
  
  tangent = as_vector([n[1], -n[0]])
  scalar_product = dot(vector, tangent)
  
  return assemble(scalar_product * ds)

def calculate_orts(W):
  er = Expression(("1", "0"), degree = 1)
  ez = Expression(("0", "1"), degree = 1)
  
  return interpolate(er, W), interpolate(ez, W)

def calculate_omega(r, gmsh): 
  # checked against circle tor
  return assemble(2 * numpy.pi * r * dx(gmsh))

def calculate_plasma_cross_surface(gmsh): 
  # checked against circle tor
  return assemble(Constant(1) * dx(gmsh))

def calculate_plasma_surface(r, ds): 
  # checked against circle tor
  return assemble(2*pi * r*ds)

def return_R0(u, V): # multiply by 1.1 and S1 = 2
  max_index = numpy.argmax(u.vector()[:])
  max_point = V.tabulate_dof_coordinates() # can be used for assigning?? from known data
  
  return max_point[max_index][0]

def return_q(r, z, R0, V, W):
  # checked all q vectors are ok!
  q = []
  R0 = interpolate(Constant(R0), V)
  zero = interpolate(Constant(0), V)
  
  q.append( project(as_vector( (r - R0, z) ), W) )
  q.append( project(as_vector( (R0, zero) ), W) )
  q.append( project(as_vector( (zero, z) ), W) )
  
  return q
  
def calculate_S_integrals(Bpa, omega, Bp, q, n, r, ds):
  S = {'S1': 0, 'S2': 0, 'S3': 0}
  for i, qi in enumerate(q):
    Si = "S%d" % (i + 1)
    S[Si] = 1 / Bpa**2 / omega * assemble( dot(Bp, Bp) * dot(qi, n) * 2*pi*r*ds )
    
  return S

def solve_SLAE(alpha, S, Rt, R0):
  matrix = numpy.array(
    [
      [3,            1, -1],
      [1,            1,  1],
      [1, -(alpha - 1), -1]
    ]
  )
  
  f1 = S['S1'] + S['S2']
  f2 = Rt/R0 * S['S2']
  f3 = S['S3']
  
  f = numpy.array([f1, f2, f3])
  
  x = numpy.linalg.solve(matrix, f)
  
  print(numpy.allclose(numpy.dot(matrix, x), f))

  return x

def form_dict():
  problem_dict = {}
  
  problem_dict['omega'] = []
  problem_dict['S_'] = []
  problem_dict['Spl'] = []
  problem_dict['Rt'] = []
  problem_dict['alpha_LB'] = []
  problem_dict['R0'] = []
  problem_dict['S1'] = []
  problem_dict['S2'] = []
  problem_dict['S3'] = []
  problem_dict['bp'] = []
  problem_dict['li'] = []
  problem_dict['mu_i'] = []
  
  return problem_dict


def plot_1D(x, y, xlabel='', ylabel=''):

  matplt.scatter(x, y)

  matplt.grid("True")
  
  if xlabel != '':
    matplt.xlabel(xlabel)
  
  if ylabel != '':
    matplt.ylabel(ylabel)

  save_contour_plot()
    
def acceptable_value(plasma_vals):
  answer = True
  for value in plasma_vals:
    if abs(value) >= 100:
      answer = False
      
  return answer