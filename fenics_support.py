from fenics import *
import time
import datetime
import numpy
import matplotlib.tri as tri
import matplotlib.pyplot as matplt
import sympy
from termcolor import colored, cprint

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

def linear_pressure(bp, Re):
  # Re ellipse center which sometimes can be R0 in your writings
  return 2*pi**2 * bp/Re**4

def linear_tor_function(q, E, Re):
  return 2*pi**2 * E**2 * q**2 / Re**2

def inverced_r_integral(Re, a, b, r, dx, gmsh):
  return Re/(pi*a*b) * assemble(1/r * dx(gmsh))

def measure_u(Re, a, b, I, bp, alph, E, q, u, V):
  psi0 = Re**3/(2*pi**3 * a*b) * M0*I/(bp + alph*E**2 * q**2)
  
  return project(psi0 * u, V), psi0

def calculate_p_psi(bp, psi0, u, Re):
  V = u.function_space()
  p0 = 2*pi**2 * psi0**2 * bp / (M0 * Re**4)
  
  return project(p0 * u/psi0, V), p0

def calculate_Fpow2_psi(E, psi0, q, u, Re):
  V = u.function_space()
  F_20 = 4*pi**2 * E**2 * psi0**2 * q**2 / Re**2
  
  return project(F_20 * u/psi0, V), F_20

def calculate_Bt(F_2_psi, r):
  V = F_2_psi.function_space()
  
  return project(sqrt(F_2_psi)/r, V)
  

def calculate_g(p_psi, Bp, Bt, Bt0 = 0):
  V = p_psi.function_space()
  
  if Bt0 == 0:
    return project(2*M0 * p_psi + dot(Bp, Bp) - dot(Bt, Bt), V)

def calculate_Rt(g, r, dx, gmsh):
  upper = assemble(g * dx(gmsh)) 
  downer = assemble(g/r * dx(gmsh)) 
  
  return upper/downer

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
                           colorbar = False,
                           note='3D countour plot saved to PATH:'):

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

    save_contour_plot(note, plot_title)

  return u_max

def calculate_d_at_boundary(u, psi0):
  gmsh = u.function_space().mesh()
  gmsh_coordinates = gmsh.coordinates().reshape((-1, 2)).T
  triang = tri.Triangulation(*gmsh_coordinates, triangles=gmsh.cells())
  u_array = u.compute_vertex_values(gmsh)

  fig = matplt.tricontour(triang, u_array, [0, psi0], zorder=2)
  
  contour = fig.allsegs[-1][0].transpose() # last contour w/ psi0 basically
  a_psi_level = 0.5 * (contour[0].max() - contour[0].min())
  b_psi_level = 0.5 * (contour[1].max() - contour[1].min())
  E_psi_level = b_psi_level/a_psi_level
  
  zero_contour = gmsh.coordinates().transpose()
  a = 0.5*(zero_contour[0].max() - zero_contour[0].min())
  b = 0.5*(zero_contour[1].max() - zero_contour[1].min())
  E = b/a
  
  matplt.close()
  der_E = (E - E_psi_level)/(a - a_psi_level)
  return a * der_E / (2*E)

def save_contour_plot(note, plot_title = ''):
  time_title = Time_name()

  path_my_file = '%s/%s' % (PICS_FOLDER, time_title)
  file_path = "%s.png" % path_my_file
  
  matplt.title(plot_title)
  matplt.savefig(file_path, dpi=DPI, bbox_inches="tight")
  matplt.close()

  print_colored(note, color='green', white_str=file_path)
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

def calculate_alpha(Bp, ez, gmsh, dx):
  numerator = 2 * assemble(dot(Bp, ez) * dot(Bp, ez) * dx(gmsh))
  denominator = assemble(dot(Bp, Bp) * dx(gmsh))
  
  return numerator/denominator

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
    
  return S['S1'], S['S2'], S['S3']

def solve_SLAE(alpha, S, Rt, R0):
  matrix = numpy.array(
    [
      [3,            1, -1],
      [1,            1,  1],
      [1, -(alpha - 1), -1]
    ]
  )
  
  f1 = S[0] + S[1]
  f2 = Rt/R0 * S[1]
  f3 = S[2]
  
  f = numpy.array([f1, f2, f3])
  
  x = numpy.linalg.solve(matrix, f)
  
  if numpy.allclose(numpy.dot(matrix, x), f):
    print_colored('SLAE Solution is:', 'green', 'True', attrs=['bold'])
  else:
    print_colored('SLAE Solution is:', 'red', 'True', attrs=['bold'])

  return x

def calculate_bp_theory(Bpa, omega, p_psi, dx, gmsh, r):
  # this gives simelar results that SLAE does
  return 2*M0/Bpa**2/omega * assemble(p_psi * 2*pi*r*dx(gmsh))

def calculate_li_theory(Bpa, omega, Bp, dx, gmsh, r):
  # this gives simelar results that SLAE does. But 6% off
  return 1/Bpa**2/omega * assemble(dot(Bp, Bp) * 2*pi*r*dx(gmsh))

def calculate_mu_i_theory(Bpa, omega, Bt, dx, gmsh, r, Bt0 = 0):
  # this gives simelar results that SLAE does
  if Bt0 == 0:
    return -1/Bpa**2/omega * assemble(dot(Bt, Bt) * 2*pi*r*dx(gmsh))

def append_problem_data(vars, problem_data):
  keys = list(problem_data.keys())
  for key in keys:
    problem_data[key].append(vars[key])
    
  return problem_data

def form_dict():
  problem_dict = {}
  
  problem_dict['omega'] = []
  problem_dict['S_'] = []
  problem_dict['Spl'] = []
  
  
  problem_dict['alpha'] = []
  problem_dict['eps_K'] = []
  problem_dict['Rt'] = []
  
  problem_dict['S1'] = []
  problem_dict['S2'] = []
  problem_dict['S3'] = []
  
  problem_dict['bp'] = []
  problem_dict['li'] = []
  problem_dict['mu_i'] = []
  
  return problem_dict

def form_dict_additions():
  problem_dict = {}
  
  problem_dict['alpha_LB'] = []
  
  problem_dict['R0'] = []
  
  problem_dict['S1_theory'] = []
  problem_dict['S2_theory'] = []
  problem_dict['S3_theory'] = []
  
  problem_dict['bp_theory'] = []
  problem_dict['li_theory'] = []
  problem_dict['mu_i_theory'] = []
  
  return problem_dict

def addition_keys():
  addition = {}
  
  addition['alpha'] = 'alpha_LB'
  
  addition['Rt'] = 'R0'
  
  addition['S1'] = 'S1_theory'
  addition['S2'] = 'S2_theory'
  addition['S3'] = 'S3_theory'
  
  addition['bp'] = 'bp_theory'
  addition['li'] = 'li_theory'
  addition['mu_i'] = 'mu_i_theory'
  
  return addition

def plot_1D(x, y, xlabel='', ylabel='', note='', additions=[]):

  matplt.scatter(x, y)
  matplt.legend([note], loc='best')
  if additions != []:
    matplt.scatter(x, additions)
    matplt.legend([note, addition_keys()[note]], loc='best')

  matplt.grid("True")
  
  if xlabel != '':
    matplt.xlabel(xlabel)
  
  if ylabel != '':
    matplt.ylabel(ylabel)

  save_contour_plot(note="2D plot of %s saved to PATH:" % note)

def acceptable_value(plasma_vals):
  answer = True
  for value in plasma_vals:
    if abs(value) >= 100:
      answer = False
      
  return answer

def print_colored(color_srt, color='white', white_str='', attrs=[]):
  print(colored(color_srt, color, attrs=attrs), white_str)