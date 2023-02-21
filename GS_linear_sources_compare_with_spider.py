from fenics import *
import matplotlib.pyplot as pyplot
import support as sup
import fenics_support as fsup
from math import pi
import numpy as np

from pathlib import Path
my_dir = "Pics/Compare_with_Spider"
Path(my_dir).mkdir(parents=True, exist_ok=True)

fsup.print_colored("Launch program for linear p(\u03C8) and F\u00b2(\u03C8)", 'red', "\n", ["bold"])
print("\n")

folder = sup.xml_files_folder()

#%% Plasma data
lao_hash = sup.lao_hash()
Re = lao_hash['Re'] # ellipse center which sometimes can be R0 in your writings
E = lao_hash['E'] # ellipse elongation
I = lao_hash['I']
ell_a = lao_hash['a']

bp_problem = 0.9 # poloidal betta from lao1985
q_problem = 1 # stability from lao1985
psi_level = 1e-3 # used to calc contours with known desired psi_level
Fpl_vs_Fvac_ratio = 1e-1

#%% Dict that will help me plot calculated data
problem_data = fsup.form_dict()
theory_data = fsup.form_dict_additions()
plot_keys = fsup.addition_keys() # do later smth with it

#%% Files
# filenames = ['a_0.370_ratio_1.000_msh_3.0e-03',
#              'a_0.370_ratio_1.100_msh_3.0e-03',
#              'a_0.370_ratio_1.200_msh_1.0e-02',
#              'a_0.370_ratio_1.300_msh_1.0e-02',
#              'a_0.370_ratio_1.400_msh_1.0e-02',
#              'a_0.370_ratio_1.500_msh_1.0e-02',
#              'a_0.370_ratio_1.600_msh_1.0e-02',
#              'a_0.370_ratio_1.700_msh_1.0e-02',
#              'a_0.370_ratio_1.800_msh_1.0e-02',
#              'a_0.370_ratio_1.900_msh_1.0e-02',
            #  'a_0.370_ratio_2.000_msh_1.0e-02']
# 1.7 is behaving strangely. maybe its integral from g
filenames = ['a_0.370_ratio_1.000_msh_3.0e-03']

#%% Program body
for i, filename in enumerate(filenames):
  fsup.print_colored("Iteration N0 %d out of %d for file:" % (i, len(filenames)), 'green', filename, attrs=['bold'])
  print("\n")
    
  mesh_size = sup.mesh_size(filename)

  xml_file = "%s/%s.xml" % (folder, filename)
  xml_file_facet = "%s/%s_facet_region.xml" % (folder, filename) # triangle surfaces
  xml_file_physical_region = "%s/%s_physical_region.xml" % (folder, filename) # domains

  gmsh = Mesh(xml_file)
  V = FunctionSpace(gmsh, 'P', 1)
#%% Define boundary condition
  u_D = Expression('0', degree=2)

  def boundary(x, on_boundary):
    return on_boundary

  bc = DirichletBC(V, u_D, boundary)
#%% Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)

  ell_b = ell_a * E[i]
  [r2, r, z] = fsup.operator_weights(V)
  
  p_part = fsup.linear_pressure(bp_problem, Re)
  F_part = fsup.linear_tor_function(q_problem, E[i], Re)
  
  a = dot(grad(u)/r, grad(r2*v))*dx
  f = Constant(p_part) * r2 + Constant(F_part) * Fpl_vs_Fvac_ratio
  L = f * r*v*dx
    
#%% Compute solution and p(psi), F(psi), J(psi)
  u = Function(V)
  solve(a == L, u, bc)

  inverced_r_integral = fsup.inverced_r_integral(Re, ell_a, ell_b, r, dx, gmsh)
  [u, psi0] = fsup.measure_u(Re, ell_a, ell_b, I, bp_problem, inverced_r_integral, E[i], q_problem, u, V, Fpl_vs_Fvac_ratio) # de de-measure solution
  
  fsup.countour_plot_via_mesh(gmsh, u, levels = 20, colorbar=True, grid=True, PATH=my_dir, plot_title="E = %.1f" % E[i])
  d = fsup.calculate_d_at_boundary(u, psi_level)
  fsup.print_colored("d", 'green', d)
  print("\n")
  
  [p_psi, p0] = fsup.calculate_p_psi(bp_problem, psi0, u, Re)
  [F2_psi, F2_0] = fsup.calculate_F2_as_small_add(E[i], psi0, q_problem, u, Re, Fpl_vs_Fvac_ratio)
  [J_psi, I] = fsup.calculate_J_psi(p0, psi0, F2_0, r, V, dx, Fpl_vs_Fvac_ratio)
