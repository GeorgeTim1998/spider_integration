from fenics import *
import matplotlib.pyplot as pyplot
import support as sup
import fenics_support as fsup
from math import pi
import numpy as np

from pathlib import Path
my_dir = "Pics/WK_linear_profs_no2pi_vs_fenics_slimmest"
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

M0 = 1.25e-6
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
  
  p0 = 10000
  psi0 = 0.5
  F2_0 = 0.25

  a = dot(grad(u)/r, grad(r2*v))*dx
  f = (Constant(M0 * p0/psi0) * r2 + Constant(0.5 * F2_0/psi0))
  L = f * r*v*dx
    
#%% Compute solution and p(psi), F(psi), J(psi)
  u = Function(V)
  solve(a == L, u, bc)
  
  ppsi = project(p0 * u/psi0, V)
  F2psi = project(F2_0*(1 + u/psi0), V)
  Fpsi = project(sqrt( F2_0*(1 + u/psi0) ), V)
  
  Jpsi = project(2*pi*(p0/psi0 * r + 1/(2*M0*r) * F2_0/psi0), V)
  I = assemble(Jpsi*dx)
  
  # fsup.plot_big_axis_profile(ppsi, yaxis='p(psi)', PATH=my_dir, grid=True)
  # fsup.plot_big_axis_profile(Fpsi, yaxis='F(psi)', PATH=my_dir, grid=True)
  # fsup.plot_big_axis_profile(Jpsi, yaxis='J(r, psi)', PATH=my_dir, grid=True)
  
  sup.print_colored("I =", color='green', white_str=I)
  sup.print_colored("\u03C8 max =", color='green', white_str=u.vector()[:].max())
  sup.print_colored("p' max =", color='green', white_str=p0/psi0)
  sup.print_colored("F2' max =", color='green', white_str=0.5 * F2_0/psi0)
  print("\n")
  
  fsup.countour_plot_via_mesh(gmsh, u, levels=[0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01], colorbar=True, grid=True, PATH=my_dir, plot_title="Fenics", xlim=[1.5, 1.52], ylim=[-0.05, 0.05])
  print("\n")
  