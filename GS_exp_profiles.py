from fenics import *
import matplotlib.pyplot as pyplot
import support as sup
import fenics_support as fsup
from math import pi
import math as mt
import numpy as np

from pathlib import Path
my_dir = "Pics/%s" % fsup.Time_name()
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
Bt0_mean = lao_hash['Bt0_mean']
p0 = lao_hash['p0'] 
psi0 = lao_hash['psi0']
F2_0 = lao_hash['F2_0']

psi_level = 1e-2 # used to calc contours with desired level
Fpl_vs_Fvac_ratio = lao_hash['Fpl_vs_Fvac_ratio']

al_pr = 0.4
kappa = 1e-1 
divider = 1 - mt.exp(al_pr) + al_pr
btor = 4e-2
bpol = 0.9

sup.create_readmi(lao_hash, my_dir)

#%% Dict that will help me plot calculated data
problem_data = fsup.form_dict()
theory_data = fsup.form_dict_additions()
plot_keys = fsup.addition_keys() # do later smth with it

#%% Files
filenames = ['a_0.370_ratio_1.000_msh_3.0e-03',
             'a_0.370_ratio_1.100_msh_3.0e-03',
             'a_0.370_ratio_1.200_msh_1.0e-02',
             'a_0.370_ratio_1.300_msh_1.0e-02',
             'a_0.370_ratio_1.400_msh_1.0e-02',
             'a_0.370_ratio_1.500_msh_1.0e-02',
             'a_0.370_ratio_1.600_msh_1.0e-02',
             'a_0.370_ratio_1.700_msh_1.0e-02',
             'a_0.370_ratio_1.800_msh_1.0e-02',
             'a_0.370_ratio_1.900_msh_1.0e-02',
             'a_0.370_ratio_2.000_msh_1.0e-02']
# 1.7 is behaving strangely. maybe its integral from g
# filenames = ['a_37.000_ratio_1.100_msh_1.0e+00']

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
  ell_b = ell_a * E[i]
  [r2, r, z] = fsup.operator_weights(V)
  
  # u = Function(V)
  u = project(Expression("1 - pow(x[0]-Re, 2)/pow(ell_a, 2) - pow(x[1], 2)/pow(ell_b, 2)", degree=2, Re=Re, ell_a=ell_a, ell_b=ell_b), V)
  v = TestFunction(V)
  
  f = Expression("2*pow(pi,2) * al_pr * bpol/pow(Re,4) * ( pow(x[0],2) + kappa/btor * pow(Re,2) ) / divider", degree=2, al_pr=al_pr, bpol=bpol, Re=Re, kappa=kappa, btor=btor, divider=divider)
  F = inner(grad(u)/r, grad(r2*v))*dx + f*(exp(al_pr*u) - 1) * r*v*dx
  
#%% Compute solution and p(psi), F(psi), J(psi)
  solve(F == 0, u, bc)

  fsup.countour_plot_via_mesh(gmsh, u, levels = 20, colorbar=True, grid=True, PATH=my_dir, plot_title="E = %.1f" % E[i])
  exit()
  d = fsup.calculate_d_at_boundary(u, psi_level)
  fsup.print_colored("d", 'green', d)
  print("\n")
  
  p_psi = fsup.calculate_p_psi_final(psi0, u, p0)
  F2_psi = fsup.calculate_F2_pow_1(F2_0, psi0, u, Fpl_vs_Fvac_ratio)
  [J_psi, I] = fsup.calculate_J_psi_final(p0, psi0, F2_0, r, V, dx, Fpl_vs_Fvac_ratio)

#%% Post solve calculus
  boundaries = MeshFunction('size_t', gmsh, gmsh.topology().dim() - 1) # get boundaries (and all marked lines???) from mesh
  ds = Measure('ds', domain=gmsh, subdomain_data=boundaries)
  n = FacetNormal(gmsh) # normal to plasma boundary

  L = fsup.boundary_length(ds)
  W = fsup.form_vector_space(u)
  Bp = fsup.calculate_Bp(u, r, W)
  Bpa = 1/L * fsup.circulation(Bp, n, ds)

  Bt = fsup.calculate_Bt(F2_psi, r)
  Bt0 = fsup.calculate_Bt0(F2_0, r, V)
  g_sol = fsup.calculate_g(p_psi, Bp, Bt, Bt0=Bt0)
  Rt = fsup.calculate_Rt(g_sol, r, dx, gmsh)

  [er, ez] = fsup.calculate_orts(W)

  omega = fsup.calculate_omega(r, gmsh)
  S_ = fsup.calculate_plasma_cross_surface(gmsh)
  Spl = fsup.calculate_plasma_surface(r, ds)
  alpha = fsup.calculate_alpha(Bp, ez, gmsh, dx)
  alpha_LB = 2 * E[i]**2 / (E[i]**2 + 1)
  eps_K = (E[i]**2 - 1) / (E[i]**2 + 1)
  R0 = fsup.return_R0(u, V)

  fsup.print_colored("\u03A9 =", 'blue', omega)
  fsup.print_colored('Spl =', 'blue', Spl)
  fsup.print_colored("S\u27C2 =", 'blue', S_)
  fsup.print_colored('L =', 'blue', L)
  fsup.print_colored("\u03B1 =", 'blue', alpha)
  fsup.print_colored("\u03B1_LB =", 'blue', alpha_LB)
  print("\n")
  
  q = fsup.return_q(r, z, R0, V, W)
  S1, S2, S3 = fsup.calculate_S_integrals(Bpa, omega, Bp, q, n, r, ds)
  
  fsup.print_colored('S1 =', 'blue', S1)
  fsup.print_colored('S2 =', 'blue', S2)
  fsup.print_colored('S3 =', 'blue', S3)
  print("\n")

#%% Plot plasma profiles
  fsup.plot_big_axis_profile(p_psi, yaxis='Pressure, Pa', grid=True, note='2D plot of p(psi) saved to PATH:', PATH=my_dir)
  fsup.plot_big_axis_profile(F2_psi, yaxis='Tor func, m*Tl', grid=True, note='2D plot of F2(psi) saved to PATH:', PATH=my_dir)
  fsup.plot_big_axis_profile(Bt, yaxis='Tor field, Tl', grid=True, note='2D plot of Bt(psi) saved to PATH:', PATH=my_dir)
  fsup.plot_big_axis_profile(J_psi, yaxis='Current Density, A/m**2', grid=True, note='2D plot of J(psi) saved to PATH:', PATH=my_dir)
  print("\n")
  
#%% Calc magnetic values
  [r_v, q_v] = fsup.return_q_1D(R0, ell_a, u, Bt, Bp)
  fsup.plot_1D(r_v, q_v, xlabel='Major radius', ylabel='q', note='q', PATH=my_dir)

  bp, li, mu_i = fsup.solve_SLAE(alpha, [S1, S2, S3], Rt, R0)
  fsup.print_colored('SLAE', 'blue', attrs=['bold'])
  fsup.print_colored('bp, li, mui', 'blue', [bp, li, mu_i])

  bp_theory = fsup.calculate_bp_theory(Bpa, omega, p_psi, dx, gmsh, r)
  li_theory = fsup.calculate_li_theory(Bpa, omega, Bp, dx, gmsh, r)
  mu_i_theory = fsup.calculate_mu_i_theory(Bpa, omega, Bt, dx, gmsh, r, Bt0=Bt0)
  
  fsup.print_colored('Theory', 'blue', attrs=['bold'])
  fsup.print_colored('bp, li, mui', 'blue', [bp_theory, li_theory, mu_i_theory])
  fsup.print_colored('-'*50, 'yellow', attrs=['bold'])
  print("\n")
  
#%% Calculate S1-S3 based on known bp, li, mu_i  
  S1_theory = 2
  S2_theory = 2*bp_theory + li_theory - 1
  S3_theory = 1 - eps_K/2 - d*(1 - eps_K**2/2)

#%% Append data
  problem_data = fsup.append_problem_data(globals(), problem_data)
  theory_data = fsup.append_problem_data(globals(), theory_data)

#%% Post problem plot
fsup.print_colored('Save 2D plots...', 'green', attrs=['bold'])
keys = list(problem_data.keys())
for key in keys:
  if key in list(plot_keys.keys()):
    fsup.plot_1D(E, problem_data[key],
                xlabel='elongation', ylabel=key, note=key,
                additions=theory_data[plot_keys[key]],
                PATH=my_dir) # maybe add special funcs for certain values
  else:
    fsup.plot_1D(E, problem_data[key],
                xlabel='elongation', ylabel=key, 
                note=key,
                PATH=my_dir) # maybe add special funcs for certain values
