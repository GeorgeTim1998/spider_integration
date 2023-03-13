import numpy as np
from math import pi
import support as sup
import matplotlib.pyplot as pyplot
from fenics import *
import fenics_support as fsup
from scipy.interpolate import interp2d, LinearNDInterpolator
from fenics_support import countour_plot_via_mesh, assign_const_to_nan_in_expression, plot_1D, interpolate_spider_data_on_function_space, M0
from helpers import eqdsk_equlx_helper as eq

#%% Prescript things
path = '/media/george/part/Spider'
working_folder = 'WK'
pic_path = "Pics/%s" % working_folder
wr_file = 'spik.wr'

path_to_file = "%s/%s/%s" % (path, working_folder, wr_file)

sup.print_colored("Launch Spider program", 'red', "\n", ["bold"])

#%% Pre solve definitions
problem_data = fsup.form_dict()
theory_data = fsup.form_dict_additions()
plot_keys = fsup.addition_keys() # do later smth with it

lao_hash = sup.lao_hash()
ell_a = lao_hash['a']
psi_level = lao_hash['psi_level'] # used to calc contours with desired level
sup.create_readmi(lao_hash, pic_path)

#%% Start calculus
msh_filenames = ["a_0.370_ratio_1.000_msh_1.0e-01",
            "a_0.370_ratio_1.100_msh_1.0e-01",
            "a_0.370_ratio_1.200_msh_1.0e-01",
            "a_0.370_ratio_1.300_msh_1.0e-01",
            "a_0.370_ratio_1.400_msh_1.0e-01",
            "a_0.370_ratio_1.500_msh_1.0e-01",
            "a_0.370_ratio_1.600_msh_1.0e-01",
            "a_0.370_ratio_1.700_msh_1.0e-01",
            "a_0.370_ratio_1.800_msh_1.0e-01",
            "a_0.370_ratio_1.900_msh_1.0e-01",
            "a_0.370_ratio_2.000_msh_1.0e-01"]
E = lao_hash['E']

for i in range(len(msh_filenames)):
  eq.del_eqdsk_equilx(working_folder)

  eq.copy_eqdsk_equilx_to_WK(working_folder, E[i])

  eq.launch_spider_fixed_b(working_folder)
  #%% Get data from file
  with open(path_to_file, 'r') as file:
    psi_size, spacial_size, size, psi_max = sup.first_line(file)
    psi_min, end_entry = sup.second_line(file)

    sqrt_psi_norm = np.zeros(psi_size)
    dpdpsi = np.zeros(psi_size)
    dfdpsi = np.zeros(psi_size)
    
    rc = np.zeros(spacial_size)
    zc = np.zeros(spacial_size)
    
    rb = np.zeros(spacial_size)
    zb = np.zeros(spacial_size)
    
    ro = np.zeros(psi_size*spacial_size)
    
    q = np.zeros(psi_size)
    fvac = 0
    
    data = []
    for line in file: # read rest of lines
      data.append([float(num) for num in line.split()])

  data = np.array(data).reshape(np.size(data))

  sqrt_psi_norm, data = sup.return_and_delete_range(data, psi_size)

  dpdpsi, data = sup.return_and_delete_range(data, psi_size) 
  dfdpsi, data = sup.return_and_delete_range(data, psi_size) 

  rc, data = sup.return_and_delete_range(data, spacial_size) 
  zc, data = sup.return_and_delete_range(data, spacial_size) 
  rb, data = sup.return_and_delete_range(data, spacial_size) 
  zb, data = sup.return_and_delete_range(data, spacial_size) 

  ro, data = sup.return_and_delete_range(data, psi_size*spacial_size) 

  q, data = sup.return_and_delete_range(data, psi_size)
  fvac, data = sup.return_and_delete_range(data, 1)
  fvac = fvac[0]

  ro = ro.reshape(psi_size, spacial_size)

  #%% Plot data needed from Spider
  psi = psi_max * (1 - sqrt_psi_norm**2) # This is magnetic flux/2pi. In spider flux is used/ Multiply by 2pi?

  dpdpsi = sup.restore_dpdpsi(dpdpsi)
  dfdpsi = sup.restore_dfdpsi(dfdpsi)

  # ppsi = sup.restore_funcpsi(2*pi*psi, dpdpsi)
  # fpsi = sup.restore_funcpsi(2*pi*psi, dfdpsi, fvac)

  ppsi, fpsi = sup.restore_pres_n_fpol(psi.max(), 0, len(psi), dpdpsi, dfdpsi, fvac)

  I = np.ones((psi_size, spacial_size))
  r_mesh = ro*(rb - rc) + I*rc
  z_mesh = ro*(zb - zc) + I*zc

  psi_mesh = (I.transpose() * psi).transpose()
  p_mesh = (I.transpose() * ppsi).transpose()
  f_mesh = (I.transpose() * fpsi).transpose()
  dfdpsi_mesh = (I.transpose() * dfdpsi).transpose()
  dpdpsi_mesh = (I.transpose() * dpdpsi).transpose()
  q_mesh = (I.transpose() * q).transpose()

  sup.print_colored("\u03C8 max =", color='green', white_str=psi_mesh.max())
  print("\n")

  sup.countour_plot_maxtrix(r_mesh, z_mesh, psi_mesh, 
                            levels=10, 
                            grid=True, 
                            colorbar=True, 
                            note='3D countour plot from Spider saved to PATH:', 
                            PATH=pic_path, 
                            plot_title="Spider. %.1f" % E[i])

#%% Import Spider solution to fenics
  msh_filename = msh_filenames[i]
  folder = sup.xml_files_folder()
  xml_file = "%s/%s.xml" % (folder, msh_filename)

  gmsh = Mesh(xml_file)
  V = FunctionSpace(gmsh, 'Lagrange', 1)
  [r2, r, z] = fsup.operator_weights(V)

  u = interpolate_spider_data_on_function_space(r_mesh, z_mesh, psi_mesh, V)
  sup.print_colored("Interpolated psi", color='green')
  ppsi_der = interpolate_spider_data_on_function_space(r_mesh, z_mesh, dpdpsi_mesh, V, boundary_val=dpdpsi_mesh[-1, -1])
  sup.print_colored("Interpolated ppsi der", color='green')
  F2psi_der = interpolate_spider_data_on_function_space(r_mesh, z_mesh, dfdpsi_mesh, V, boundary_val=dfdpsi_mesh[-1, -1])
  sup.print_colored("Interpolated F2 der", color='green')

  p_psi = interpolate_spider_data_on_function_space(r_mesh, z_mesh, p_mesh, V)
  sup.print_colored("Interpolated pressure", color='green')
  F2_psi = interpolate_spider_data_on_function_space(r_mesh, z_mesh, f_mesh**2, V, boundary_val=f_mesh[-1, -1]**2)
  sup.print_colored("Interpolated toroidal func", color='green')
  J_psi, I = fsup.define_J_from_spider_data(ppsi_der, F2psi_der, r, V)
  sup.print_colored("Defined current density. I =", color='green', white_str=I)

  countour_plot_via_mesh(gmsh, u, 
                        levels=10, 
                        colorbar=True, 
                        grid=True,
                        plot_title="Fenics. %.1f" % E[i],
                        PATH=pic_path)

#%% Post solve
  d = fsup.calculate_d_at_boundary(u, psi_level)
  fsup.print_colored("d", 'green', d)

  boundaries = MeshFunction('size_t', gmsh, gmsh.topology().dim() - 1) # get boundaries (and all marked lines???) from mesh
  ds = Measure('ds', domain=gmsh, subdomain_data=boundaries)
  n = FacetNormal(gmsh) # normal to plasma boundary

  L = fsup.boundary_length(ds)
  W = fsup.form_vector_space(u)
  Bp = fsup.calculate_Bp(u, r, W)
  Bpa = 1/L * fsup.circulation(Bp, n, ds)
  
  sup.print_colored("Bpa vs M0*I/L =", color='yellow', white_str=Bpa / (M0*I/L))

  Bt = fsup.calculate_Bt(F2_psi, r)
  Bt0 = fsup.calculate_Bt0(fvac**2, r, V)
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
  # Raxis = fsup.calc_big_axis_R(r_mesh, z_mesh)
  fsup.plot_big_axis_profile(p_psi, yaxis='Pressure, Pa', grid=True, note='2D plot of p(psi) saved to PATH:', PATH=pic_path)
  fsup.plot_big_axis_profile(F2_psi, yaxis='Tor func, m*Tl', grid=True, note='2D plot of F2(psi) saved to PATH:', PATH=pic_path)
  fsup.plot_big_axis_profile(Bt, yaxis='Tor field, Tl', grid=True, note='2D plot of Bt(psi) saved to PATH:', PATH=pic_path)
  fsup.plot_big_axis_profile(J_psi, yaxis='Current Density, A/m**2', grid=True, note='2D plot of J(psi) saved to PATH:', PATH=pic_path)
  print("\n")

#%% Calc magnetic values
  [r_v, q_v] = fsup.return_q_1D(R0, ell_a, u, Bt, Bp)
  fsup.plot_1D(r_v, q_v, xlabel='Major radius', ylabel='q', note='q', PATH=pic_path)

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
  # S2_theory = (2*bp_theory + li_theory - 1) * (1 + eps_K/2) 
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
                xlabel="$\it{K}$", ylabel=key, note=key,
                additions=theory_data[plot_keys[key]],
                PATH=pic_path) # maybe add special funcs for certain values
  else:
    fsup.plot_1D(E, problem_data[key],
                xlabel="$\it{K}$", ylabel=key, 
                note=key,
                PATH=pic_path) # maybe add special funcs for certain values
