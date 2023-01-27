from fenics import *
import matplotlib.pyplot as pyplot
import support as sup
import fenics_support as fsup
from math import pi
import numpy as np

folder = sup.xml_files_folder()
lao_hash = sup.lao_hash()
Re = lao_hash['Rt'][0] # ellipse center which sometimes can be R0 in your writings
E = lao_hash['E'] # ellipse elongation
I = lao_hash['I']
ell_a = lao_hash['a']

bp_problem = 0.9 # poloidal betta from lao1985
q_problem = 1 # stability from lao1985

problem_data = fsup.form_dict()

filenames = ['a_37.000_ratio_1.100_msh_1.0e+00', 
             'a_37.000_ratio_1.200_msh_1.0e+00', 
             'a_37.000_ratio_1.300_msh_1.0e+00', 
             'a_37.000_ratio_1.400_msh_1.0e+00', 
             'a_37.000_ratio_1.500_msh_1.0e+00', 
             'a_37.000_ratio_1.600_msh_1.0e+00', 
             'a_37.000_ratio_1.700_msh_1.0e+00', 
             'a_37.000_ratio_1.800_msh_1.0e+00', 
             'a_37.000_ratio_1.900_msh_1.0e+00', 
             'a_37.000_ratio_2.000_msh_1.0e+00']
# filenames = ['a_37.000_ratio_1.100_msh_1.0e+00']
for i, filename in enumerate(filenames):
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
    [r_2, r, z] = fsup.operator_weights(V)
    
    p_part = fsup.linear_pressure(bp_problem, Re)
    F_part = fsup.linear_tor_function(q_problem, E[i], Re)
    
    # [bp, f0_2, psi0] = fsup.linear_profiles()

    a = dot(grad(u)/r, grad(r_2*v))*dx
    # f = r_2 * Constant(pi * bp / Re**2)
    f = Constant(p_part) * r_2 + Constant(F_part)
    L = f * r*v*dx
    
#%% Compute solution and p(psi), F(psi)
    u = Function(V)
    solve(a == L, u, bc)

    inverced_r_integral = fsup.inverced_r_integral(Re, ell_a, ell_b, r, dx, gmsh)
    [u, psi0] = fsup.measure_u(Re, ell_a, ell_b, I, bp_problem, inverced_r_integral, E[i], q_problem, u, V) # de de-measure solution
    
    fsup.countour_plot_via_mesh(gmsh, u, levels = 30, colorbar=True, grid=True)

    [p_sol, p0] = fsup.calculate_p_psi(bp_problem, psi0, u, Re)
    [F_2_sol, F_20] = fsup.calculate_Fpow2_psi(E[i], psi0, q_problem, u, Re)

#%% Post solve calculus
    boundaries = MeshFunction('size_t', gmsh, gmsh.topology().dim() - 1) # get boundaries (and all marked lines???) from mesh
    ds = Measure('ds', domain=gmsh, subdomain_data=boundaries)
    n = FacetNormal(gmsh) # normal to plasma boundary

    L = fsup.boundary_length(ds)
    W = fsup.form_vector_space(u)
    Bp = fsup.calculate_Bp(u, r, W)
    Bpa = 1/L * fsup.circulation(Bp, n, ds)

    [er, ez] = fsup.calculate_orts(W)

    omega = fsup.calculate_omega(r, gmsh)
    S_ = fsup.calculate_plasma_cross_surface(gmsh)
    Spl = fsup.calculate_plasma_surface(r, ds)
    Rt = 1/(2*pi) * omega / S_
    alpha = fsup.calculate_alpha(Bp, ez, gmsh, dx)
    alpha_LB = 2 * E[i]**2 / (E[i]**2 + 1)
    eps_K = (E[i]**2 - 1) / (E[i]**2 + 1)
    R0 = fsup.return_R0(u, V)

    print('Omega =', omega)
    print('Spl =', Spl)
    print('S_ =', S_)
    print('L =', L)
    print('alpha = ', alpha)
    print('alpha_LB =', alpha_LB)

    q = fsup.return_q(r, z, R0, V, W)
    S1, S2, S3 = fsup.calculate_S_integrals(Bpa, omega, Bp, q, n, r, ds)

    print(mesh_size, S1)
    print(mesh_size, S2)
    print(mesh_size, S3)
    
#%% Calc magnetic values
    bp, li, mu_i = fsup.solve_SLAE(alpha, [S1, S2, S3], Rt, R0)
    print([bp, li, mu_i])

    vars = globals()
    keys = list(problem_data.keys())
    for key in keys:
        problem_data[key].append(vars[key])

#%% post problem plot
keys = list(problem_data.keys())
for key in keys:
    fsup.plot_1D(E, problem_data[key], xlabel='elongation', ylabel=key) # maybe add special funcs for certain values
