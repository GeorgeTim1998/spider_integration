from fenics import *
import matplotlib.pyplot as pyplot
import support as sup
import fenics_support as fsup
from math import pi
import numpy as np

folder = sup.xml_files_folder()
lao_hash = sup.lao_hash()
Re = lao_hash['Rt'][0] # ellipse center
E = lao_hash['E'] # ellipse elongation
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

    [r_2, r, z] = fsup.operator_weights(V)
    [bp, f0_2, psi0] = fsup.linear_profiles()

    a = dot(grad(u)/r, grad(r_2*v))*dx
    f = r_2 * Constant(pi * bp / Re**2)
    L = f * r*v*dx
#%% Compute solution
    u = Function(V)
    solve(a == L, u, bc)

    fsup.countour_plot_via_mesh(gmsh, u, levels = 30, colorbar=True, grid=True)

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
    alpha_LB = 2 * E[i]**2 / (E[i]**2 + 1)
    R0 = fsup.return_R0(u, V)

    print('Omega =', omega)
    print('Spl =', Spl)
    print('S_ =', S_)
    print('L =', L)
    print('[1,0] circ = ', fsup.circulation(er, n, ds))
    print('[z,0] circ = ', fsup.circulation(z*er, n, ds))

    q = fsup.return_q(r, z, R0, V, W)
    S1, S2, S3 = fsup.calculate_S_integrals(Bpa, omega, Bp, q, n, r, ds)

    print(mesh_size, S1)
    print(mesh_size, S2)
    print(mesh_size, S3)
    
#%% Calc magnetic values
    bp, li, mu_i = fsup.solve_SLAE(alpha_LB, [S1, S2, S3], Rt, R0)
    print([bp, li, mu_i])
    
    vars = globals()
    keys = list(problem_data.keys())
    for key in keys:
        problem_data[key].append(vars[key])

fsup.plot_1D(E, problem_data['omega'], xlabel='elongation', ylabel='Plasma volume, m**3')
fsup.plot_1D(E, problem_data['S_'], xlabel='elongation', ylabel='Plasma cross surface, m**2')
fsup.plot_1D(E, problem_data['Spl'], xlabel='elongation', ylabel='Plasma surface, m**2')

fsup.plot_1D(E, problem_data['Rt'], xlabel='elongation', ylabel='Rt')
fsup.plot_1D(E, problem_data['alpha_LB'], xlabel='elongation', ylabel='2E**2 / (E**2 + 1)')
fsup.plot_1D(E, problem_data['R0'], xlabel='elongation', ylabel='Magnetic axis, m')

fsup.plot_1D(E, problem_data['S1'], xlabel='elongation', ylabel='S1')
fsup.plot_1D(E, problem_data['S2'], xlabel='elongation', ylabel='S2')
fsup.plot_1D(E, problem_data['S3'], xlabel='elongation', ylabel='S3')

fsup.plot_1D(E, problem_data['bp'], xlabel='elongation', ylabel='bp')
fsup.plot_1D(E, problem_data['li'], xlabel='elongation', ylabel='li')
fsup.plot_1D(E, problem_data['mu_i'], xlabel='elongation', ylabel='mu_i')
