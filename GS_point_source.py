from fenics import *
import matplotlib.pyplot as pyplot
import support as sup
import fenics_support as fsup
from math import pi
import numpy as np

folder = sup.xml_files_folder()

# filenames = ['a_37.000_ratio_1.600_msh_1.0e-01', 'a_37.000_ratio_1.600_msh_2.5e-01', 'a_37.000_ratio_1.600_msh_5.0e-01', 'a_37.000_ratio_1.600_msh_1.0e+00']
# filenames = ['a_37.000_ratio_1.600_msh_2.5e-01', 'a_37.000_ratio_1.600_msh_5.0e-01', 'a_37.000_ratio_1.600_msh_1.0e+00']
filenames = ['a_37.000_ratio_1.000_msh_1.0e+00']
# for filename in filenames:
filename = filenames[-1]
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

a = dot(grad(u)/r, grad(r_2*v))*dx
I = 1
sigma = 1
R0pi=142.5
f = fsup.point_source(I=I, sigma=sigma, R0=R0pi)
L = fsup.M0*r * f * r*v*dx
#%% Compute solution
u = Function(V)
solve(a == L, u, bc)

# fsup.countour_plot_via_mesh(gmsh, u, levels = 30, colorbar=True, grid=True)

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
R0 = fsup.return_R0(u, V)

a = 10
print("Bp(top) = %e" % Bp(142.5, a)[0])
print("Bp(top) = %e" % ( fsup.M0*I / (2*pi*a) ))

# print('Omega =', omega)
# print('Spl =', Spl)
# print('S_ =', S_)
# print('L =', L)
# print("Bp*dl = %e" % fsup.circulation(Bp, n, ds))
# print("Bp*dl = %e" % (fsup.M0 * I))
# print("Bp*dl = %e" % ( fsup.circulation(Bp, n, ds) / (fsup.M0 * I) ))
# print("Bp*dl = %e" % ( fsup.circulation(Bp, n, ds) / (fsup.M0 * I) )**-1)

q = fsup.return_q(r, z, R0, V, W)
S = fsup.calculate_S_integrals(Bpa, omega, Bp, q, n, r, ds)

# print(mesh_size, S['S1'])
# print(mesh_size, S['S2'])
# print(mesh_size, S['S3'])
