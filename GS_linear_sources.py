from fenics import *
import matplotlib.pyplot as pyplot
import support as sup
import fenics_support as fsup
from math import pi
import numpy as np

folder = sup.xml_files_folder()

# filenames = ['a_37.000_ratio_1.600_msh_1.0e-01', 'a_37.000_ratio_1.600_msh_2.5e-01', 'a_37.000_ratio_1.600_msh_5.0e-01', 'a_37.000_ratio_1.600_msh_1.0e+00']
filenames = ['a_37.000_ratio_1.600_msh_2.5e-01', 'a_37.000_ratio_1.600_msh_5.0e-01', 'a_37.000_ratio_1.600_msh_1.0e+00']
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
[p0, f0_2, psi0] = fsup.linear_profiles()

a = dot(grad(u)/r, grad(r_2*v))*dx
f = Constant( (4 * pi * p0 + 0.5 * f0_2) / psi0 )
L = f*v*dx
#%% Compute solution
u = Function(V)
solve(a == L, u, bc)

# fsup.countour_plot_via_mesh(gmsh, u, levels = 30, colorbar=True, grid=True)

#%% Post solve calculus
boundaries = MeshFunction('size_t', gmsh, gmsh.topology().dim() - 1) # get boundaries (and all marked lines???) from mesh
ds = Measure('ds', domain=gmsh, subdomain_data=boundaries)
n = FacetNormal(gmsh) # normal to plasma boundary

L = fsup.boundary_length(ds)
[Bp, W] = fsup.calculate_Bp(u, r)
Bpa = 1/L * fsup.circulation(Bp, n, ds)

[er, ez] = fsup.calculate_orts(W)

omega = fsup.calculate_omega(r, gmsh)
S_plasma = fsup.calculate_plasma_cross_surface(gmsh)
Rt = 1/(2*pi) * omega / S_plasma

R0 = fsup.return_R0(u, V)

q = fsup.return_q(r, z, R0, V)
S = fsup.calculate_S_integrals(Bpa, omega, Bp, q, n, r, ds)

print(mesh_size, S['S1'])
print(mesh_size, S['S2'])
print(mesh_size, S['S3'])
