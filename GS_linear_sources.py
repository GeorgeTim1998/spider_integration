from fenics import *
import matplotlib.pyplot as pyplot
import support as sup
import fenics_support as fsup
from math import pi

folder = sup.xml_files_folder()
filename = 'a_37.000_ratio_1.600_msh_1.0e+00'
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

[r_2, r] = fsup.operator_weights(V)
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

lao_hash = sup.lao_hash()
z = interpolate(Expression('x[1]', degree = 1), V) # interpolation is needed so that 'a' could evaluate deriviations and such
q = []
er = as_vector((interpolate(Constant(1), V), interpolate(Constant(0), V)))
# pick one
# q.append( (r - interpolate(Constant(lao_hash['Rt'][0]), V))*er + z*ez )
q.append( as_vector((r - interpolate(Constant(lao_hash['Rt'][0]), V), z))) 
# q.append( as_vector((r - interpolate(Constant(lao_hash['Rt'][0]), V), 0))) 

S1 = 1 / Bpa**2 / omega * assemble( dot(Bp, Bp) * dot(q[0], n) * 2*pi*r*ds )

print(mesh_size, S1)