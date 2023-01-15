from fenics import *
import matplotlib.pyplot as pyplot
import support as sup
import fenics_support as fsup
from math import pi

folder = sup.xml_files_folder()
xml_file = folder + '/' + 't1.xml'
xml_file_facet = folder + '/' + 't1_facet_region.xml'
xml_file_physical_region = folder + '/' + 't1_physical_region.xml'

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

plot(u)
pyplot.show()

error_L2 = errornorm(u_D, u, 'L2')

# Compute maximum error at vertices
vertex_values_u_D = u_D.compute_vertex_values(gmsh)
vertex_values_u = u.compute_vertex_values(gmsh)
import numpy as np
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

print(error_max)