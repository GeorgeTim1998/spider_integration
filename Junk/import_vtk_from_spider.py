from fenics import *
import fenics_support as fsup
import matplotlib.pyplot as pyplot

gmsh = UnitSquareMesh(8, 8)
V = FunctionSpace(gmsh, 'P', 1)

u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
def boundary(x, on_boundary):
  return on_boundary
bc = DirichletBC(V, u_D, boundary)
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx
# Compute solution
u = Function(V)
solve(a == L, u, bc)

fsup.countour_plot_via_mesh(gmsh, u, levels=5, colorbar=True)

gmsh2 = RectangleMesh(Point(0.1, 0.1), Point(0.7, 0.7), 10, 10)
VV = FunctionSpace(gmsh2, 'P', 1)
uu = project(u, VV)

fsup.countour_plot_via_mesh(gmsh2, uu, levels=[1, 1.5, 2, 2.5, 3, 3.5, 4], colorbar=True)