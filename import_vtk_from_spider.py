from fenics import *
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

plot(u)
pyplot.show()