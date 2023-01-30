"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.

  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary

  u_D = 1 + x^2 + 2y^2
    f = -6
"""

from fenics import *
import matplotlib.pyplot as plt
import numpy as np
# from termcolor import colored

# Create mesh and define function space
gmsh = UnitSquareMesh(8, 8)
V = FunctionSpace(gmsh, 'P', 1)

# Define boundary condition
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

max_index = np.argmax(u.vector()[:]) # return index of a max item
max_point = V.tabulate_dof_coordinates() # can be used for assigning?? from known data

# Plot solution and gmsh
plot(u)
# plot(gmsh)

# Compute error in L2 norm
plt.show()

# короче схема рабочая. тебе в set_local надо по нормальному просто присвоить. надо это делать в последовательности, dof коэффициентов