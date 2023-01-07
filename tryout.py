from mpi4py import MPI
from dolfinx import mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 8, 8, mesh.CellType.triangle) # 8*8 squares with triangle mesh

from dolfinx.fem import FunctionSpace
V = FunctionSpace(domain, ("CG", 1))

from dolfinx import fem
uD = fem.Function(V)
uD.interpolate(lambda x: 1 + x[0]**2 + 2 * x[1]**2)

import numpy
# Create facet to cell connectivity required to determine boundary facets
tdim = domain.topology.dim # mesh dimensions. For 2d problem equals 2
fdim = tdim - 1 # facets dimensions. for 2d problem equals 2
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology) # find line segments of the boundary of the domain

boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets) # indexes of line segments
bc = fem.dirichletbc(uD, boundary_dofs)

import ufl
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

from petsc4py.PETSc import ScalarType
f = fem.Constant(domain, ScalarType(-6))

a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx

problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

from dolfinx import plot
import pyvista
# topology, cell_types, geometry = plot.create_vtk_mesh(domain, tdim)
# grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

u_topology, u_cell_types, u_geometry = plot.create_vtk_mesh(V)
u_grid = pyvista.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
u_grid.point_data["u"] = uh.x.array.real
u_grid.set_active_scalars("u")
u_plotter = pyvista.Plotter()
u_plotter.add_mesh(u_grid, show_edges=True)

u_plotter.view_xy() # view curtain surface
u_plotter.show()
