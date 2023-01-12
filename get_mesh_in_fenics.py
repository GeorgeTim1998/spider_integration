import gmsh
from mpi4py import MPI
from dolfinx.io import XDMFFile, gmshio
from dolfinx import plot

gmsh.initialize()

model = gmsh.model()
sphere = model.occ.addEllipse(0, 0, 0, 2, 1, tag=5)

dim=2
face=model.occ.addCurveLoop([sphere]) # say that curve is closed
surface = model.occ.addPlaneSurface([face]) # say that curve delimits surface
model.occ.synchronize()

model.addPhysicalGroup(dim, [surface]) # add physical group to surface (physical group needed for fenics)
model.addPhysicalGroup(dim-1, [sphere]) # add physical group to curve (physical group needed for fenics)

model.mesh.generate(dim)
gmsh.fltk.run()
#%% 
msh, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_SELF, 0)
tdim = msh.topology.dim # mesh dimensions. For 2d problem equals 2
fdim = tdim - 1 # facets dimensions. for 2d problem equals 2
import pyvista
topology, cell_types, geometry = plot.create_vtk_mesh(msh, tdim)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

plotter = pyvista.Plotter()
plotter.add_mesh(grid, show_edges=True)
plotter.view_xy()
plotter.show()

gmsh.finalize()