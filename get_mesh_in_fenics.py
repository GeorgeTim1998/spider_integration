import gmsh
from mpi4py import MPI
from dolfinx.io import XDMFFile, gmshio
from dolfinx import plot

# Initialize gmsh:
# gmsh.initialize()
# model = gmsh.model()

# # points:
# lc = 5e-2
# point1 = model.geo.add_point(0, 0, 0)
# point2 = model.geo.add_point(1, 0, 0)
# point3 = model.geo.add_point(1, 1, 0)
# point4 = model.geo.add_point(0, 1, 0)

# # lines
# line1 = model.geo.add_line(point1, point2)
# line2 = model.geo.add_line(point2, point3)
# line3 = model.geo.add_line(point3, point4)
# line4 = model.geo.add_line(point4, point1) 
# # faces 
# face1 = model.geo.add_curve_loop([line1, line2, line3, line4])
# # surfaces of cube:
# model.geo.add_plane_surface([face1])
# # Create the relevant Gmsh data structures
# # from Gmsh model.
# model.geo.synchronize()
# # Generate mesh:
# model.mesh.generate(2)
# # It finalize the Gmsh API
# gmsh.fltk.run()
# gmsh.finalize()
gmsh.initialize()

# Choose if Gmsh output is verbose
gmsh.option.setNumber("General.Terminal", 0)
model = gmsh.model()
model.add("Sphere")
model.setCurrent("Sphere")
sphere = model.occ.addSphere(0, 0, 0, 1, tag=1)

# Synchronize OpenCascade representation with gmsh model
model.occ.synchronize()

# Add physical marker for cells. It is important to call this function
# after OpenCascade synchronization
model.add_physical_group(3, [sphere])


# Generate the mesh
model.mesh.generate(3)

msh, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_SELF, 0)
msh.name = "Sphere"
cell_markers.name = f"{msh.name}_cells"
facet_markers.name = f"{msh.name}_facets"
#%% 
msh, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, MPI.COMM_SELF, 0)
tdim = msh.topology.dim # mesh dimensions. For 2d problem equals 2
fdim = tdim - 1 # facets dimensions. for 2d problem equals 2
import pyvista
topology, cell_types, geometry = plot.create_vtk_mesh(msh, tdim)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

plotter = pyvista.Plotter()
plotter.add_mesh(grid, show_edges=True)
plotter.view_xy()
plotter.show()
