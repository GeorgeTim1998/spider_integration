import gmsh
from dolfinx.io import XDMFFile, gmshio
from mpi4py import MPI
import sys

# Initialize gmsh:
gmsh.initialize()
 
# points:
lc = 5e-2
model = gmsh.model()
point1 = model.geo.add_point(0, 0, 0)
point2 = model.geo.add_point(1, 0, 0)
point3 = model.geo.add_point(1, 1, 0)
point4 = model.geo.add_point(0, 1, 0)

# lines
line1 = model.geo.add_line(point1, point2)
line2 = model.geo.add_line(point2, point3)
line3 = model.geo.add_line(point3, point4)
line4 = model.geo.add_line(point4, point1) 

# faces 
face1 = model.geo.add_curve_loop([line1, line2, line3, line4])

# surfaces of cube:
model.geo.add_plane_surface([face1])

# Create the relevant Gmsh data structures
# from Gmsh model.
model.geo.synchronize()

# Generate mesh:
model.mesh.generate(2)

# Write mesh data:
# gmsh.write("Mesh/FirstMesh.msh")

# # Creates graphical user interface
# if 'close' not in sys.argv:
# gmsh.fltk.run()
 
msh, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_SELF, 0)
# msh.name = "Sphere"
# cell_markers.name = f"{msh.name}_cells"
# facet_markers.name = f"{msh.name}_facets"
#%% 
# msh, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_SELF, 0)
# It finalize the Gmsh API
gmsh.finalize()