import gmsh
import os
import support as sup
import numpy as np

FOLDER = sup.msh_files_folder()
DIM = 2
MESH_SIZE = 5e-1
LC = 0.5e-1

def create_gmsh_mesh_from_points(boundary, ell_a=0, elongation=0, msh_size=MESH_SIZE, show=False):
  boundary = np.reshape(boundary, (-1, 2))

  file_name = "a_%.3f_ratio_%.3f_msh_%.1e" % (ell_a, elongation, msh_size)
  
  if not mesh_exists(file_name):
    gmsh.initialize()
    gmsh.option.setNumber('Mesh.MeshSizeMax', msh_size) # you can set general options via strings in set number. see https://gmsh.info/doc/texinfo/gmsh.pdf
    # gmsh.option.setNumber('Mesh.Algorithm', 3) # you can set general options via strings in set number. see https://gmsh.info/doc/texinfo/gmsh.pdf
    model = gmsh.model()
    
    point_tags = []
    for point in boundary:
      point_tags.append(model.occ.addPoint(point[0], point[1], 0, LC)) # it is looped in SPIDER
    
    line_tags = []
    for i in range(len(point_tags) - 1):
      line_tags.append(model.occ.addLine(point_tags[i], point_tags[i+1]))
    
    curve = model.occ.addCurveLoop(line_tags)
    surface = model.occ.addPlaneSurface([curve])
    model.occ.synchronize()  
    
    model.addPhysicalGroup(DIM, [surface]) # add physical group to surface (physical groups needed for fenics)
    model.addPhysicalGroup(DIM-1, [curve]) # add physical group to curve (physical groups needed for fenics)
    
    model.mesh.generate(DIM)
    if show:
      gmsh.fltk.run()
    
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2) # setting ascii 2 version so fenics understands whats up
    gmsh.write("%s/%s.msh" % (FOLDER, file_name))
    gmsh.finalize()
  
def mesh_exists(filename):
  return (filename + '.msh') in os.listdir(FOLDER)