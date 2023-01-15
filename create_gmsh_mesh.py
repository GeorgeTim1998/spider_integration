import gmsh
import math

folder = 'Mesh/msh'

gmsh.initialize()

model = gmsh.model()
sphere = model.occ.addEllipse(0, 0, 0, 2, 1, tag=5, zAxis=[0, 0, 1], xAxis=[0, 1, 0])

dim=2

face=model.occ.addCurveLoop([sphere]) # say that curve is closed
surface = model.occ.addPlaneSurface([face]) # say that curve delimits surface

model.occ.synchronize()

model.addPhysicalGroup(dim, [surface]) # add physical group to surface (physical group needed for fenics)
model.addPhysicalGroup(dim-1, [sphere]) # add physical group to curve (physical group needed for fenics)

model.mesh.generate(dim)

gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)   
gmsh.write(folder + '/' + "filename.msh")

gmsh.fltk.run()
gmsh.finalize()
