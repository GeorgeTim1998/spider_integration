import gmsh

folder = 'Mesh/msh'

dim=2
mesh_size = 0.01
ellipse_center = [0, 0, 0]
ellipse_ratio = 3
ellipse_a = 3
ellipse_b = 1

gmsh.initialize()
gmsh.option.setNumber('Mesh.MeshSizeMax', mesh_size) # you can set general options via strings in set number. see https://gmsh.info/doc/texinfo/gmsh.pdf

model = gmsh.model()
ellipse = model.occ.addEllipse(*ellipse_center, ellipse_a, ellipse_b, zAxis=[0, 0, 1], xAxis=[0, 1, 0])
loop = model.occ.addCurveLoop([ellipse]) # say that curve is closed/looped
surface = model.occ.addPlaneSurface([loop]) # say that curve delimits surface
model.occ.synchronize()

model.addPhysicalGroup(dim, [surface]) # add physical group to surface (physical groups needed for fenics)
model.addPhysicalGroup(dim-1, [ellipse]) # add physical group to curve (physical groups needed for fenics)

model.mesh.generate(dim)

gmsh.option.setNumber('Mesh.MshFileVersion', 2.2) # setting ascii 2 version so fenics understands whats up
gmsh.write(folder + '/' + "filename.msh")

gmsh.fltk.run()
gmsh.finalize()
