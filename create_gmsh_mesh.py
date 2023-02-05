import gmsh
import support as sup

folder = sup.msh_files_folder()
lao_hash = sup.lao_hash()

file_names = []
dim = 2
mesh_size = 1e-2
ellipse_center = [lao_hash['Re'][0], 0, 0]
# E_array = lao_hash['E']
E_array = [1.0, 1.1]

for ellipse_ratio in E_array:
  ellipse_a = lao_hash['a']
  ellipse_b = ellipse_ratio * ellipse_a

  gmsh.initialize()
  gmsh.option.setNumber('Mesh.MinimumCircleNodes', 200) # you can set general options via strings in set number. see https://gmsh.info/doc/texinfo/gmsh.pdf

  model = gmsh.model()
  ellipse = model.occ.addEllipse(*ellipse_center, ellipse_b, ellipse_a, zAxis=[0, 0, 1], xAxis=[0, 1, 0])
  loop = model.occ.addCurveLoop([ellipse]) # say that curve is closed/looped
  surface = model.occ.addPlaneSurface([loop]) # say that curve delimits surface
  point = model.occ.addPoint(*ellipse_center, meshSize=mesh_size*3)
  model.occ.synchronize()

  model.mesh.embed(0, [point], 2, surface) # embed point 5 (dim 0) in surface 6 (dim 2)
  
  model.addPhysicalGroup(dim, [surface]) # add physical group to surface (physical groups needed for fenics)
  model.addPhysicalGroup(dim-1, [ellipse]) # add physical group to curve (physical groups needed for fenics)
  
  model.mesh.generate(dim)
  # gmsh.fltk.run()

  file_name = "a_%.3f_ratio_%.3f_msh_%.1e" % (ellipse_a, ellipse_ratio, mesh_size)
  file_names.append(file_name)
  gmsh.option.setNumber('Mesh.MshFileVersion', 2.2) # setting ascii 2 version so fenics understands whats up
  gmsh.write("%s/%s.msh" % (folder, file_name))

  gmsh.finalize()
  print("\n")

for file_name in file_names:
  print(file_name)