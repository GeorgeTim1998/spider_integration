import os 
import support as sup
import math

msh_folder = sup.msh_files_folder()
msh_files = sup.drop_files_ext(os.listdir(msh_folder))

xml_folder = sup.xml_files_folder()
xml_files = sup.drop_files_ext(os.listdir(xml_folder))
xml_files = sup.drop_physical_n_dacet_files(xml_files)

msh_files.sort()
xml_files.sort()

if msh_files != xml_files:
  convert_indexes = []
  
  if len(xml_files) == 0:
    convert_indexes = range(len(msh_files))
  
  for i, xml in enumerate(xml_files):
    for j, msh in enumerate(msh_files):
      if msh != xml:
        convert_indexes.append(j)     
      
  print(len(convert_indexes))
  
  for i, index in enumerate(convert_indexes):
    print("\n")
    print(math.floor(i/len(convert_indexes)*100))
    bash_command = "dolfin-convert %s/%s.msh %s/%s.xml" % (msh_folder, msh_files[index], xml_folder, msh_files[index])
    os.system(bash_command)