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

diff = list(set(msh_files) - set(xml_files))

if diff != []:
  for i, file in enumerate(diff):
    print("\n")
    print(i + 1, 'out of', len(diff))
    bash_command = "dolfin-convert %s/%s.msh %s/%s.xml" % (msh_folder, file, xml_folder, file)
    print(bash_command)
    os.system(bash_command)
  
  print('Converted:')
  for file_name in diff:
    if file_name == diff[-1]:
      print(file_name)
    else:
      print(file_name + ',')
else:
  print('Nothing to convert')