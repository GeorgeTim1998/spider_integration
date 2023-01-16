import os
def msh_files_folder():
  return 'Mesh/msh'

def xml_files_folder():
  return 'Mesh/xml'

def lao_hash():
  Rt = [142.5, 142.5, 142.6, 142.3] # sm
  E =  [1.60, 1.59, 1.59, 1.60] # elongation
  betta_p = [0.91, 0.91, 0.90, 0.92] # poloidal betta
  a = 37
  
  return {
    'Rt': Rt,
    'E': E,
    'betta_p': betta_p,
    'a': a
  }

def drop_files_ext(files):
  for i, file in enumerate(files):
    files[i] = os.path.splitext(file)[0]
    
  return files

def drop_physical_n_dacet_files(xml_files):
  physical = 'physical'
  facet = 'facet'
  delete_indexes = []
  
  for i, file in enumerate(xml_files):
    if physical in file or facet in file:
      delete_indexes.append(i)
  
  for i in list(reversed(delete_indexes)):
    xml_files.remove(xml_files[i])
    
  return xml_files