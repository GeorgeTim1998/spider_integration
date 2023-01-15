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
  