from fenics import *

def operator_weights(self, V):
  r_2 = interpolate(Expression('x[0]*x[0]', degree = 2), V) # interpolation is needed so that 'a' could evaluate deriviations and such
  r = Expression('x[0]', degree = 1) # interpolation is needed so that 'a' could evaluate deriviations and such
  
  return r_2, r

def linear_profiles():
  p0 = 1
  f0_2 = 1
  psi0 = 1
  
  return p0, f0_2, psi0

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