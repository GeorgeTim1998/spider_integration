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