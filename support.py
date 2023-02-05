import os
import math
import numpy
import regex
from termcolor import colored, cprint

M0 = 1.25e-6

def msh_files_folder():
  return 'Mesh/msh'

def xml_files_folder():
  return 'Mesh/xml'

def lao_hash():
  Re = 142.5/100 # sm
  # E =  [1.60, 1.59, 1.59, 1.60] # elongation
  E =  numpy.linspace(1, 2, 11)
  betta_p = 0.9 # poloidal betta
  
  a = 37/100
  I = 345e3 # Amperes
  
  Bt0_mean = 2
  p0 = 6400
  Bp = math.sqrt(2*M0*p0 / betta_p)
  
  psi0 = Re**2 * Bp
  F2_0 = Re**2 * Bt0_mean**2
  psi_level = 1e-3 # used to calc contours with known desired psi_level
  Fpl_vs_Fvac_ratio = 1e-3
  
  return {
    'Re': Re,
    'E': E,
    'a': a,
    'I': I,
    'Bt0_mean': Bt0_mean,
    'p0': p0,
    'psi0': psi0,
    'F2_0': F2_0,
    'psi_level': psi_level,
    'Fpl_vs_Fvac_ratio': Fpl_vs_Fvac_ratio 
  }

def create_readmi(lao_hash, folder):
  file = open("%s/readme.txt" % folder, "w")
  keys = list(lao_hash.keys())
  for key in keys:
    file.write("%s = %s\n" % (key, lao_hash[key]))
    
  file.close()
  
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

def mesh_size(filename):
  e_number_pattern = r"\d+\.\d+e\+\d+|\d+\.\d+e\-\d+"
  numbers = regex.findall(e_number_pattern, filename)
  
  return numbers[-1]

def print_colored(color_srt, color='white', str=''):
  print(colored(color_srt, color), str)