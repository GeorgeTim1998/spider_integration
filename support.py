import os
import numpy
import regex
from termcolor import colored, cprint

def msh_files_folder():
  return 'Mesh/msh'

def xml_files_folder():
  return 'Mesh/xml'

def lao_hash():
  Re = [142.5/100, 142.5/100, 142.6/100, 142.3/100] # sm
  # E =  [1.60, 1.59, 1.59, 1.60] # elongation
  E =  numpy.linspace(1, 2, 11)
  betta_p = [0.91, 0.91, 0.90, 0.92] # poloidal betta
  a = 37/100
  I = 345e3 # Amperes
  
  return {
    'Re': Re,
    'E': E,
    'betta_p': betta_p,
    'a': a,
    'I': I
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

def mesh_size(filename):
  e_number_pattern = r"\d+\.\d+e\+\d+|\d+\.\d+e\-\d+"
  numbers = regex.findall(e_number_pattern, filename)
  
  return numbers[-1]

def print_colored(color_srt, color='white', str=''):
  print(colored(color_srt, color), str)