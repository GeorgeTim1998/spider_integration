import os
import time
from matplotlib import pyplot
import math
import numpy
import regex
from termcolor import colored, cprint
from math import pi, sqrt
import datetime

PICS_FOLDER = 'Pics'
DPI = 240

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
  Fpl_vs_Fvac_ratio = 1e-3
  psi_level = 1e-3
  
  return {
    'Re': Re,
    'E': E,
    'a': a,
    'I': I,
    'Bt0_mean': Bt0_mean,
    'p0': p0,
    'psi0': psi0,
    'psi_level': psi_level,
    'F2_0': F2_0,
    'Fpl_vs_Fvac_ratio': Fpl_vs_Fvac_ratio 
  }

def create_readmi(lao_hash, folder, p_pow=1, F_pow=1):
  file = open("%s/readme.txt" % folder, "w")
  keys = list(lao_hash.keys())
  for key in keys:
    file.write("%s = %s\n" % (key, lao_hash[key]))
  
  file.write("p_pow = %s\n" % p_pow)
  file.write("F_pow = %s" % F_pow)
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
  
def first_line(file):
  psi_size, spacial_size, size, psi_max = next(file).split()
  
  return int(psi_size), int(spacial_size), int(size), float(psi_max)

def second_line(file):
  psi_min, end_entry = next(file).split()
  
  return float(psi_min), int(end_entry)

def return_and_delete_range(data, data_range):
  return_array = data[0:data_range]
  data = numpy.delete(data, numpy.s_[0:data_range])
  
  return return_array, data

def restore_funcpsi(psi, dfuncdpsi, funcb=0):
  psi = numpy.flip(psi)
  dfuncdpsi = numpy.flip(dfuncdpsi)

  funcpsi = numpy.zeros(len(psi))
  
  funcpsi[0] = funcb
  for i in range(1, len(psi)):
    funcpsi[i] = dfuncdpsi[i] * (psi[i]-psi[i-1]) + funcpsi[i-1]
    
  return numpy.flip(funcpsi)

def restore_dpdpsi(dpdpsi):
  pscale = 2 * 2*pi * 1e-7 # it is M0
  # maybe in eqdsk_equlx pressure derivative is not by full flux but by flux/2pi
  return dpdpsi / pscale

def restore_dfdpsi(dfdpsi):
  # maybe in eqdsk_equlx pressure derivative is not by full flux but by flux/2pi
  fscale = 1
  return dfdpsi / fscale

def restore_pres_n_fpol(um, up, meshr, pprime, ffprim, fvac):
  um = -um
  up = -up
  dpsi = (um - up)/(meshr-1)

  pres = numpy.zeros(len(pprime))
  fpol = numpy.zeros(len(pprime))
  fpol[-1] = fvac

  my_steps = numpy.flip(numpy.array(range(0, meshr-1)))
  for i in my_steps:
    pprim_c = -0.5 * (pprime[i+1] + pprime[i])
    fprim_c = -0.5 * (ffprim[i+1] + ffprim[i])
    
    pres[i]=(pres[i+1] + pprim_c * dpsi) # for some reason pres needs 2pi multiplicator to be actual pressure
    fpol[i]=math.sqrt(fpol[i+1]**2 + 2*fprim_c * dpsi) # for some reason pres needs 2pi multiplicator to be actual pressure
    # fpol[i]=sqrt(fpol[i+1]**2 + 2*fprim_c * 2*pi*dpsi)

  return pres, fpol

def Time_name():
  ttime = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
  time_title = str(ttime)  # get current time to make figure name unique
  return time_title

def save_contour_plot(note, PATH=''):
  time_title = Time_name()

  if PATH == '':
    path_my_file = '%s/%s' % (PICS_FOLDER, time_title)
  else:
    path_my_file = '%s/%s' % (PATH, time_title)
  file_path = "%s.png" % path_my_file
  
  pyplot.savefig(file_path, dpi=DPI, bbox_inches="tight")
  pyplot.close()

  print_colored(note, color='green', white_str=file_path)
  time.sleep(1)
  
def countour_plot_maxtrix(r_mesh, z_mesh, psi_mesh, 
                          levels,
                          grid = False,
                          colorbar = False,
                          note='3D countour plot saved to PATH:',
                          PATH='',
                          plot_title='',
                          xlim=[],
                          ylim=[]):

  figure = pyplot.contour(r_mesh, z_mesh, psi_mesh, levels)
  pyplot.gca().set_aspect("equal")
  
  if colorbar:
    pyplot.colorbar(figure).set_label("\u03C8(r, z), Вб")
  
  if xlim != []:
    pyplot.xlim(*xlim)
      
  if ylim != []:
    pyplot.ylim(*ylim)
  
  if grid:
    pyplot.grid(True)
    
  if plot_title != '':
    pyplot.title(plot_title)
  
  pyplot.xlabel("r, cм")
  pyplot.ylabel("z, cм")
  save_contour_plot(note, PATH=PATH)

def print_colored(color_srt, color='white', white_str='', attrs=[]):
  print(colored(color_srt, color, attrs=attrs), white_str)