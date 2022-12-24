import numpy as np

class LoadDataFromMatlabScripts:
  ARRAYS_FILES = ['matlab_scripts_output/20221223/r.txt',
                   'matlab_scripts_output/20221223/z.txt',
                   'matlab_scripts_output/20221223/rxb.txt',
                   'matlab_scripts_output/20221223/zxb.txt']
  MATRIXES_FILES = ['matlab_scripts_output/20221223/psi.txt']
  
  def __init__(self):
    self.r, self.z, self.rxb, self.zxb = self.load_arrays_from_m()
    self.psi = self.load_matrixes_from_m()[0]
    
  def load_arrays_from_m(self):
    data = []
    for file_name in LoadDataFromMatlabScripts.ARRAYS_FILES:
      file = open(file_name, "r")
      data.append(np.fromfile(file, sep=','))
      file.close()
    
    return data
  
  def load_matrixes_from_m(self):
    data = []
    for file_name in LoadDataFromMatlabScripts.MATRIXES_FILES:
      file = open(file_name, "r")
      data.append(np.loadtxt(file, delimiter=','))
      file.close()
    
    return data
    
