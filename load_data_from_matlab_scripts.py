class LoadDataFromMatlabScripts:
  MATLAB_SCRIPTS_OUTPUT = ['matlab_scripts_output/20221223/r.txt',
                           'matlab_scripts_output/20221223/z.txt',
                           'matlab_scripts_output/20221223/rxb.txt',
                           'matlab_scripts_output/20221223/zxb.txt',
                           'matlab_scripts_output/20221223/psi.txt']
  
  def __init__(self):
    text_file1 = open(LoadDataFromMatlabScripts.MATLAB_SCRIPTS_OUTPUT[0], "r")
    text_file2 = open(LoadDataFromMatlabScripts.MATLAB_SCRIPTS_OUTPUT[1], "r")
    text_file3 = open(LoadDataFromMatlabScripts.MATLAB_SCRIPTS_OUTPUT[2], "r")
    text_file4 = open(LoadDataFromMatlabScripts.MATLAB_SCRIPTS_OUTPUT[3], "r")
    
    self.r = text_file1.readlines()
    self.z = text_file2.readlines()
    self.rxb = text_file3.readlines()
    self.zxb = text_file4.readlines()
    self.psi = self.read_matrix(LoadDataFromMatlabScripts.MATLAB_SCRIPTS_OUTPUT[4])
    
    text_file1.close()
    text_file2.close()
    text_file3.close()
    text_file4.close()
    
  def read_matrix(self, file_name):
    with open(file_name, 'r') as f:
      matrix = [[num for num in line.split(',')] for line in f]
      
    
