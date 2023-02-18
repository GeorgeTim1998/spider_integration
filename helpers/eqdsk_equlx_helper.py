import numpy as np
import math as mt

def write_np_array_to_file(file, array, number_format):
  columns = 5
  rows = mt.floor(len(array)/5)
  separator = columns*rows

  arr5 = array[0: separator]
  arrleftover = array[separator: :]

  arr5 = np.resize(arr5, (rows, columns))
  arrleftover = np.resize(arrleftover, (1, len(arrleftover)))
  np.savetxt(file, arr5, fmt=number_format, delimiter='')
  np.savetxt(file, arrleftover, fmt=number_format, delimiter='')