import numpy

def solve_SLAE(alpha):
  matrix = numpy.array(
    [
      [3,            1, -1],
      [1,            1,  1],
      [1, -(alpha - 1), -1]
    ]
  )
  
  f = numpy.array([3.28681356, 1.25700438, 0.99833036])
  
  x = numpy.linalg.solve(matrix, f)
  
  return x

# alpha = 1.0226596521091218
alpha = 1.0216596521091218
print(solve_SLAE(alpha), 'alpha = ', alpha)