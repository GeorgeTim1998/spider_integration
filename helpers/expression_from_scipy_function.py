from fenics import UserExpression

class ExpressionFromScipyFunction(UserExpression):
  def __init__(self, f, **kwargs):
      self._f = f
      UserExpression.__init__(self, **kwargs)
  def eval(self, values, x):
      values[:] = self._f(*x)