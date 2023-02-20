from load_data_from_matlab_scripts import LoadDataFromMatlabScripts
import matplotlib.pyplot as matplt
import numpy as np
import math as mh
import datetime
import time
from fenics import *
from scipy.interpolate import interp2d
from fenics_support import countour_plot_via_mesh

class ExpressionFromScipyFunction(UserExpression):
    def __init__(self, f, **kwargs):
        self._f = f
        UserExpression.__init__(self, **kwargs)
    def eval(self, values, x):
        values[:] = self._f(*x)

matlab_data = LoadDataFromMatlabScripts()

#%% create mesh
min_point = Point(matlab_data.r.min(), matlab_data.z.min())
max_point = Point(matlab_data.r.max(), matlab_data.z.max())
gmsh = RectangleMesh(min_point, max_point, len(matlab_data.z) - 1, len(matlab_data.r) - 1)

#%% create psi
V   = FunctionSpace(gmsh, 'Lagrange', 1)

r = matlab_data.r
z = matlab_data.z

interpolant = interp2d(r, z, matlab_data.psi, kind='linear', copy=False, bounds_error=True)

expression = ExpressionFromScipyFunction(interpolant, element=V.ufl_element())
expression = interpolate(expression, V) 

countour_plot_via_mesh(gmsh, expression, levels=np.linspace(0.15, 0.68, 10), colorbar=True, xlim=[0.8, 2.4], ylim=[-1, 1], grid=True)
