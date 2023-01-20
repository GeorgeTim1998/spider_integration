import matplotlib.pyplot as mp
import fenics_support as fsup

a = [2.5e-01,
5.0e-01,
1.0e+00]

b = [1.869797004195355,
1.869802181325442,
1.8698486942603931]

mp.plot(a, b)
mp.xlabel("mesh size, sm")
mp.ylabel("S1")
mp.show()

fsup.save_contour_plot('Pics/support/S1_vs_mesh.png')