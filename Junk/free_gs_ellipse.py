import numpy
import freegs
import matplotlib.pyplot as plt

tokamak = freegs.machine.TestTokamak()

eq = freegs.Equilibrium(tokamak=tokamak,
                     Rmin=0.1, Rmax=2.0,    # Radial domain
                     Zmin=-1.0, Zmax=1.0,   # Height range
                     nx=65, ny=65,          # Number of grid points
                     boundary=freegs.boundary.fixedBoundary)

profiles = freegs.jtor.ConstrainBetapIp(0.5, # Poloidal beta
                                        1e6, # Plasma current [Amps]
                                        1.0) # Vacuum f=R*Bt

freegs.solve(eq,          # The equilibrium to adjust
           profiles,    # The toroidal current profile function
           psi_bndry=0.0)  # Because no X-points, specify the separatrix psi

eq.plot(show=True)