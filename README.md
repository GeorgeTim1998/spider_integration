## Project name: SPIDER integration
- Integrates solutions obtained by fortran code SPIDER and fenics python package for calcucating viral relations along plasma boundary
- Allow to create eqdsk_equilx files used as input files for SPIDER solver
- Has interface that allows to launch SPIDER solver from python console 
- Creates meshes with gmsh mesh generator to import SPIDER solutions to
- Calculates viral relations over plasma boundary for SPIDER code solutions imported to gmsh mesh with fenics integrals package

## How to set up project
- Clone the repository and navigate to created folder
- `conda env create -f environment.yml`

## How to launch
- `conda activate fenicsproject`
