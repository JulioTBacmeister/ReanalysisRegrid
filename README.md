# ReanalysisRegrid
Minimal (?) set of functions to regrid reanalyses (currently only ERA5 ... working on ERAi, MERRA ...)

How to use (Here we assume './' is the 'root' directory of your clone)

1) Go to ./Drivers

2) cp config_ERA5regrid_template.yaml config_ERA5regrid.yaml

3) edit config_ERA5regrid.yaml file according to instructions given in file

4) qsub PyBatch_ERA5regrid.csh

Note: As the jobs run your config_ERA5regrid.yaml file will updated. All comments in template will be stripped off etc. 


# Background information
This code uses ESMF routines to regrid horizontally and scipy.interpolate functions to interpolate vertically.
In order to accomplish these tasks, grid information is needed.  Horizontal grids are described via SCRIP files. 
Vertical grids are described via "hybrid a and b" coefficients for SE and FV dycores, and via the "zgrid" 
variable for MPAS grids  (found in ncdata or history files). In addition, for correct regridding in the vertical,
topography is also needed.

This information is currently obtained using grid 'nicknames' such as 'ne30pg3' for the horizontal grids, and 'L93' 
for vertical grids.  These nicknames are supplied to a function 'gridInfo' in the module ./Utils/GridUtils.py, which 
returns the correct SCRIP, vertical grid, and topography information to the regridding functions.