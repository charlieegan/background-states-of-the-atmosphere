import numpy as np
import atmosphere_bgs
import xarray as xr

# Compare lifecycle background states: zonal mean, OT, ELIPVI
lc1_OT = xr.open_dataset('./output/lc1low_2001010000_OT.nc')
lc1_zm = xr.open_dataset('./output/lc1_2001010000_zonal_mean.nc')
lc1_ELIPVI = xr.open_dataset('./output/lc1_2001010000_ELIPVI.nc')