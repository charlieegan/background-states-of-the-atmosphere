# Function to reconstruct the optimal Laguerre diagram fromt the saved netCDF file
import numpy as np
import _atmosphere_bgs
import xarray as xr

def reconstruct_laguerre_diagram(filepath):
    '''Function to reconstruct the optimal Laguerre diagram from the saved netCDF file.'''
    # Open the netCDF file
    ds = xr.open_dataset(filepath)
    
    # Get default physical and simulation parameters
    pp = _atmosphere_bgs.PhysicalParameters()
    sp = _atmosphere_bgs.SimulationParameters()
    
    # Overwrite domain boundary parameters
    sp.smin = ds['smin'].values
    sp.smax = ds['smax'].values
    pp.p00 = ds['p00'].values

    sp.pmin = ds['p00'].values
    
    # Create a new instance of the LaguerreDiagram class
    ld = _atmosphere_bgs.LaguerreDiagram(ds['ld_seeds'].values, ds['ld_duals'].values,pp, sp)
    
    return ld

data_types = ['lc1low','lc2low','ERA5']
dates = ['200101','200101','200911']
dtimes =['0000','0000','0100'] 

for data_type, date, dtime in zip(data_types,dates,dtimes):

    filepath = '/home/users/cq934523/ENMs/background-states-of-the-atmosphere/output/'+data_type+'_'+date+dtime
    ld = reconstruct_laguerre_diagram(filepath+'.nc')
    