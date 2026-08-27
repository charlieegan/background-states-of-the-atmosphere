# This script obtains gridded MLM states and saves them as netCDF files.
# It uses the optimal transport method followed by post-processing steps.

import numpy as np
import atmosphere_bgs
import xarray as xr
import sys
import importlib.util

# 0. Define the simulation parameters
ot_tol = 1e-4       # tolerance for optimal transport problem
eps = None          # smoothing parameter, if None then set depending on data type later on
res = [2**8,2**8]   # resolution of gridded data in latitude and log(pressure)

# 1. Load the input data and define where to save data
print('Loading input data')
data_types = ['lc1low']
dates = ['200101']
dtimes = ['0000']

data_type = 'lc2low'
date = '200101'
dtime = '0000'

data_type = 'ERA5'
date = '200911'
dtime = '0100'

data_types = ['lc1low','lc2low','ERA5']
dates = ['200101','200101','200911']
dtimes =['0000','0000','0100'] 

for data_type, date, dtime in zip(data_types,dates,dtimes):

    save_name = '/home/users/cq934523/ENMs/background-states-of-the-atmosphere/output/'+data_type+'_'+date+dtime

    if 'lc' in data_type:
        filepath = '/storage/research/diamet/swrmethn/inv3/tradv/LC/wadiag_highres/'
        filename = 'bs_'+data_type+date+dtime
        p00 = 1e+5
        pmin = 10
        
    elif data_type == 'ERA5':
        filepath = '/storage/research/s2senm/wd914288/ENM_data/step2_integrals/'+str(date)[:4]+'/nbs/'    
        filename = 'bs_diout'+date+dtime
        p00 = None
        pmin = 1

    input_data = atmosphere_bgs.DataLoader(filepath+filename,p00=p00,load_all=True,pmin=pmin)

    # 2. Solve the optimal transport problem to get Laguerre diagram
    print('Solving optimal transport problem')
    solv = atmosphere_bgs.OTSolver(input_data)
    solv.ot_tol = ot_tol
    solv.get_bgs(verbose=False)

    # 3. Apply post processing steps to obtain gridded smooth solution in geophysical (s,p) and isentropic (lat,theta) coordinates
    print('Obtaining smooth gridded solution')
    bscirc = input_data.bscirc
    pvlevs = input_data.pvlevs
    thlevs = input_data.thlevs
    pvmaxth = input_data.data_dict['MAX PV ON THETA LEVELS']

    if eps is None:
        if 'lc' in data_type:
            eps = 1e-3
        elif data_type == 'ERA5':
            eps = 5e-3
        
    mlm = atmosphere_bgs.PostProcessor(solv,bscirc,pvlevs,thlevs,pvmaxth,eps=eps,nb_type='Lag_cell_8',res=res,make_surface_adjustments=True)

    # 4. Put gridded data into an xarray Dataset
    print('Saving gridded data in netCDF format')
    
    ## Define data dictionary
    data_dict = dict()
    
    ## Coordinates
    data_dict['latitude_levels'] = np.unique(mlm.lg)
    data_dict['pressure_levels'] = np.unique(mlm.pg)
    data_dict['potential_temperature_levels'] = np.sort(mlm.thlevs)
    
    ## Interior variables in latitude/pressure coordinates
    data_dict['zonal_angular_momentum'] = surf_mask*mlm.Z_eps.T
    data_dict['zonal_wind'] = surf_mask*mlm.u_eps.T
    data_dict['potential_temperature'] = surf_mask*mlm.th_eps.T
    data_dict['geopotential'] = surf_mask*mlm.Phi_eps.T
    
    ## Interior variables in latitude/theta coordinates
    data_dict['isentropic_zonal_angular_momentum'] = np.flip(mlm.Z_isentropic.T,axis=1)
    data_dict['isentropic_zonal_wind'] = np.flip(mlm.u_isentropic.T,axis=1)
    data_dict['isentropic_ertel_potential_vorticity'] = np.flip(mlm.q_isentropic.T,axis=1)
    data_dict['isentropic_density'] = np.flip(mlm.r_isentropic.T,axis=1)
    data_dict['isentropic_pressure'] = np.flip(mlm.p_isentropic.T,axis=1)
    
    ## Surface variables
    data_dict['surface_zonal_angular_momentum'] = mlm.Z_lower_eps
    data_dict['surface_zonal_wind'] = mlm.u_lower_eps
    data_dict['surface_pressure'] = mlm.p_lower_eps
    data_dict['surface_potential_temperature'] = mlm.th_lower_eps
    data_dict['surface_isentropic_density'] = mlm.r_lower_eps

    ## Physical and simulation parameters needed for reconstructing Laguerre diagrams
    data_dict['ld_slims'] = np.array([mlm.sp.smin,mlm.sp.smax])
    data_dict['ld_plims'] = np.array([mlm.sp.pmin,mlm.pp.p00])
    
    ## Seeds and weights for Laguerre diagram
    data_dict['ld_seeds'] = input_data.y
    data_dict['ld_duals'] = solv.ld.duals
    
    ds = atmosphere_bgs.data_io.save_bgs_to_netCDF(data_dict,file_path=save_name+'.nc')