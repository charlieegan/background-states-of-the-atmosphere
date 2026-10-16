import numpy as np
import atmosphere_bgs
import os

# This script applies the optimal transport method to one ERA5 dataset and saves the output in netCDF format.

# 0. Define the simulation parameters
ot_tol = 1e-4       # tolerance for optimal transport problem
eps = 5e-3          # smoothing parameter, if None then set depending on data type later on
res = [2**8,2**8]   # resolution of gridded data in latitude and log(pressure)

dates = ['200101','200101','201001']
times = ['0000','0000','2218']
exp_types = ['lc1low','lc2low','era5']

base_dir = os.path.dirname(os.path.abspath(__file__))

for exp_type, date, time in zip(exp_types,dates,times):

    # 1. Load the input data and define where to save data
    print('Loading input data')
    filepath = os.path.join(base_dir, 'input_data') + os.sep
    save_name = os.path.join(base_dir, 'output', 'otbgs_'+exp_type+date+time)
    
    if exp_type == 'era5':
        filename = 'bs_diout'+date+time
        p00 = None
        pmin = 1
    else:    
        filename = 'bs_'+exp_type+date+time
        p00 = 1e+5
        pmin = 10

    input_data = atmosphere_bgs.DataLoader(filepath+filename,p00=p00,load_all=True,pmin=pmin,mode='legacy')

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
    pvmaxth = input_data.data_dict['pvmaxth']   # previously 'MAX PV ON THETA LEVELS'
            
    mlm = atmosphere_bgs.PostProcessor(solv,bscirc,pvlevs,thlevs,pvmaxth,eps=eps,nb_type='Lag_cell_8',res=res,make_surface_adjustments=True)

    # 4. Put gridded data into an xarray Dataset
    print('Saving gridded data in netCDF format')

    data_dict = atmosphere_bgs.data_io.get_mlm_data_dict(mlm)
    ds = atmosphere_bgs.data_io.save_bgs_to_netCDF(data_dict,file_path=save_name+'.nc')