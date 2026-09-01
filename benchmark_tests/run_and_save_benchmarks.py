import numpy as np
import atmosphere_bgs

# This script runs all benchmark tests for the optimal transport method and saves the outputs in netCDF format.
# This is done for both the smoothed optimal transport solution and for the raw semi-discrete solution.
# The corresponding outputs of the ELIPVI method are also parsed and saved in the same format after interpolation
# onto the same grids. Solutions are represented in latitude and pressure coordinates, as well as latitude and
# potential temperature coordinates.

################ Benchmark tests for optimal transport method ##################
# Run all benchmark tests and save results (both smoothed solution and semi-discrete optimal transport solution)

# 0. Define the simulation parameters
ot_tol = 1e-4       # tolerance for optimal transport problem
eps = None          # smoothing parameter, if None then set depending on data type later on
res = [2**8,2**8]   # resolution of gridded data in latitude and log(pressure)

# 1. Load the input data and define where to save data
print('Loading input data')

data_types = ['lc1low','lc2low','ERA5']
dates = ['200101','200101','201001']
dtimes =['0000','0000','2218'] 

filepath = './input_data/'

for data_type, date, dtime in zip(data_types,dates,dtimes):

    save_name = './output/'+data_type+'_'+date+dtime+'_OT'

    if 'lc' in data_type:
        filename = 'bs_'+data_type+date+dtime
        p00 = 1e+5
        pmin = 10
        
    elif data_type == 'ERA5':
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
    
    # 5. Define data dictionary and save to netCDF
    data_dict = atmosphere_bgs.data_io.get_mlm_data_dict(mlm)
    ds = atmosphere_bgs.data_io.save_bgs_to_netCDF(data_dict,file_path=save_name+'.nc')

################## Parse lifecycle zonal average solutions ##################

# Define the file paths
filepath = './input_data/'
filenames = ['ZM_IGCM_LC1_day0.txt','ZM_IGCM_LC2_day0.txt']

# Parse and format LC1 and LC2 data
for i, filename in enumerate(filenames):
    input_path = filepath+filename
    experiment_type = 'lc' + str(i+1) + '_zonal_average'
    output_path = './output/' + experiment_type + '.nc'
    atmosphere_bgs.data_io.save_IGCM_zonal_mean_to_netCDF(input_path=input_path,experiment_type=experiment_type,output_path=output_path)
    
################## Parse lifecycle ELIPVI output ##################
lc_date = 2001010000
for i in np.arange(2)+1:
    experiment_type = 'lc'+str(i)+'_ELIPVI_'+str(lc_date)
    input_file_path = './input_data/GRID_LC'+str(i)+'_'+str(lc_date)
    output_file_path = './output/'+experiment_type+'.nc'
    data_dict = atmosphere_bgs.data_io.read_ELIPVI_output(input_file_path)
    ds = atmosphere_bgs.data_io.save_bgs_to_netCDF(data_dict,experiment_type=experiment_type,file_path=output_file_path)
    
################## Parse ELIPVI output for ERA5 data ##################
# Define file paths
idatadate = dates[-1]+dtimes[-1]
output_file_path = './output/ELIPVI_output_'+str(idatadate)+'.nc'
input_file_path = './input_data/GRID_NHANMW_'+str(idatadate)
data_dict = atmosphere_bgs.data_io.read_ELIPVI_output(input_file_path)
ds = atmosphere_bgs.data_io.save_bgs_to_netCDF(data_dict,experiment_type='ERA5_ELIPVI_'+str(idatadate),file_path=output_file_path)
