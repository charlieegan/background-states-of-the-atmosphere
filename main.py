# This script obtains gridded MLM states and saves them as netCDF files.
# It uses the optimal transport method followed by post-processing steps.

import numpy as np
import atmosphere_bgs
import xarray as xr

# 0. Define the simulation parameters
ot_tol = 1e-4       # tolerance for optimal transport problem
eps = None          # smoothing parameter, if None then set depending on data type later on
res = [2**8,2**8]   # resolution of gridded data in latitude and log(pressure)

# 1. Load the input data and define where to save data
print('Loading input data')
# data_type = 'lc1low'
# date = '200101'
# dtime = '0000'

# data_type = 'lc2low'
# date = '200101'
# dtime = '0000'

# data_type = 'ERA5'
# date = '200911'
# dtime = '0100'

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

    surf_mask = mlm.surf_mask.T
    lats = np.unique(mlm.lg)
    plevs = np.unique(mlm.pg)
    thlevs = np.sort(mlm.thlevs)

    dims_lat_p = ("latitude","pressure")
    coords_lat_p = (lats,plevs)

    dims_lat_th = ("latitude","theta")
    coords_lat_th = (lats,thlevs)

    ## Interior variables in latitude/pressure coordinates
    Phi_lat_p = xr.DataArray(surf_mask*mlm.Phi_eps.T,dims=dims_lat_p,coords=coords_lat_p)
    Z_lat_p = xr.DataArray(surf_mask*mlm.Z_eps.T,dims=dims_lat_p,coords=coords_lat_p)
    th_lat_p = xr.DataArray(surf_mask*mlm.th_eps.T,dims=dims_lat_p,coords=coords_lat_p)
    u_lat_p = xr.DataArray(surf_mask*mlm.u_eps.T,dims=dims_lat_p,coords=coords_lat_p)

    ## Interior variables in latitude/theta coordinates
    p_lat_th = xr.DataArray(np.flip(mlm.p_isentropic.T,axis=1),dims=dims_lat_th,coords=coords_lat_th)
    Z_lat_th = xr.DataArray(np.flip(mlm.Z_isentropic.T,axis=1),dims=dims_lat_th,coords=coords_lat_th)
    u_lat_th = xr.DataArray(np.flip(mlm.u_isentropic.T,axis=1),dims=dims_lat_th,coords=coords_lat_th)
    q_lat_th = xr.DataArray(np.flip(mlm.q_isentropic.T,axis=1),dims=dims_lat_th,coords=coords_lat_th)
    r_lat_th = xr.DataArray(np.flip(mlm.r_isentropic.T,axis=1),dims=dims_lat_th,coords=coords_lat_th)

    ## Surface variables
    Z_lower = xr.DataArray(mlm.Z_lower_eps,dims=('latitude',),coords=(lats,))
    p_lower = xr.DataArray(mlm.p_lower_eps,dims=('latitude',),coords=(lats,))
    r_lower = xr.DataArray(mlm.r_lower_eps,dims=('latitude',),coords=(lats,))
    th_lower = xr.DataArray(mlm.th_lower_eps,dims=('latitude',),coords=(lats,))
    u_lower = xr.DataArray(mlm.u_lower_eps,dims=('latitude',),coords=(lats,))
    
    ## Physical and simulation parameters needed for reconstructing Laguerre diagrams
    smin = xr.DataArray(mlm.sp.smin)
    smax = xr.DataArray(mlm.sp.smax)
    pmin = xr.DataArray(mlm.sp.pmin)
    p00 = xr.DataArray(mlm.pp.p00)
    
    ## Seeds and weights for Laguerre diagram
    ld_seeds = xr.DataArray(input_data.y,dims=('Zonal angular momentum','Potential temperature'))
    ld_duals = xr.DataArray(solv.ld.duals,dims=('Seed index',))
    
    mlm_gridded = xr.Dataset(dict(Phi_lat_p=Phi_lat_p,
                                Z_lat_p=Z_lat_p,
                                th_lat_p=th_lat_p,
                                u_lat_p=u_lat_p,
                                p_lat_th=p_lat_th,
                                Z_lat_th=Z_lat_th,
                                u_lat_th=u_lat_th,
                                q_lat_th=q_lat_th,
                                r_lat_th=r_lat_th,
                                Z_lower=Z_lower,
                                p_lower=p_lower,
                                r_lower=r_lower,
                                th_lower=th_lower,
                                u_lower=u_lower,
                                smin=smin,
                                smax=smax,
                                pmin=pmin,
                                p00=p00,
                                ld_seeds=ld_seeds,
                                ld_duals=ld_duals
                                ))

    # 5. Save gridded data in netCDF format
    mlm_gridded.to_netcdf(save_name+'.nc')