import numpy as np
import xarray as xr
import os
import _atmosphere_bgs

def save_bgs_to_netCDF(data_dict,experiment_type=None,file_path=None):
    """Creates an xarray.Dataset for the background state.

    All numpy data arrays are passed as explicit inputs.
 
    file_path : str
        Optional path to save the dataset as a NetCDF file. If None, the dataset is returned without saving.
    """
    if experiment_type is None:
        experiment_type = 'unspecified'
        
    # Helper function to fall back to zeros if an input array is missing
    def _get_data(arr, shape=0):
        return np.nan*np.zeros(shape) if arr is None else arr
    
    # -------------------------------------------------------------------------
    # 0. Extract coordinate levels and Laguerre cell seed indices from the dictionary (Defaults to None if missing)
    # -------------------------------------------------------------------------
    latitude_levels = data_dict.get("latitude_levels")
    pressure_levels = data_dict.get("pressure_levels")
    potential_temperature_levels = data_dict.get("potential_temperature_levels")
    eta_levels = data_dict.get("eta_levels")
    ld_duals = data_dict.get("ld_duals")

    # Calculate sizes from vector lengths (safely defaults to None if vector is missing)
    n_lat = len(latitude_levels) if latitude_levels is not None else 0
    n_pres = len(pressure_levels) if pressure_levels is not None else 0
    n_theta = len(potential_temperature_levels) if potential_temperature_levels is not None else 0
    n_eta = len(eta_levels) if eta_levels is not None else 0
    n_duals = len(ld_duals) if ld_duals is not None else 0

    # -------------------------------------------------------------------------
    # 1. Define Coordinates Dictionary
    # -------------------------------------------------------------------------
    coords_dict = dict()
    if latitude_levels is not None:
        coords_dict["latitude"] = latitude_levels
    if pressure_levels is not None:
        coords_dict["pressure"] = pressure_levels
    if potential_temperature_levels is not None:
        coords_dict["potential_temperature"] = potential_temperature_levels
    if eta_levels is not None:
        coords_dict["eta"] = eta_levels
    if ld_duals is not None:
        coords_dict["ld_seed_index"] = np.arange(n_duals-1)
        coords_dict["ld_dual_index"] = np.arange(n_duals)
        coords_dict["ld_lim_index"] = np.arange(2)

    # -------------------------------------------------------------------------
    # 2. Define Data Variables
    # -------------------------------------------------------------------------
    data_vars_dict = {
        # Zonal Angular Momentum
        "zonal_angular_momentum": (
            ("latitude", "pressure"),
            _get_data(data_dict.get("zonal_angular_momentum"), shape=(n_lat, n_pres)),
            {"units": "m2 s-1"},
        ),
        "isentropic_zonal_angular_momentum": (
            ("latitude", "potential_temperature"),
            _get_data(data_dict.get("isentropic_zonal_angular_momentum"), shape=(n_lat, n_theta)),
            {"units": "m2 s-1"},
        ),
        "eta_level_zonal_angular_momentum": (
                    ("latitude", "eta"),
                    _get_data(data_dict.get("eta_level_zonal_angular_momentum"), shape=(n_lat, n_eta)),
                    {"units": "m2 s-1"},
                ),
        "surface_zonal_angular_momentum": (
            ("latitude",),
            _get_data(data_dict.get("surface_zonal_angular_momentum"), shape=(n_lat,)),
            {"units": "m2 s-1"},
        ),
        # Zonal Wind
        "zonal_wind": (
            ("latitude", "pressure"),
            _get_data(data_dict.get("zonal_wind"), shape=(n_lat, n_pres)),
            {"units": "m s-1"},
        ),
        "isentropic_zonal_wind": (
            ("latitude", "potential_temperature"),
            _get_data(data_dict.get("isentropic_zonal_wind"), shape=(n_lat, n_theta)),
            {"units": "m s-1"},
        ),
        "eta_level_zonal_wind": (
            ("latitude", "eta"),
            _get_data(data_dict.get("eta_level_zonal_wind"), shape=(n_lat, n_eta)),
            {"units": "m s-1"},
        ),
        "surface_zonal_wind": (
            ("latitude",),
            _get_data(data_dict.get("surface_zonal_wind"), shape=(n_lat,)),
            {"units": "m s-1"},
        ),
        # Potential temperature
        "potential_temperature_field": (
            ("latitude", "pressure"),
            _get_data(data_dict.get("potential_temperature"), shape=(n_lat, n_pres)),
            {"units": "Kelvin"},
        ),
        "eta_level_potential_temperature_field": (
            ("latitude", "eta"),
            _get_data(data_dict.get("eta_level_potential_temperature"), shape=(n_lat, n_eta)),
            {"units": "Kelvin"},
        ),
        "surface_potential_temperature": (
            ("latitude",),
            _get_data(data_dict.get("surface_potential_temperature"), shape=(n_lat,)),
            {"units": "Kelvin"},
        ),
        # Potential Vorticity
        "ertel_potential_vorticity": (
            ("latitude", "pressure"),
            _get_data(data_dict.get("ertel_potential_vorticity"), shape=(n_lat, n_pres)),
            {"units": "PVU"},
        ),
        "isentropic_ertel_potential_vorticity": (
            ("latitude", "potential_temperature"),
            _get_data(data_dict.get("isentropic_ertel_potential_vorticity"), shape=(n_lat, n_theta)),
            {"units": "PVU"},
        ),
        "eta_level_ertel_potential_vorticity": (
            ("latitude", "eta"),
            _get_data(data_dict.get("eta_level_ertel_potential_vorticity"), shape=(n_lat, n_eta)),
            {"units": "PVU"},
        ),
        # Isentropic Density
        "isentropic_density": (
            ("latitude", "potential_temperature"),
            _get_data(data_dict.get("isentropic_density"), shape=(n_lat, n_theta)),
            {"units": "kg m-3"},
        ),
        "surface_isentropic_density": (
            ("latitude",),
            _get_data(data_dict.get("surface_isentropic_density"), shape=(n_lat,)),
            {"units": "kg m-3"},
        ),
        # Geopotential
        "geopotential": (
            ("latitude", "pressure"),
            _get_data(data_dict.get("geopotential"), shape=(n_lat, n_pres)),
            {"units": "m2 s-2"},
        ),
        "isentropic_geopotential": (
            ("latitude", "potential_temperature"),
            _get_data(data_dict.get("isentropic_geopotential"), shape=(n_lat, n_theta)),
            {"units": "m2 s-2"},
        ),
        # Pressure
        "isentropic_pressure": (
            ("latitude", "potential_temperature"),
            _get_data(data_dict.get("isentropic_pressure"), shape=(n_lat, n_theta)),
            {"units": "Pascals"},
        ),
        "surface_pressure": (
            ("latitude",),
            _get_data(data_dict.get("surface_pressure"), shape=(n_lat,)),
            {"units": "Pascals"},
        ),
        "eta_level_pressure": (
            ("latitude", "eta"),
            _get_data(data_dict.get("eta_level_pressure"), shape=(n_lat, n_eta)),
            {"units": "Pascals"},
        ),
        # Variables for reconstructing Laguerre diagram
        "ld_duals": (
                    ("ld_dual_index",),
                    _get_data(data_dict.get("ld_duals")),
                    {"units": "m2 s-2"}
                ),
        "ld_zonal_angular_momentum": (
            ("ld_seed_index",),
            _get_data(data_dict.get("ld_seeds"))[:,0],
            {"units": "m2 s-1"}
        ),
        "ld_potential_temperature": (
                    ("ld_seed_index",),
                    _get_data(data_dict.get("ld_seeds"))[:,1],
                    {"units": "K"}
                ),
        "ld_slims": (
            ("ld_lim_index",),
            _get_data(data_dict.get("ld_slims")),
            {"units": "dimensionless"}
        ),
        "ld_plims": (
            ("ld_lim_index",),
            _get_data(data_dict.get("ld_plims")),
            {"units": "Pascals"}
        )
    }
    
    if data_dict.get("ld_ertel_potential_vorticity") is not None:
        data_vars_dict["ld_ertel_potential_vorticity"] = (
                                                            ("ld_seed_index",),
                                                            _get_data(data_dict.get("ld_ertel_potential_vorticity")),
                                                            {"units": "PVU"}
                                                        )
    
    # Delete data variables that have missing dimensions (i.e., if the corresponding coordinate is None)
    data_vars_dict_filtered = dict()
    for var_name, var_tuple in data_vars_dict.items():
        if all(dim in coords_dict for dim in var_tuple[0]):
            data_vars_dict_filtered[var_name] = var_tuple

    # -------------------------------------------------------------------------
    # 3. Assemble Dataset
    # -------------------------------------------------------------------------
    ds = xr.Dataset(data_vars=data_vars_dict_filtered, coords=coords_dict)
    ds.attrs["experiment_type"] = experiment_type
    
    if file_path is not None:
        # Create parent directories if they do not exist
        output_dir = os.path.dirname(file_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Export to netCDF
        ds.to_netcdf(file_path)
        print(f"Successfully saved dataset to: {file_path}")
    else:
        print("No file path specified. Dataset generated in memory only.")

    return ds


def reconstruct_laguerre_diagram(filepath):
    '''Function to reconstruct the optimal Laguerre diagram from the saved netCDF file.'''
    # Open the netCDF file
    ds = xr.open_dataset(filepath)
    
    # Get default physical and simulation parameters
    pp = _atmosphere_bgs.PhysicalParameters()
    sp = _atmosphere_bgs.SimulationParameters()
    
    # Overwrite domain boundary parameters
    sp.smin, sp.smax = ds['ld_slims'].values
    sp.pmin = ds['ld_plims'].values[0]
    pp.p00 = ds['ld_plims'].values[1]
    
    # Create a new instance of the LaguerreDiagram class
    y = np.vstack([ds['ld_zonal_angular_momentum'],ds['ld_potential_temperature']]).T
    ld = _atmosphere_bgs.LaguerreDiagram(y,ds['ld_duals'].values,pp,sp)
    #ld_epv = ds['ld_ertel_potential_vorticity']
    
    return ld #, ld_epv


def get_u(Z,s,pp):
    '''Compute zonal wind for zonal angular momentum Z and sin-of-lat s'''
    Omega = pp.Omega
    a = pp.a
    return Z/a/np.sqrt(1-s**2) - Omega*a*np.sqrt(1-s**2)

def save_rasterized_sdot_bgs_to_netCDF(input_path,experiment_type=None,output_path=None,res=[256,256]):
    '''Function to rasterize the semi-discrete optimal transport solution (before smoothing)
    and save the rasterized (i.e. gridded) solution in netCDF format. Grid is regular rectangular
    grtid in coordinates latitude and pseudo-height (-H*log(p/p00))
    
    input_path: string; path of input data containing generators and parameters of Laguerre tessellation
    experiment_type: string; name for experiment
    output_path: string; path on which to save output
    res: list; length 2 list specifying resolution of rasterized result 
    '''
    # reconstruct the Laguerre tessellation from the generators and paramters specified by the input path
    ld = reconstruct_laguerre_diagram(input_path)
    
    # define relevant latitude levels and pressure levels using pseudo-height coordinate
    H = 6500
    p00 = ld.phys.p00
    
    s2lat = lambda s : np.rad2deg(np.arcsin(s))
    p2h = lambda p : -H*np.log(p/p00)
    lat2s = lambda lat : np.sin(np.deg2rad(lat))
    h2p = lambda h : p00*np.exp(-h/H)

    transform = lambda x : [s2lat(x[0]), p2h(x[1])]
    rast = ld.get_rasterizer(transform)
            
    smin = ld.sim.smin; smax = ld.sim.smax
    pmin = ld.sim.pmin; pmax = ld.sim.pmax
    lg = np.linspace(s2lat(smin),s2lat(smax),res[0]+1)
    hg = np.linspace(p2h(pmin),p2h(pmax),res[1]+1)
    
    latitude_levels = lg[:-1] + np.diff(lg)/2 # midpoints
    pressure_levels = h2p(hg[:-1] + np.diff(hg)/2) # midpoints
    
    # get rasterised zonal angular momentum and potential temperature
    def rasterize_centroidal_values(val):
        rv = rast.rasterize(val, res).copy()
        fill = rast.fill
        rv = np.where(fill > 0.5, np.divide(rv, fill, where=(fill > 0.5)), np.nan)
        return rv

    zonal_angular_momentum = np.flip(rasterize_centroidal_values(ld.ys[:,0]),axis=1)
    potential_temperature = np.flip(rasterize_centroidal_values(ld.ys[:,1]),axis=1)

    # get rasterised zonal wind from zonal angular momentum
    lg, _ = np.meshgrid(latitude_levels,np.flip(p2h(pressure_levels)))
    sg = lat2s(lg).T
    zonal_wind = get_u(zonal_angular_momentum,sg,ld.phys)
    #ertel_potential_vorticity = rast.rast(ld_epv,res)             ##### need cellwise epv saved to netCDF file along with seeds and weights

    # extract surface pressure and interpolate onto latitude levels
    n = len(ld.ys)
    edges = [e for e in ld.edglist if e.pj >= n+6]
    edge_pts, edge_idx = np.unique(np.vstack([e.ls.x for e in edges]),axis=0,return_index=True)
    idx = np.argsort(edge_pts[:,0])
    s_surf, p_surf = edge_pts[idx].T # Laguerre cell surface pressure
    lat_surf = s2lat(s_surf)
    surface_pressure = np.interp(latitude_levels,lat_surf,p_surf)
    
    # get surface potential temperature and zonal wind and evaluate on grid
    surf_cell_idx = np.hstack([np.ones(len(e.ls.x),dtype=int)*(e.pi-6) for e in edges])[edge_idx][idx]
    zam_surf, th_surf = ld.ys[surf_cell_idx].T
    surface_zonal_angular_momentum = np.interp(latitude_levels,lat_surf,zam_surf)
    surface_potential_temperature = np.interp(latitude_levels,lat_surf,th_surf)
    surface_zonal_wind = get_u(surface_zonal_angular_momentum,lat2s(latitude_levels),ld.phys)

    ##### need pv values on cells to be defined or need to do interpolation of pv as function of circulation so need input data...
    ##### probably best to save values of pv on cells as obtained by interpolation of pv from circulation at point of solution

    data_dict = {
        'latitude_levels' : latitude_levels,
        'pressure_levels' : pressure_levels,
        'zonal_angular_momentum' : zonal_angular_momentum,
        'potential_temperature' : potential_temperature,
        'zonal_wind' : zonal_wind,
        #'ertel_potential_vorticity' : ertel_potential_vorticity,
        'surface_pressure' : surface_pressure,
        'surface_potential_temperature' : surface_potential_temperature,
        'surface_zonal_wind' : surface_zonal_wind,
        'ld_duals' : ld.duals,
        'ld_seeds' : ld.ys,
        'ld_slims' : np.array([ld.sim.smin,ld.sim.smax]),
        'ld_plims' : np.array([ld.sim.pmin,ld.phys.p00])
        }
    
    ds = save_bgs_to_netCDF(data_dict,experiment_type=experiment_type,file_path=output_path)
    return ds