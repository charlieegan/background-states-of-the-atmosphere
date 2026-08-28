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