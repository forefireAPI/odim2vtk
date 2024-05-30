#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:39:13 2024

@author: filippi_j
"""
 

import xarray as xr
import numpy as np
import pandas as pd



def create_and_save_gridded_dataset(ds, bounds, startdate, endDate, timedelta, integrationdelta, outputFileName, resolutions=(0.01, 0.01, 100)):
    lon_min, lon_max, lat_min, lat_max = bounds
    res_lat, res_lon, res_alt = resolutions

    # Convert startdate and endDate to datetime64
    startdate = np.datetime64(startdate)
    endDate = np.datetime64(endDate)

    # Create the desired grid
    alt_min, alt_max = 0, 20000  # altitude in meters
    lat_bins = np.arange(lat_min, lat_max, res_lat)
    lon_bins = np.arange(lon_min, lon_max, res_lon)
    alt_bins = np.arange(alt_min, alt_max + res_alt, res_alt)  # +res_alt to include the top value

    # Define the dimensions and coordinates
    lat_coords = (lat_bins[:-1] + lat_bins[1:]) / 2
    lon_coords = (lon_bins[:-1] + lon_bins[1:]) / 2
    alt_coords = (alt_bins[:-1] + alt_bins[1:]) / 2

    # Time-related calculations
    time_coords = np.arange(startdate, endDate, np.timedelta64(timedelta, 's'))
    n_time_steps = len(time_coords)

    # Initialize an empty array for the 4D grid
    grid_shape = (n_time_steps, len(alt_coords), len(lat_coords), len(lon_coords))
    print("data is ",grid_shape)
    grid_array = np.zeros(grid_shape, dtype=np.int16)  

    # Ensure the time variable is in datetime64 format
    ds['time'] = ds['time'].astype('datetime64[ns]')

    # Bin the data into the 4D grid
    for t_index, current_time in enumerate(time_coords):
        print(f"Handling '{current_time}'")
        
        integration_end_time = current_time + np.timedelta64(integrationdelta, 's')
        time_filtered_ds = ds.where((ds['time'] >= current_time) & (ds['time'] < integration_end_time), drop=True)
        
        # Create a pandas DataFrame from the filtered dataset
        df = time_filtered_ds.to_dataframe()

        # Drop rows with NaN values in critical columns
        df = df.dropna(subset=['lat', 'lon', 'alt_m'])

        # Digitize lat, lon, and alt into bins
        df['lat_bin'] = np.digitize(df['lat'], lat_bins) - 1
        df['lon_bin'] = np.digitize(df['lon'], lon_bins) - 1
        df['alt_bin'] = np.digitize(df['alt_m'], alt_bins) - 1

        # Filter out-of-bounds values (those that fall outside the bins)
        df = df[(df['lat_bin'] >= 0) & (df['lat_bin'] < len(lat_bins) - 1)]
        df = df[(df['lon_bin'] >= 0) & (df['lon_bin'] < len(lon_bins) - 1)]
        df = df[(df['alt_bin'] >= 0) & (df['alt_bin'] < len(alt_bins) - 1)]

        # Assign values to the 4D grid (e.g., count of points in each bin)
        for _, row in df.iterrows():
            lat_idx = row['lat_bin']
            lon_idx = row['lon_bin']
            alt_idx = row['alt_bin']
            grid_array[t_index, alt_idx, lat_idx, lon_idx] += 1

    # Create the xarray.DataArray for the 4D data
    data_array_4d = xr.DataArray(grid_array, coords=[time_coords, alt_coords, lat_coords, lon_coords], dims=['time',  'alt', 'lat', 'lon'])

    # Aggregate data for all altitudes into a 3D array for each time step
    data_array_3d = data_array_4d.sum(dim='alt')

    # Create the xarray.Dataset
    result_ds = xr.Dataset({
        'point_count_4d': data_array_4d,
        'point_count_3d': data_array_3d
    })

    # Save the result as a compressed NetCDF file with int8 storage
    encoding = {
        'point_count_4d': {'dtype': 'int16', 'zlib': True},
        'point_count_3d': {'dtype': 'int16', 'zlib': True}
    }
    result_ds.to_netcdf(outputFileName, encoding=encoding)

    print(f"Dataset saved as '{outputFileName}'")

import numpy as np
import xarray as xr
import pandas as pd

from pyevtk.hl import gridToVTK

def create_empty_gridded_dataset(bounds, resolutions=(0.01, 0.01, 100), alt_min=0, alt_max=20000):
    lon_min, lon_max, lat_min, lat_max = bounds
    res_lat, res_lon, res_alt = resolutions

    # Create the desired grid
    lat_bins = np.arange(lat_min, lat_max, res_lat)
    lon_bins = np.arange(lon_min, lon_max, res_lon)
    alt_bins = np.arange(alt_min, alt_max + res_alt, res_alt)  # +res_alt to include the top value

    # Define the dimensions and coordinates
    lat_coords = (lat_bins[:-1] + lat_bins[1:]) / 2
    lon_coords = (lon_bins[:-1] + lon_bins[1:]) / 2
    alt_coords = (alt_bins[:-1] + alt_bins[1:]) / 2

    # Initialize an empty array for the 4D grid
    grid_shape = (len(alt_coords), len(lat_coords), len(lon_coords))
    grid_array = np.zeros(grid_shape, dtype=np.int64)

    # Create the xarray.DataArray for the 4D data
    data_array_3d = xr.DataArray(grid_array, coords=[alt_coords, lat_coords, lon_coords], dims=['alt', 'lat', 'lon'])
    data_array_2d = data_array_3d.sum(dim='alt')

    # Create the xarray.Dataset
    data_array_lightning = xr.Dataset({
        'point_count_3d': data_array_3d,
        'point_count_surface': data_array_2d
    })

    return data_array_lightning, lat_bins, lon_bins, alt_bins

def integrate_filepoints_in_dataset(ncfilepath, data_array_lightning, lat_bins, lon_bins, alt_bins):
    ds = xr.open_dataset(ncfilepath)
   

    # Convert the dataset to a DataFrame
    df = ds.to_dataframe().reset_index() 

    # Drop rows with NaN values in critical columns
    df = df.dropna(subset=['lat', 'lon', 'alt_m'])

    # Digitize lat, lon, and alt into bins
    lat_idx = np.digitize(df['lat'], lat_bins) - 1
    lon_idx = np.digitize(df['lon'], lon_bins) - 1
    alt_idx = np.digitize(df['alt_m'], alt_bins) - 1

    # Ensure the bins are integers
    lat_idx = lat_idx.astype(int)
    lon_idx = lon_idx.astype(int)
    alt_idx = alt_idx.astype(int)

    # Filter out-of-bounds values (those that fall outside the bins)
    valid = (lat_idx >= 0) & (lat_idx < len(lat_bins) - 1) & \
            (lon_idx >= 0) & (lon_idx < len(lon_bins) - 1) & \
            (alt_idx >= 0) & (alt_idx < len(alt_bins) - 1)

    lat_idx = lat_idx[valid]
    lon_idx = lon_idx[valid]
    alt_idx = alt_idx[valid]

    total_rows = len(lat_idx)
    print("Read ", ncfilepath, "Number of valid rows:", total_rows)

    # Use NumPy to create a 4D histogram
    np.add.at(data_array_lightning['point_count_3d'].values, (alt_idx, lat_idx, lon_idx), 1)

    # Update the 3D data array after processing
    data_array_lightning['point_count_surface'][:] = data_array_lightning['point_count_3d'].sum(dim='alt')

def save_to_compressed_nc(data_array_lightning, outputFileName):
    # Save the result as a compressed NetCDF file with int16 storage
    encoding = {
        'point_count_3d': {'dtype': 'int64', 'zlib': True},
        'point_count_surface': {'dtype': 'int64', 'zlib': True}
    }
    data_array_lightning.to_netcdf(outputFileName, encoding=encoding)

    print(f"Dataset saved as '{outputFileName}'")

from pyevtk.hl import gridToVTK
def save_to_vtp(data_array_lightning, outputFileName):
    # Extract the data and coordinates
    point_count_3d = data_array_lightning['point_count_3d'].values
    alt_coords = data_array_lightning.coords['alt'].values
    lat_coords = data_array_lightning.coords['lat'].values
    lon_coords = data_array_lightning.coords['lon'].values

    # Create the VTK file using pyevtk
    gridToVTK(outputFileName, lon_coords, lat_coords, alt_coords, pointData={"point_count": point_count_3d})

    print(f"Dataset saved as VTP file '{outputFileName}'")

# Example usage:

def read_all_lighth():
    bounds = (5, 13, 38, 47)
    resolutions = (0.01, 0.01, 100)
    alt_min = 0
    alt_max = 20000
    data_array_lightning, lat_bins, lon_bins, alt_bins = create_empty_gridded_dataset(bounds, resolutions, alt_min, alt_max)
    
    import os
    
    root_directory = '/Users/filippi_j/Volumes/dataorsu/LAERO/lma/stageMIA2024/STAGE/dataNETCDF4'
    
    # Example usage:
    bounds = (5, 13, 38, 47)
    resolutions = (0.01, 0.01, 100)
    alt_min = 0
    alt_max = 20000
    
    data_array_lightning, lat_bins, lon_bins, alt_bins = create_empty_gridded_dataset(bounds, resolutions, alt_min, alt_max)
    
    for subdir, _, files in os.walk(root_directory):
        for file in files:
            if file.endswith('.nc'):
                ncfilepath = os.path.join(subdir, file)
                integrate_filepoints_in_dataset(ncfilepath, data_array_lightning, lat_bins, lon_bins, alt_bins)
    
    
    #save_to_vtp(data_array_lightning, 'outputXF')
    save_to_compressed_nc(data_array_lightning, 'allSaetaInAFile.nc')
    
from pyevtk.hl import gridToVTK
import numpy as np

def save_to_vts(data_array_lightning, outputFileName):
    # Extract the data and coordinates
    point_count_3d = data_array_lightning['point_count_3d'].values
    alt_coords = data_array_lightning.coords['alt'].values
    lat_coords = data_array_lightning.coords['lat'].values
    lon_coords = data_array_lightning.coords['lon'].values

    # Ensure the coordinates are 1D arrays
    if len(alt_coords.shape) != 1 or len(lat_coords.shape) != 1 or len(lon_coords.shape) != 1:
        raise ValueError("Coordinates should be 1D arrays")

    # Create a meshgrid for the coordinates
    lon_grid, lat_grid, alt_grid = np.meshgrid(lon_coords, lat_coords, alt_coords, indexing='ij')

    # Create the VTK file using pyevtk
    gridToVTK(outputFileName, lon_coords, lat_coords, alt_coords, pointData={"point_count_3d": point_count_3d})

    print(f"Dataset saved as VTS file '{outputFileName}'")

ds = xr.open_dataset("allSaetaInAFile.nc")
save_to_vts(ds, "output")

# Example usage:
# ds = xr.open_dataset("/Users/filippi_j/data/2024/rain18082022/20220818.nc")
# create_and_save_gridded_dataset(
#     ds,
#     bounds=(5, 11, 40, 44),
#     startdate='2022-08-18T02:00:00',
#     endDate='2022-08-18T10:00:00',
#     timedelta=60, 
#     integrationdelta=180,  
#     outputFileName='outputTAS2.nc',
#     resolutions=(0.01, 0.01, 100)
# )

