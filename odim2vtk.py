import numpy as np
import xarray as xr
import vtk
import os
import glob
from vtk.util.numpy_support import numpy_to_vtk


def polar_to_cartesian(azimuth, elevation, range_data):
    azimuth_rad = np.deg2rad(azimuth)
    elevation_rad = np.deg2rad(elevation)
    x = range_data * np.sin(azimuth_rad) * np.cos(elevation_rad)
    y = range_data * np.cos(azimuth_rad) * np.cos(elevation_rad)
    z = range_data * np.sin(elevation_rad)
    return x, y, z

def write_to_bin(file_path, array, mode='wb'):
    array = np.ascontiguousarray(array, dtype=np.float32)
    print("adding ", len(array), " elements to ", file_path)
    with open(file_path, mode) as file:
        array.tofile(file)

def odimNC2bin(ncinpath_time_sorted_array, binoutpath, keys=['RHOHV', 'TH', 'DBZH', 'VRADH', 'WRADH']):
    data_keys = {key: [] for key in keys}  # Dictionary to hold data for each key

    # Load the first file to get the initial time
    with xr.open_dataset(ncinpath_time_sorted_array[0]) as data:
        first_time = data['time'].values[0]

    # Process each file
    for ncpath in ncinpath_time_sorted_array:
        print("Reading", ncpath)
        with xr.open_dataset(ncpath) as data:
            all_x, all_y, all_z, all_time_diffs = [], [], [], []
            
            for i in range(len(data['time'])):
                azimuth = data['azimuth'][i].values
                elevation = data['elevation'][i].values
                range_data = data['range'].values
                current_time = data['time'][i].values
                
                x, y, z = polar_to_cartesian(azimuth, elevation, range_data)
                dbzh = data['DBZH'][i].values.flatten()
                valid_indices = dbzh > 0  # Filter condition
                
                all_x.extend(x[valid_indices])
                all_y.extend(y[valid_indices])
                all_z.extend(z[valid_indices])
                
                time_diff = (current_time - first_time) / np.timedelta64(1, 's')
                # Append the time_diff for each valid index, but only once per valid index
                all_time_diffs.extend([time_diff] * len(x[valid_indices]))
          
                
                for key in keys:
                    data_keys[key].extend(data[key][i].values.flatten()[valid_indices])
                    
            # Append all collected data to binary files
            write_to_bin(f"{binoutpath}/all_x.bin", all_x, 'ab')
            write_to_bin(f"{binoutpath}/all_y.bin", all_y, 'ab')
            write_to_bin(f"{binoutpath}/all_z.bin", all_z, 'ab')
            write_to_bin(f"{binoutpath}/all_time.bin", np.array(all_time_diffs, dtype=np.float32), 'ab')
            
            for key in keys:
                write_to_bin(f"{binoutpath}/all_{key}.bin", data_keys[key], 'ab')

            # Clear lists after writing to file
            all_x, all_y, all_z, all_time_diffs = [], [], [], []
            for key in keys:
                data_keys[key] = []

def read_bin(file_path, dtype=np.float32):
    """Reads a binary file into a NumPy array."""
    return np.fromfile(file_path, dtype=dtype)

def rawBin2vtk(binoutpath, vtkoutpath, keys=['RHOHV', 'TH', 'DBZH','VRADH', 'WRADH']):
    # Read binary data for coordinates
    all_x = read_bin(f"{binoutpath}/all_x.bin")
    all_y = read_bin(f"{binoutpath}/all_y.bin")
    all_z = read_bin(f"{binoutpath}/all_z.bin")

    # Stack arrays and convert them to a VTK-compatible format
    flat_positions = np.vstack((all_x, all_y, all_z)).T.ravel()
    vtk_positions = vtk.vtkFloatArray()
    vtk_positions.SetNumberOfComponents(3)
    vtk_positions.SetArray(flat_positions, len(flat_positions), 1)

    particle_points = vtk.vtkPoints()
    particle_points.SetData(vtk_positions)

    particle_polydata = vtk.vtkPolyData()
    particle_polydata.SetPoints(particle_points)

    # Read binary data for time and convert it to a VTK scalar array
    all_time = read_bin(f"{binoutpath}/all_time.bin")
    add_scalar_array(particle_polydata, all_time, "airtime")

    # Read binary data for other keys and convert them to VTK scalar arrays
    for key in keys:
        key_data = read_bin(f"{binoutpath}/all_{key}.bin")
        add_scalar_array(particle_polydata, key_data, key)

    # Write the VTK file
    write_vtk(particle_polydata, vtkoutpath)
    

def add_scalar_array(polydata, data_array, name):
    vtk_scalar_array = vtk.vtkFloatArray()
    vtk_scalar_array.SetName(name)
    vtk_scalar_array.SetArray(np.ascontiguousarray(data_array), len(data_array), 1)
    polydata.GetPointData().AddArray(vtk_scalar_array)

def write_vtk(polydata, file_path):
    particle_writer = vtk.vtkXMLPolyDataWriter()
    particle_writer.SetFileName(file_path)
    particle_writer.SetInputData(polydata)
    particle_writer.SetDataModeToBinary()
    particle_writer.Write()
    
def get_sorted_nc_files(nc_path):
    # Create a pattern to match all .nc files in the specified directory
    pattern = os.path.join(nc_path, '*.nc')
    # Find all files matching the pattern
    nc_files = glob.glob(pattern)
    # Sort files based on the timestamp in the filename
    
    nc_files_sorted = sorted(nc_files, key=lambda x: x.split('_')[-1])
    return nc_files_sorted


def main():
    ncinpath = '/Users/filippi_j/data/2024/ballarat/2_20240222_033000.pvol.nc'
    vtkoutpath = "/Users/filippi_j/data/2024/ballarat/output_radar_data.vtp"
    
    binOutpath = "/Users/filippi_j/data/2024/ballarat/tmp/"
    nc_dir_path = '/Users/filippi_j/data/2024/ballarat/radar/AURA_2_20240222_nc 2/'
    
    all_ncs = get_sorted_nc_files(nc_dir_path)
    selected = all_ncs[:]
    
    print("\n".join(selected))
    odimNC2bin(selected, binOutpath)
    
    rawBin2vtk(binOutpath,vtkoutpath)
   

if __name__ == "__main__":
    main()
