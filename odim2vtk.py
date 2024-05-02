import numpy as np
import xarray as xr
import vtk
from vtk.util.numpy_support import numpy_to_vtk

def polar_to_cartesian(azimuth, elevation, range_data):
    azimuth_rad = np.deg2rad(azimuth)
    elevation_rad = np.deg2rad(elevation)
    x = range_data * np.sin(azimuth_rad) * np.cos(elevation_rad)
    y = range_data * np.cos(azimuth_rad) * np.cos(elevation_rad)
    z = range_data * np.sin(elevation_rad)
    return x, y, z

def odim2vtk(ncinpath, vtkoutpath):
    file_path = ncinpath
    data = xr.open_dataset(file_path)
    first_time = data['time'].values[0]
    all_x, all_y, all_z, all_time_diffs, th_data, rhohv_data, dbzh_data = [], [], [], [], [], [], []

    for i in range(len(data['time'])):
        azimuth = data['azimuth'][i].values
        elevation = data['elevation'][i].values
        range_data = data['range'].values
        current_time = data['time'][i].values
        rhohv = data['RHOHV'][i].values.flatten()
        dbzh = data['DBZH'][i].values.flatten()
        
        x, y, z = polar_to_cartesian(azimuth, elevation, range_data)
        
        valid_indices = np.where(dbzh > 0)[0]
        all_x.extend(x[valid_indices])
        all_y.extend(y[valid_indices])
        all_z.extend(z[valid_indices])
        
        th_data.extend(data['TH'][i].values.flatten()[valid_indices])
        rhohv_data.extend(rhohv[valid_indices])
        dbzh_data.extend(dbzh[valid_indices])
        
        time_diff = (current_time - first_time) / np.timedelta64(1, 's')
        all_time_diffs.extend([time_diff] * len(valid_indices))

    num_particles = len(all_time_diffs)
    particle_points = vtk.vtkPoints()

    # Convert lists to numpy arrays, ensure they are float32 and contiguous
    all_x = np.ascontiguousarray(np.array(all_x, dtype=np.float32))
    all_y = np.ascontiguousarray(np.array(all_y, dtype=np.float32))
    all_z = np.ascontiguousarray(np.array(all_z, dtype=np.float32))
    all_time_diffs = np.ascontiguousarray(np.array(all_time_diffs, dtype=np.float32))
    th_data = np.ascontiguousarray(np.array(th_data, dtype=np.float32))
    rhohv_data = np.ascontiguousarray(np.array(rhohv_data, dtype=np.float32))
    dbzh_data = np.ascontiguousarray(np.array(dbzh_data, dtype=np.float32))
    
    flat_positions = np.ascontiguousarray(np.vstack((all_x, all_y, all_z)).T.ravel(), dtype=np.float32)
    vtk_positions = vtk.vtkFloatArray()
    vtk_positions.SetNumberOfComponents(3)
    vtk_positions.SetArray(flat_positions, len(flat_positions), 1)
    particle_points.SetData(vtk_positions)
    particle_polydata = vtk.vtkPolyData()
    particle_polydata.SetPoints(particle_points)

    # Adding scalar arrays
    add_scalar_array(particle_polydata, all_time_diffs, "airtime")
    add_scalar_array(particle_polydata, rhohv_data, "RHOHV")
    add_scalar_array(particle_polydata, th_data, "TH")
    add_scalar_array(particle_polydata, dbzh_data, "DBZH")

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
    
  
def main():
    ncinpath = '/Users/filippi_j/data/2024/ballarat/2_20240222_033000.pvol.nc'
    vtkoutpath = "/Users/filippi_j/data/2024/ballarat/output_radar_data.vtp"
    odim2vtk(ncinpath, vtkoutpath)
   

if __name__ == "__main__":
    main()
