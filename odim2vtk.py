import numpy as np
import xarray as xr
import vtk
import os
import glob
from vtk.util.numpy_support import numpy_to_vtk
from pyevtk.vtk import *  

import numpy as np
import xarray as xr

def _addDataToFile(vtkFile, cellData, pointData):
    # Point data
    if pointData is not None:
        keys = list(pointData.keys())

        vtkFile.openData("Point", scalars = keys[0])
        for key in keys:
            data = pointData[key]
            vtkFile.addData(key, data)
        vtkFile.closeData("Point")

    # Cell data
    if cellData is not None:
        keys = list(cellData.keys())
        vtkFile.openData("Cell", scalars = keys[0])
        for key in keys:
            data = cellData[key]
            vtkFile.addData(key, data)
        vtkFile.closeData("Cell")

def polar_to_cartesian3D(azimuth, elevation, range_data):
    # Ensure inputs are arrays
    azimuth = np.atleast_1d(azimuth)
    elevation = np.atleast_1d(elevation)
    range_data = np.atleast_1d(range_data)

    # Convert azimuth and elevation from degrees to radians and expand dimensions
    azimuth_rad = np.deg2rad(azimuth)[:, np.newaxis, np.newaxis]  # shape (nazimuth, 1, 1)
    elevation_rad = np.deg2rad(elevation)[np.newaxis, :, np.newaxis]  # shape (1, nelevation, 1)
    range_data = range_data[np.newaxis, np.newaxis, :]  # shape (1, 1, nrange)

    # Calculate Cartesian coordinates for each combination of azimuth, elevation, and range
    x = range_data * np.sin(azimuth_rad) * np.cos(elevation_rad)
    y = range_data * np.cos(azimuth_rad) * np.cos(elevation_rad)
    z = range_data * np.sin(elevation_rad)

    return x, y, z

def gridToVTK(path, x, y, z, cellData=None, pointData=None, vectData=None, start=None, end=None):
    if start is None:
        start = (0, 0, 0)

    ftype = VtkStructuredGrid

    s = x.shape
    nx, ny, nz = s[0] - 1, s[1] - 1, s[2] - 1

    if end is None:
        end = (nx, ny, nz)

    w = VtkFile(path, ftype)
    w.openGrid(start=start, end=end)
    w.openPiece(start=start, end=end)

    w.openElement("Points")
    w.addData("points", (x, y, z))
    w.closeElement("Points")

    _addDataToFile(w, cellData, pointData)
    w.closePiece()
    w.closeGrid()
    # Write coordinates
    w.appendData((x, y, z))
    # Write data
    _appendDataToFile(w, cellData, pointData)
    w.save()
    return w.getFileName()


def odim2MeshVTK(ncinpath_time_sorted_array, outvtk_path, filter=None, keys=['RHOHV', 'TH', 'DBZH', 'VRADH', 'WRADH']):
    # Open the first dataset to get the reference time
    with xr.open_dataset(ncinpath_time_sorted_array[0]) as data:
        first_time = data['time'].values[0]

    data_keys = {key: None for key in keys}
    all_x, all_y, all_z = None, None, None
    all_time_diffs = []

    # Process each file
    for ncpath in ncinpath_time_sorted_array:
        with xr.open_dataset(ncpath) as data:
            print("Reading", ncpath)

            for i in range(len(data['time'])):
                azimuth = data['azimuth'][i].values
                elevation = data['elevation'][i].values
                range_data = data['range'].values
                current_time = data['time'][i].values

                # Convert polar to Cartesian coordinates
                x, y, z = polar_to_cartesian3D(azimuth, elevation, range_data)

                if all_x is None:
                    # Initialize the 3D arrays if this is the first iteration
                    all_x, all_y, all_z = x.copy(), y.copy(), z.copy()
                    for key in keys:
                        data_keys[key] = data[key][i].values.copy()
                else:
                    # Append new data to existing arrays
                    all_x = np.concatenate((all_x, x), axis=0)
                    all_y = np.concatenate((all_y, y), axis=0)
                    all_z = np.concatenate((all_z, z), axis=0)
                    for key in keys:
                        data_keys[key] = np.concatenate((data_keys[key], data[key][i].values), axis=0)

                # Calculate time difference in seconds
                time_diff = (current_time - first_time) / np.timedelta64(1, 's')
                time_diff_array = np.full(x.shape, time_diff)
                all_time_diffs.append(time_diff_array)

        # Combine all time differences into a single array
        all_time_diffs = np.concatenate(all_time_diffs, axis=0)
    
        # Prepare and save the VTK file
        # Note: You'll need to ensure gridToVTK can handle these full 3D datasets
        gridToVTK(f"{outvtk_path}/output{int(np.min(all_time_diffs))}.vtk",
                  all_x, all_y, all_z, 
                  pointData={**{key: data_keys[key] for key in keys}, 'time': all_time_diffs})
        
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
        print(first_time)
        return 0
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

def rawBin2vtk(binoutpath, vtkoutpath, keys=['RHOHV', 'TH', 'DBZH','VRADH', 'WRADH'], origin =(0,0,0),timeOrigin=0):
    # Read binary data for coordinates
    all_x = read_bin(f"{binoutpath}/all_x.bin")+origin[0]
    all_y = read_bin(f"{binoutpath}/all_y.bin")+origin[1]
    all_z = read_bin(f"{binoutpath}/all_z.bin")+origin[2]
    print("data read")
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
    all_time = read_bin(f"{binoutpath}/all_time.bin")+timeOrigin
    print("time read",all_time[0])
    add_scalar_array(particle_polydata, all_time, "airtime")

    # Read binary data for other keys and convert them to VTK scalar arrays
    for key in keys:
        key_data = read_bin(f"{binoutpath}/all_{key}.bin")
        print("adding ",key)
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



def transform_dataset(ds):
    # Extraire les valeurs uniques d'élévation
    elevations = np.unique(ds.elevation.data)
    
    # Dictionnaire pour stocker les résultats par élévation
    subsets = {}
    
    # Boucler sur chaque élévation
    for elevation in elevations:
        # Filtrer les données pour une élévation spécifique
        subset = ds.where(ds.elevation == elevation, drop=True)
        
        # Extraire les données 'azimuth' et les assigner comme coordonnée
        azimuth_data = subset.azimuth.data  # Utiliser .data ici pour éviter l'erreur
        subset = subset.assign_coords(azimuth=('time', azimuth_data))
        
        # Échanger la dimension 'time' avec 'azimuth' seulement si tous les 'time' sont uniques par 'azimuth'
        if subset.sizes['time'] == len(np.unique(azimuth_data)):
            subset = subset.swap_dims({'time': 'azimuth'})
        
        # Trier les données par 'azimuth' pour garantir l'ordre croissant
        subset = subset.sortby('azimuth')
        
        # Enregistrer le subset transformé
        subsets[elevation] = subset

    return subsets

import vtk
from vtk.util import numpy_support  

def create_vtk_file(ds, elevation, filename):
    # Preparation of data
    azimuths = ds.azimuth.data
    ranges = ds.range.data
    values = ds.DBZH.data  # Example for 'DBZH' variable, adapt as needed

    # Convert azimuths from -180 to +180 to 0 to 360 degrees
    azimuths = np.mod(azimuths + 360, 360)
    
    # Polar coordinates calculations and elevation conversion to radians
    theta = np.deg2rad(azimuths)
    elevation_rad = np.deg2rad(elevation)

    # Adjust range data for the cosine of elevation
    adjusted_ranges = ranges * np.cos(elevation_rad)

    npts = np.shape(values)
    Points = vtk.vtkPoints()
    Vertices = vtk.vtkCellArray()
    values_vtk = []
    ids = np.zeros((npts[1], npts[0]))

    # Insert points in VTK
    for i_radius in range(npts[1]):
        _r = adjusted_ranges[i_radius]
        _z = ranges[i_radius] * np.sin(elevation_rad)  # Calculate Z based on actual range and sine of elevation
        for i_theta in range(npts[0]):
            _theta = theta[i_theta]
            _x = _r * np.cos(_theta)
            _y = _r * np.sin(_theta)
            _id = Points.InsertNextPoint(_x, _y, _z)
            ids[i_radius, i_theta] = _id
            values_vtk.append(values[i_theta, i_radius])

    values_vtk = np.array(values_vtk)

    # Create connectivity of the cells
    for i_radius in range(npts[1]-1):
        for i_theta in range(npts[0]):
            i_next_theta = (i_theta + 1) % npts[0]  # Wrap around to connect end to start
            Vertices.InsertNextCell(4)
            Vertices.InsertCellPoint(int(ids[i_radius, i_theta]))
            Vertices.InsertCellPoint(int(ids[i_radius+1, i_theta]))
            Vertices.InsertCellPoint(int(ids[i_radius+1, i_next_theta]))
            Vertices.InsertCellPoint(int(ids[i_radius, i_next_theta]))

    # Create the PolyData
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(Points)
    polyData.SetPolys(Vertices)

    # Add data
    pointData = polyData.GetPointData()
    values_vtk2 = numpy_support.numpy_to_vtk(values_vtk)
    values_vtk2.SetName("InvDFT")
    pointData.AddArray(values_vtk2)

    # Write to a file
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polyData)
    writer.Write()
 

def create_vtk_volume_points(subsets, filename):
    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()
    values_vtk = []

    elevations = sorted(subsets.keys())[:-1]  # Sorted elevations

    # Dictionary to keep track of point IDs across elevations
    point_ids = {elevation: {} for elevation in elevations}

    # Generate points and values for all elevations
    for elevation in elevations:
        print(elevation)
        ds = subsets[elevation]
        azimuths = ds.azimuth.data
        ranges = ds.range.data
        values = ds.DBZH.data  # Assuming data variable 'DBZH'

        theta = np.deg2rad(azimuths)
        elevation_rad = np.deg2rad(elevation)
        adjusted_ranges = ranges * np.cos(elevation_rad)

        npts = (len(azimuths), len(ranges))
        ids = np.zeros((npts[1], npts[0]), dtype=int)

        # Insert points for current elevation
        for i_radius in range(npts[1]):
            _r = adjusted_ranges[i_radius]
            _z = ranges[i_radius] * np.sin(elevation_rad)
            for i_theta in range(npts[0]):
                _theta = theta[i_theta]
                _x = _r * np.cos(_theta)
                _y = _r * np.sin(_theta)
                _id = points.InsertNextPoint(_x, _y, _z)
                ids[i_radius, i_theta] = _id
                values_vtk.append(values[i_theta, i_radius])

        point_ids[elevation] = ids

    # Create connectivity within and between elevations
    for elevation_index, elevation in enumerate(elevations[:-1]):
        print(elevation)
        next_elevation = elevations[elevation_index + 1]
        ids_current = point_ids[elevation]
        ids_next = point_ids[next_elevation]
        
        for i_radius in range(npts[1]-1):
            for i_theta in range(npts[0]):
                i_next_theta = (i_theta + 1) % npts[0]
                
                # Connect current level
                polys.InsertNextCell(4)
                polys.InsertCellPoint(ids_current[i_radius, i_theta])
                polys.InsertCellPoint(ids_current[i_radius+1, i_theta])
                polys.InsertCellPoint(ids_current[i_radius+1, i_next_theta])
                polys.InsertCellPoint(ids_current[i_radius, i_next_theta])

                # Correct vertical connections
                polys.InsertNextCell(4)
                polys.InsertCellPoint(ids_current[i_radius, i_theta])
                polys.InsertCellPoint(ids_current[i_radius, i_next_theta])
                polys.InsertCellPoint(ids_next[i_radius, i_next_theta])
                polys.InsertCellPoint(ids_next[i_radius, i_theta])

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetPolys(polys)

    # Add data
    values_vtk_array = numpy_support.numpy_to_vtk(values_vtk, deep=True, array_type=vtk.VTK_FLOAT)
    values_vtk_array.SetName("DBZH")
    polyData.GetPointData().AddArray(values_vtk_array)

    # Write to a file
    writer = vtk.vtkXMLPolyDataWriter()  # For .vtp file, or use vtkXMLUnstructuredGridWriter for .vtu file
    writer.SetFileName(filename)
    writer.SetInputData(polyData)
    writer.Write()
 

def create_vtk_full_volume(subsets, filename):
    points = vtk.vtkPoints()
    ugrid = vtk.vtkUnstructuredGrid()

    elevations = sorted(subsets.keys())[:-1]  # Sorted elevations
    values_vtk = []

    # Generate points for all elevations and store their point IDs
    point_ids = {}  # To hold point IDs for each elevation and position
    azimuths = subsets[elevations[0]].azimuth.data
    ranges = subsets[elevations[0]].range.data
    for elevation in elevations:
        print(elevation)
        ds = subsets[elevation]
        values = ds.DBZH.data

        theta = np.deg2rad(azimuths)
        elevation_rad = np.deg2rad(elevation)
        adjusted_ranges = ranges * np.cos(elevation_rad)

        npts = (len(azimuths), len(ranges))
        ids = np.zeros((npts[1], npts[0]), dtype=int)

        for i_radius in range(npts[1]):
            _r = adjusted_ranges[i_radius]
            _z = ranges[i_radius] * np.sin(elevation_rad)
            for i_theta in range(npts[0]):
                _theta = theta[i_theta]
                _x = _r * np.cos(_theta)
                _y = _r * np.sin(_theta)
                _id = points.InsertNextPoint(_x, _y, _z)
                ids[i_radius, i_theta] = _id
                values_vtk.append(values[i_theta, i_radius])

        point_ids[elevation] = ids

    # Create hexahedrons linking points from consecutive elevations
    for e_index, elevation in enumerate(elevations[:-1]):
        print(elevation)
        next_elevation = elevations[e_index + 1]
        for i_radius in range(npts[1] - 1):
            for i_theta in range(npts[0] - 1):
                hex_cell = vtk.vtkHexahedron()
                for k, (ei, ri, ti) in enumerate([
                    (elevation, i_radius, i_theta),
                    (elevation, i_radius+1, i_theta),
                    (elevation, i_radius+1, i_theta+1),
                    (elevation, i_radius, i_theta+1),
                    (next_elevation, i_radius, i_theta),
                    (next_elevation, i_radius+1, i_theta),
                    (next_elevation, i_radius+1, i_theta+1),
                    (next_elevation, i_radius, i_theta+1)
                ]):
                    hex_cell.GetPointIds().SetId(k, point_ids[ei][ri, ti])
                ugrid.InsertNextCell(hex_cell.GetCellType(), hex_cell.GetPointIds())

    ugrid.SetPoints(points)

    # Convert values to VTK array and add to grid
    values_vtk_array = numpy_support.numpy_to_vtk(np.array(values_vtk), deep=True, array_type=vtk.VTK_FLOAT)
    values_vtk_array.SetName("DBZH")
    ugrid.GetPointData().AddArray(values_vtk_array)

    # Write data to file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(ugrid)
    writer.SetDataModeToBinary()
    writer.Write()
 

def create_vtk_filtered_volume(subsets, filename):
    points = vtk.vtkPoints()
    ugrid = vtk.vtkUnstructuredGrid()

    elevations = sorted(subsets.keys())  # Sorted elevations

    # Determine variables to include from the first dataset
    first_ds = subsets[elevations[0]]
    variables = [var for var in first_ds.data_vars if var in ['DBZH', 'SNRH', 'DBZH_CLEAN']]
    
    # Prepare a dictionary to hold the data arrays
    data_arrays = {var: [] for var in variables}

    # Use azimuths from the first elevation
    first_azimuths = first_ds.azimuth.data
    theta = np.deg2rad(first_azimuths)

    for elevation in elevations:
        ds = subsets[elevation]
        azimuths = ds.azimuth.data
        ranges = ds.range.data
        elevation_rad = np.deg2rad(elevation)
        adjusted_ranges = ranges * np.cos(elevation_rad)

        for i_radius in range(len(ranges)):
            _r = adjusted_ranges[i_radius]
            _z = ranges[i_radius] * np.sin(elevation_rad)
            for i_theta in range(len(azimuths)):
                _theta = theta[i_theta]
                _x = _r * np.cos(_theta)
                _y = _r * np.sin(_theta)
                _id = points.InsertNextPoint(_x, _y, _z)
                
                # Check the condition to store data
                if ds['DBZH'].data[i_theta, i_radius] > 0:
                    for var in variables:
                        data_arrays[var].append(ds[var].data[i_theta, i_radius])
                else:
                    for var in variables:
                        data_arrays[var].append(np.nan)  # Use NaN for excluded data

    ugrid.SetPoints(points)

    # Add data arrays to the grid
    for var in variables:
        array_data = np.array(data_arrays[var], dtype=np.float32)  # Ensure type consistency
        vtk_data = numpy_support.numpy_to_vtk(num_array=array_data, deep=True, array_type=vtk.VTK_FLOAT)
        vtk_data.SetName(var)
        ugrid.GetPointData().AddArray(vtk_data)

    # Write the data to a file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(ugrid)
    writer.SetDataModeToBinary()  # Use binary format for compactness and speed
    writer.Write()
 
def main():
    ncinpath = '/Users/filippi_j/data/2024/ballarat/2_20240222_033000.pvol.nc'
    vtkoutpath = "/Users/filippi_j/data/2024/ballarat/output_radar_data_SimuCentered.vtp"
    
    binOutpath = "/Users/filippi_j/data/2024/ballarat/tmp/"
    nc_dir_path = '/Users/filippi_j/data/2024/ballarat/radar/AURA_2_20240222_nc 2/'
    laverton_bom_insim = (242500.0,80074.0,18.0)
    laverton_bom_timorigin = 3600*24+20.0

ncinpath = '/Users/filippi_j/data/2024/ballarat/2_20240222_033000.pvol.nc'
ds = xr.open_dataset(ncinpath)
ss = transform_dataset(ds)
# Example use case, assuming `subsets` is already loaded with your datasets
 
create_vtk_filtered_volume(ss, '3d_radar_volume.vtu')
 
#create_vtk_file(ss[elevation],elevation, "test.vtk")
#for elevation, subset in ss.items():
#   filename = f"elevation_{elevation}_data.vtk"
#   create_vtk_file(subset, elevation, filename)   
   # all_ncs = get_sorted_nc_files(nc_dir_path)
   # selected = all_ncs[:1]
    
   # print("\n".join(selected))
    #odimNC2bin(selected, binOutpath)
    
   # rawBin2vtk(binOutpath,vtkoutpath,origin=laverton_bom_insim,timeOrigin=laverton_bom_timorigin)
    #odim2MeshVTK(selected,"/Users/filippi_j/data/2024/ballarat/")

#if __name__ == "__main__":
#    main()
