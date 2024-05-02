
# Timed VTK Converter for Radar Data

This project provides a Python script that converts radar data from ODIM HDF5 or NetCDF format into VTK format, incorporating timing information. This conversion facilitates advanced visualization and analysis in software that supports VTK, such as ParaView or VTK-based custom tools.

## Features

- Converts ODIM HDF5 or NetCDF radar data into VTK format.
- Includes time difference calculation to track the changes over different scans.
- Extracts and processes radar parameters like reflectivity (DBZH), differential reflectivity (ZDR), and copolar correlation coefficient (RHOHV).

## Requirements

- Python 3.x
- NumPy
- xarray
- VTK
- HDF5 or NetCDF4 libraries installed (depending on your data format)


## Usage

1. Place your ODIM HDF5 or NetCDF radar data file in a known directory.
2. Modify the `file_path` in the script to point to your input data file.
3. Run the script to perform the conversion:

```bash
python radar_to_vtk.py
```

4. The output will be a VTK file with timing information, which can be visualized in ParaView or any VTK-compatible software.

## Customization

- You can adjust the radar parameters extracted by modifying the `main()` function in the script.
- To handle large datasets efficiently, consider implementing parallel processing or optimizing memory usage.

## Output

The script will generate a `.vtp` file containing the radar data with associated scalar values like time differences and radar measurements, structured for easy 3D visualization.

## Authors

- Name: A.Guyot, J.B. Filippi

## License

This project is licensed under the MIT License - see the LICENSE file for details.
