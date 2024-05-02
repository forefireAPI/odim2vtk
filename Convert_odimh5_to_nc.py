import pyart
import os
import glob

## Level 1 data in h5
path_lvl1 = r'/Users/uqaguyot_local/Documents/python-master/2-Current/GoogleProject/data/AURA/2/Level_1/2_20240222'
os.chdir(path_lvl1)
file_list_l1 = sorted(glob.glob('*.h5*'))

## Export folder to netCDF
path_export = r'/Users/uqaguyot_local/Documents/python-master/2-Current/GoogleProject/data/AURA/2/export_to_netCDF'
os.chdir(path_export)

for i in np.arange(0,len(file_list_l1),1):
    
    os.chdir(path_lvl1)
    my_radar_l1 = pyart.aux_io.read_odim_h5(file_list_l1[i], file_field_names=True)
    
    os.chdir(path_export)
    pyart.io.write_cfradial(file_list_l1[i][0:22]+".nc", my_radar_l1, format='NETCDF4', time_reference=None, arm_time_variables=False)