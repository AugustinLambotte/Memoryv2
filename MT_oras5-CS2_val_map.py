import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import xarray as xr
from scipy import stats
import numpy.ma as ma
from ORAS5.useful_func_oras5 import plot as plot_oras
from ORAS5.useful_func_oras5 import show as show_ORAS
from ORAS5.useful_func_oras5 import interpol_obs_to_oras5
from obs_data.useful_func import plot as plot_obs
from obs_data.useful_func import show as show_obs

from datetime import date


month_name = ['Ja','Fe','Mar','Ap','May','Ju','July','Au','Se','Oc','No','De']

#################################
########### - ORAS5 - ###########
#################################
sit_oras5 = []
sic_oras5 = []

EGC_path_oras5 = np.nan_to_num(np.loadtxt('ORAS5/Data/geo/EGC_path.txt'))
surface_oras5  = np.loadtxt('ORAS5/Data/geo/surface.txt')
for file in os.listdir('ORAS5/Data/ke'):
    
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        sit_oras5.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        sic_oras5.append(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))

sit_oras5 = np.array(sit_oras5)
sic_oras5 = np.array(sic_oras5)
sit_oras5[:,EGC_path_oras5 == 0] = np.nan
sic_oras5[:,EGC_path_oras5 == 0] = np.nan

###############################
########### - OBS - ###########
###############################
sit_obs = []
sic_obs = []

EGC_path_obs = np.nan_to_num(np.loadtxt('obs_data/Data/geo/EGC_path.txt'))
surface_obs  = np.loadtxt('obs_data/Data/geo/surface.txt')
print(f'surface ORAS 5 - CS2 = {(np.nansum(np.where(EGC_path_oras5 == 1, surface_oras5,0)) - np.nansum(np.where(EGC_path_obs == 1, surface_obs,0)))*1e-6} km^2')
for file in os.listdir('obs_data/Data/gos/ke/month_mean'):
    
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        sit_obs.append(interpol_obs_to_oras5(np.loadtxt(f'obs_data/Data/SI/sit/month_mean/{file}' )))
        sic_obs.append(interpol_obs_to_oras5(np.loadtxt(f'obs_data/Data/SI/sic/month_mean/{file}' )))

sit_obs = np.array(sit_obs)
sic_obs = np.array(sic_obs)

sit_obs[:,EGC_path_oras5 == 0] = np.nan
sic_obs[:,EGC_path_oras5 == 0] = np.nan

print(f'mean error SIT with month resolution : {round(np.mean(sit_oras5 - sit_obs),2)} m ')
print(f'mean error SIC with month resolution : {round(np.mean(sic_oras5 - sic_obs),2)} %')


##################################
###### - PLOT ANNUALY AV. - ######
##################################

error_sit = sit_oras5-sit_obs
error_sic = sic_oras5-sic_obs

plot_oras(np.nanmean(error_sit,axis = 0))
plot_oras(np.nanmean(error_sic,axis = 0))
