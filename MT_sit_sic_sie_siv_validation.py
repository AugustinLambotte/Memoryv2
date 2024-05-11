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
lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

EGC_path_oras5 = np.where((lat >= 66) & 
                          (lat < 79) & 
                          (lon > -30) & 
                          (lon < 12)  & 
                          (lat > 60 + (80-60)/(15+22) * (lon + 22)), 1, 0)

#################################
########### - ORAS5 - ###########
#################################
sit_oras5 = []
sic_oras5 = []

surface_oras5  = np.loadtxt('ORAS5/Data/geo/surface.txt')
for file in os.listdir('ORAS5/Data/ke'):
    
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        sit_oras5.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        sic_oras5.append(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))

sit_oras5 = np.array(sit_oras5)
sic_oras5 = np.array(sic_oras5)


sit_oras5[:,EGC_path_oras5 == 0] = np.nan
sic_oras5[:,EGC_path_oras5 == 0] = np.nan

#sit_oras5[(sit_oras5 == 0) | (sic_oras5 < 0.15)] = np.nan
#sic_oras5[(sic_oras5 == 0) | (sic_oras5 < 0.15)] = np.nan

#sit_oras5[(sic_oras5 < 0.15)] = 0
#sic_oras5[(sic_oras5 < 0.15)] = 0

siv_oras5 = []
sie_oras5 = []
for i in range(len(sit_oras5)):
    siv_oras5.append(np.nansum(sit_oras5[i] * sic_oras5[i] * surface_oras5) *1e-9)
    sie_oras5.append(np.nansum(np.where(sic_oras5[i]>0.15,surface_oras5,np.nan))*1e-6)
siv_oras5 = np.array(siv_oras5)
sie_oras5 = np.array(sie_oras5)



sit_oras5 = np.nanmean(sit_oras5,axis = (1,2))
sic_oras5 = np.nanmean(sic_oras5,axis = (1,2))*1e2

sit_oras5_ra = np.convolve(sit_oras5, np.ones(12)/12, mode = 'valid')
sic_oras5_ra = np.convolve(sic_oras5, np.ones(12)/12, mode = 'valid')
siv_oras5_ra = np.convolve(siv_oras5, np.ones(12)/12, mode = 'valid')
sie_oras5_ra = np.convolve(sie_oras5, np.ones(12)/12, mode = 'valid')




###############################
########### - OBS - ###########
###############################
sit_obs = []
sic_obs = []

lon_obs = np.nan_to_num(np.loadtxt('obs_data/Data/geo/longitude.txt'))
lat_obs = np.nan_to_num(np.loadtxt('obs_data/Data/geo/latitude.txt'))

EGC_path_obs = np.where((lat_obs >= 66) & 
                          (lat_obs < 79) & 
                          (lon_obs > -30) & 
                          (lon_obs < 12)  & 
                          (lat_obs > 60 + (80-60)/(15+22) * (lon_obs + 22)), 1, 0)
surface_obs  = np.loadtxt('obs_data/Data/geo/surface.txt')
print(f'surface ORAS 5 - CS2 = {(np.nansum(np.where(EGC_path_oras5 == 1, surface_oras5,0)) - np.nansum(np.where(EGC_path_obs == 1, surface_obs,0)))*1e-6} km^2')
for file in os.listdir('obs_data/Data/gos/ke/month_mean'):
    
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2019:
        sit_obs.append(np.loadtxt(f'obs_data/Data/SI/sit/month_mean/{file}' ))
        sic_obs.append(np.loadtxt(f'obs_data/Data/SI/sic/month_mean/{file}' ))

sit_obs = np.array(sit_obs)
sic_obs = np.array(sic_obs)

sit_obs[:,EGC_path_obs == 0] = np.nan
sic_obs[:,EGC_path_obs == 0] = np.nan

#sit_obs[(sit_obs == 0)|(sic_obs <0.15)] = np.nan
#sic_obs[(sic_obs == 0)|(sic_obs <0.15)] = np.nan

#sit_obs[(sic_obs <0.15)] = 0
#sic_obs[(sic_obs <0.15)] = 0

siv_obs = []
sie_obs = []
for i in range(len(sit_obs)):
    siv_obs.append(np.nansum(sit_obs[i] * sic_obs[i] * surface_obs) *1e-9)
    sie_obs.append(np.nansum(np.where(sic_obs[i]>0.15,surface_obs,np.nan))*1e-6)
siv_obs = np.array(siv_obs)
sie_obs = np.array(sie_obs)

""" sit_obs = interpol_obs_to_oras5(sit_obs)
sit_obs = interpol_obs_to_oras5(sit_obs) """

sit_obs = np.nanmean(sit_obs,axis = (1,2))
sic_obs = np.nanmean(sic_obs,axis = (1,2))*1e2

sit_obs_ra = np.convolve(sit_obs, np.ones(12)/12, mode = 'valid')
sic_obs_ra = np.convolve(sic_obs, np.ones(12)/12, mode = 'valid')
siv_obs_ra = np.convolve(siv_obs, np.ones(12)/12, mode = 'valid')
sie_obs_ra = np.convolve(sie_obs, np.ones(12)/12, mode = 'valid')



######## - Mean error - ##########

print(f'mean error SIT on running average    : {round(np.mean(sit_oras5_ra - sit_obs_ra[:-12]),2)} m ')
print(f'mean error SIC on running average    : {round(np.mean(sic_oras5_ra - sic_obs_ra[:-12]),2)} %')
print(f'mean error SIV on running average    : {round(np.mean(siv_oras5_ra - siv_obs_ra[:-12]),2)} km^3')
print(f'mean error SIE on running average    : {round(np.mean(sie_oras5_ra - sie_obs_ra[:-12]),2)} km^2')
print('\n..............\n')
print(f'mean error SIT with month resolution : {round(np.mean(sit_oras5 - sit_obs[:-12]),2)} m ')
print(f'mean error SIC with month resolution : {round(np.mean(sic_oras5 - sic_obs[:-12]),2)} %')
print(f'mean error SIV with month resolution : {round(np.mean(siv_oras5 - siv_obs[:-12]),2)} km^3')
print(f'mean error SIE with month resolution : {round(np.mean(sie_oras5 - sie_obs[:-12]),2)} km^2')


##################################
###### - PLOT ANNUALY AV. - ######
##################################
labelsize = 25
titlesize = 25
ticklabelsize = 21
legendsize = 15    

fig,axes = plt.subplots(nrows = 2, ncols = 2,figsize = (18,10),sharex = True)

axes[0,0].plot(np.linspace(2011,2019,len(sit_oras5))   ,sit_oras5   ,color = 'red',alpha = 0.4)
axes[0,0].plot(np.linspace(2012,2019,len(sit_oras5_ra)),sit_oras5_ra,color = 'red',linewidth = 3,label = 'ORAS5')
axes[0,0].set_title(f'SIT',fontsize = titlesize)
axes[0,0].grid()
axes[0,0].tick_params(labelsize = ticklabelsize)
axes[0,0].set_ylabel('m',size = labelsize)
axes[0,0].set_xlim(2011,2020)

axes[0,0].plot(np.linspace(2011,2020,len(sit_obs))   ,sit_obs   ,color = 'blue',alpha = 0.4)
axes[0,0].plot(np.linspace(2012,2020,len(sit_obs_ra)),sit_obs_ra,color = 'blue',linewidth = 3,label = 'CryoSat-2')
axes[0,0].legend(fontsize = legendsize)

################
axes[0,1].plot(np.linspace(2011,2019,len(sic_oras5))   ,sic_oras5   ,color = 'red',alpha = 0.4)
axes[0,1].plot(np.linspace(2012,2019,len(sic_oras5_ra)),sic_oras5_ra,color = 'red',linewidth = 3,label = 'ORAS5')
axes[0,1].set_title(f'SIC',fontsize = titlesize)
axes[0,1].grid()
axes[0,1].tick_params(labelsize = ticklabelsize)
axes[0,1].set_ylabel('%',size = labelsize)
axes[0,1].set_xlim(2011,2020)

axes[0,1].plot(np.linspace(2011,2020,len(sic_obs))   ,sic_obs   ,color = 'blue',alpha = 0.4)
axes[0,1].plot(np.linspace(2012,2020,len(sic_obs_ra)),sic_obs_ra,color = 'blue',linewidth = 3,label = 'OSI-SAF')
axes[0,1].legend(fontsize = legendsize)
#################
axes[1,0].plot(np.linspace(2011,2019,len(siv_oras5))   ,siv_oras5   ,color = 'red',alpha = 0.4)
axes[1,0].plot(np.linspace(2012,2019,len(siv_oras5_ra)),siv_oras5_ra,color = 'red',linewidth = 3,label = 'ORAS5')
axes[1,0].set_title(f'SIV',fontsize = titlesize)
axes[1,0].grid()
axes[1,0].tick_params(labelsize = ticklabelsize)
axes[1,0].set_ylabel(r'$km^3$',size = labelsize)
axes[1,0].set_xlim(2011,2020)

axes[1,0].plot(np.linspace(2011,2020,len(siv_obs))   ,siv_obs   ,color = 'blue',alpha = 0.4)
axes[1,0].plot(np.linspace(2012,2020,len(siv_obs_ra)),siv_obs_ra,color = 'blue',linewidth = 3,label = 'CryoSat-2')
axes[1,0].legend(fontsize = legendsize)
################
axes[1,1].plot(np.linspace(2011,2019,len(sie_oras5))   ,sie_oras5*1e-5   ,color = 'red',alpha = 0.4)
axes[1,1].plot(np.linspace(2012,2019,len(sie_oras5_ra)),sie_oras5_ra*1e-5,color = 'red',linewidth = 3,label = 'ORAS5')
axes[1,1].set_title(f'SIE',fontsize = titlesize)
axes[1,1].grid()
axes[1,1].tick_params(labelsize = ticklabelsize)
axes[1,1].set_ylabel(r'$10^5km^2$',size = labelsize)
axes[1,1].set_xlim(2011,2020)

axes[1,1].plot(np.linspace(2011,2020,len(sie_obs))   ,sie_obs*1e-5   ,color = 'blue',alpha = 0.4)
axes[1,1].plot(np.linspace(2012,2020,len(sie_obs_ra)),sie_obs_ra*1e-5,color = 'blue',linewidth = 3,label = 'OSI-SAF')
axes[1,1].legend(fontsize = legendsize)
################


plt.savefig('validation/MT_sit_sic_siv_sie_averaged.png')
plt.close()

##################################
######## - Monthly plot - ########
##################################
sie_obs_monthly = sie_obs[:-12].reshape((8,12)).T
siv_obs_monthly = siv_obs[:-12].reshape((8,12)).T
sit_obs_monthly = sit_obs[:-12].reshape((8,12)).T
sic_obs_monthly = sic_obs[:-12].reshape((8,12)).T

sit_oras5_monthly = sit_oras5.reshape((8,12)).T
sic_oras5_monthly = sic_oras5.reshape((8,12)).T
sie_oras5_monthly = sie_oras5.reshape((8,12)).T
siv_oras5_monthly = siv_oras5.reshape((8,12)).T

sie_obs_monthly_std = np.std(sie_obs_monthly,axis = 1)
siv_obs_monthly_std = np.std(siv_obs_monthly,axis = 1)
sit_obs_monthly_std = np.std(sit_obs_monthly,axis = 1)
sic_obs_monthly_std = np.std(sic_obs_monthly,axis = 1)

sit_oras5_monthly_std = np.std(sit_oras5_monthly,axis = 1)
sic_oras5_monthly_std = np.std(sic_oras5_monthly,axis = 1)
sie_oras5_monthly_std = np.std(sie_oras5_monthly,axis = 1)
siv_oras5_monthly_std = np.std(siv_oras5_monthly,axis = 1)

sie_obs_monthly_mean = np.mean(sie_obs_monthly,axis = 1)
siv_obs_monthly_mean = np.mean(siv_obs_monthly,axis = 1)
sit_obs_monthly_mean = np.mean(sit_obs_monthly,axis = 1)
sic_obs_monthly_mean = np.mean(sic_obs_monthly,axis = 1)

sit_oras5_monthly_mean = np.mean(sit_oras5_monthly,axis = 1)
sic_oras5_monthly_mean = np.mean(sic_oras5_monthly,axis = 1)
sie_oras5_monthly_mean = np.mean(sie_oras5_monthly,axis = 1)
siv_oras5_monthly_mean = np.mean(siv_oras5_monthly,axis = 1)

fig,axes = plt.subplots(nrows = 2, ncols = 2,figsize = (18,10),sharex = True)
plt.subplots_adjust(right=0.87)

axes[0,0].errorbar(month_name,sit_obs_monthly_mean,sit_obs_monthly_std,color = 'blue',label = 'CryoSat-2',marker = 's',capsize=5,linestyle = 'dashed')
axes[0,0].errorbar(month_name,sit_oras5_monthly_mean,sit_oras5_monthly_std,color = 'red',label = 'ORAS5',marker = 's',capsize=5,linestyle = 'dashed')

axes[0,0].set_title(f'SIT',fontsize = titlesize)
axes[0,0].grid()
axes[0,0].legend(fontsize = legendsize)
axes[0,0].tick_params(labelsize = ticklabelsize)
axes[0,0].set_ylabel('m',size = labelsize)

axes[0,1].errorbar(month_name,sic_obs_monthly_mean,sic_obs_monthly_std,color = 'blue',label = 'OSI-SAF',marker = 's',capsize=5,linestyle = 'dashed')
axes[0,1].errorbar(month_name,sic_oras5_monthly_mean,sic_oras5_monthly_std,color = 'red',label = 'ORAS5',marker = 's',capsize=5,linestyle = 'dashed')

axes[0,1].set_title(f'SIC',fontsize = titlesize)
axes[0,1].legend(fontsize = legendsize)
axes[0,1].grid()
axes[0,1].tick_params(labelsize = ticklabelsize)
axes[0,1].set_ylabel('%',size = labelsize)

axes[1,1].errorbar(month_name,sie_obs_monthly_mean*1e-5,sie_obs_monthly_std*1e-5,color = 'blue',label = 'OSI-SAF',marker = 's',capsize=5,linestyle = 'dashed')
axes[1,1].errorbar(month_name,sie_oras5_monthly_mean*1e-5,sie_oras5_monthly_std*1e-5,color = 'red',label = 'ORAS5',marker = 's',capsize=5,linestyle = 'dashed')
axes[1,1].set_title(f'SIE',fontsize = titlesize)
axes[1,1].grid()
axes[1,1].legend(fontsize = legendsize)
axes[1,1].tick_params(labelsize = ticklabelsize)
axes[1,1].set_ylabel(r'$10^5km^2$',size = labelsize)

axes[1,0].errorbar(month_name,siv_obs_monthly_mean,siv_obs_monthly_std,color = 'blue',label = 'CryoSat-2',marker = 's',capsize=5,linestyle = 'dashed')
axes[1,0].errorbar(month_name,siv_oras5_monthly_mean,siv_oras5_monthly_std,color = 'red',label = 'ORAS5',marker = 's',capsize=5,linestyle = 'dashed')
axes[1,0].set_title(f'SIV',fontsize = titlesize)
axes[1,0].legend(fontsize = legendsize)
axes[1,0].grid()
axes[1,0].tick_params(labelsize = ticklabelsize)
axes[1,0].set_ylabel(r'$km^3$',size = labelsize)

ax2 = fig.add_axes([0.75,0.4,0.3,0.3],projection = ccrs.LambertConformal(central_longitude = -6))
xlim = [-33, 20]
ylim = [65, 80]
lower_space = 3 
rect = mpath.Path([[xlim[0], ylim[0]],
                [xlim[1], ylim[0]],
                [xlim[1], ylim[1]],
                [xlim[0], ylim[1]],
                [xlim[0], ylim[0]],
                ]).interpolated(20)
proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax2) - ax2.transData
rect_in_target = proj_to_data.transform_path(rect)
ax2.set_boundary(rect_in_target)
ax2.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

ax2.coastlines()
ax2.gridlines()        

ds_bathy = xr.open_dataset('ORAS5/Data/geo/ETOPO_2022_v1_60s_N90W180_surface.nc')
lat_bathy = np.array(ds_bathy['lat'].where((ds_bathy.lat < 81) & (ds_bathy.lat > 60) & (ds_bathy.lon> -40) & (ds_bathy.lon < 20)))
lon_bathy = np.array(ds_bathy['lon'].where((ds_bathy.lat < 81) & (ds_bathy.lat > 60) & (ds_bathy.lon> -40) & (ds_bathy.lon < 20))).T
z         = np.array(ds_bathy['z'].where(  (ds_bathy.lat < 81) & (ds_bathy.lat > 60) & (ds_bathy.lon> -40) & (ds_bathy.lon < 20)))
z         = z[        9000:10200,8430:11990]
lat_bathy = lat_bathy[9000:10200,8430:11990]
lon_bathy = lon_bathy[9000:10200,8430:11990]

#Reducing the size of the array
for i in range(0):
    lat_bathy = np.delete(lat_bathy,[i for i in range(0,np.shape(lat_bathy)[0], 2)], axis = 0)
    lat_bathy = np.delete(lat_bathy,[i for i in range(0,np.shape(lat_bathy)[1], 2)], axis = 1)

    lon_bathy = np.delete(lon_bathy,[i for i in range(0,np.shape(lon_bathy)[0], 2)], axis = 0)
    lon_bathy = np.delete(lon_bathy,[i for i in range(0,np.shape(lon_bathy)[1], 2)], axis = 1)

    z = np.delete(z,[i for i in range(0,np.shape(z)[0], 2)], axis = 0)
    z = np.delete(z,[i for i in range(0,np.shape(z)[1], 2)], axis = 1)

lat_bathy = np.nan_to_num(lat_bathy)
lon_bathy = np.nan_to_num(lon_bathy)
z = np.nan_to_num(z)
z[z>0] = 0
z *= -1
cs = ax2.pcolormesh(lon_bathy, lat_bathy, z,cmap = "cmo.dense", transform=ccrs.PlateCarree())

cs = ax2.contour(lon, lat, EGC_path_oras5, levels = [0,1],linestyle = 'bold',linewidths = 3, colors = ['red'],transform=ccrs.PlateCarree())

plt.savefig('MT_plot/validation/MT_sit_sic_siv_sie_month.png')
