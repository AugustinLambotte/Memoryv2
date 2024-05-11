import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os 
from datetime import timedelta, date
from scipy import interpolate, stats
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import math
import seawater as sw
from ORAS5.useful_func_oras5 import plot, show
from dateutil.relativedelta import *

lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

surface   = np.loadtxt('ORAS5/Data/geo/Surface.txt')
X_dist    = np.loadtxt('ORAS5/Data/geo/X_dist.txt')
Y_dist    = np.loadtxt('ORAS5/Data/geo/Y_dist.txt')

def FS_MB(file):
    iicevelv = np.loadtxt(f'ORAS5/Data/iicevelv/{file}')
    iicevelu = np.loadtxt(f'ORAS5/Data/iicevelu/{file}')
    sit      = np.loadtxt(f'ORAS5/Data/iicethic/{file}')
    sic      = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')
    """ plot(Y_dist)
    plot(np.where(((lat == 82.5) & (lon<= 17) & (lon >= -12)) | ((lon == 18) & (lat <= 81.5) & (lat >= 80.5)),
                                    1, np.nan)) """
    
    fromNorth = np.nansum(np.where((lat == 82.5) & (lon<= 17) & (lon >= -14),
                                    sit*sic*X_dist*iicevelv, 0)) #m^3/s
    fromEast  = np.nansum(np.where((lon == 18) & (lat <= 81.5) & (lat >= 80.5),
                                    sit*sic*Y_dist*iicevelu, 0)) #m^3/s
    MB = (fromNorth + fromEast)*60*60*24*30 #m^3/month
    return MB


###### - MB - ######

day_btwn_2012_start_and_2018_end = []
MB_list = []
for file in os.listdir('ORAS5/Data/iicevelu'): 
    file_year  = int(file[:4])
    file_month = int(file[5:7])
    if file_year < 2011:
        continue
    print(f'>>>   Computing FS_MB for month {file_year} - {file_month}   <<<')
    if file_year >= 2012 or file_month == 12:
        day_btwn_2012_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)
    MB_list.append(FS_MB(file))

MB_list = -1*np.array(MB_list)*1e-9 # m^3 -> km^3
M1 = -1*np.loadtxt('ORAS5/Data/M1.txt').flatten()
M2 = -1*np.loadtxt('ORAS5/Data/M2.txt').flatten()
R  = -1*np.loadtxt('ORAS5/Data/Ricker.txt').flatten()

print(np.shape(M1))

fig,axes = plt.subplots(nrows = 2,figsize = (10,9))
labelsize = 20
ticklabelsize = 20
legendsize = 13

plt.subplots_adjust(
                    bottom=0.05,  
                    top=0.93,)
#####################
##### - TOTAL - #####
#####################
axes[0].plot(np.linspace(2011,2017,len(MB_list[:-24])+1),np.append(MB_list[:-24],np.nan),label = 'ORAS5', marker = 'x' )
axes[0].plot(np.linspace(2011,2017,len(M1)+1),np.append(M1,np.nan),label = 'M1', marker = 'o',alpha = 0.8 )
axes[0].plot(np.linspace(2011,2017,len(M1)+1),np.append(M2,np.nan),label = 'M2', marker = 's',alpha = 0.8)
axes[0].plot(np.linspace(2011,2017,len(M1)+1),np.append(R,np.nan),label = 'R', marker = 'v',alpha = 0.8)
axes[0].set_ylabel(r'$km^3\cdot month^{-1}$',size = labelsize)
axes[0].tick_params(labelsize = ticklabelsize)
axes[0].legend(fontsize = legendsize)
axes[0].grid()

##############################
##### - Seasonal Cycle - #####
##############################

axes[1].errorbar(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],np.mean(np.reshape(MB_list[:-24],(6,12)),axis = 0),np.std(np.reshape(MB_list[:-24],(6,12)),axis = 0),label = 'ORAS5', marker = 's',capsize=5,linestyle = 'dashed')
axes[1].errorbar(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],np.mean(np.reshape(M1,           (6,12)),axis = 0),np.std(np.reshape(M1,           (6,12)),axis = 0),label = 'M1',    marker = 's',capsize=5,linestyle = 'dashed')
axes[1].errorbar(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],np.mean(np.reshape(M2,           (6,12)),axis = 0),np.std(np.reshape(M2,           (6,12)),axis = 0),label = 'M2',    marker = 's',capsize=5,linestyle = 'dashed')
axes[1].errorbar(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],np.mean(np.reshape(R,            (6,12)),axis = 0),np.std(np.reshape(R,            (6,12)),axis = 0),label = 'R',     marker = 's',capsize=5,linestyle = 'dashed')
axes[1].legend(fontsize = legendsize)
axes[1].set_ylabel(r'$km^3\cdot month^{-1}$',size = labelsize)
axes[1].tick_params(labelsize = ticklabelsize)
axes[1].grid()


ax2 = fig.add_axes([0.455,0.22,0.25,0.25],projection = ccrs.LambertConformal(central_longitude = -6))
xlim = [-28, 23]
ylim = [75, 84]
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
lat_bathy = np.array(ds_bathy['lat'].where((ds_bathy.lat > 70) & (ds_bathy.lon> -30) & (ds_bathy.lon < 35)))
lon_bathy = np.array(ds_bathy['lon'].where((ds_bathy.lat > 70) & (ds_bathy.lon> -30) & (ds_bathy.lon < 35))).T
z         = np.array(ds_bathy['z'].where(  (ds_bathy.lat > 70) & (ds_bathy.lon> -30) & (ds_bathy.lon < 35)))
z         = z[        9000:10500,8430:12500]
lat_bathy = lat_bathy[9000:10500,8430:12500]

lon_bathy = lon_bathy[9000:10500,8430:12500]

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

gate = np.concatenate(([[82,i] for i in range(-12,18)],[[81,17],[80,17]]))
ax2.plot(np.array(gate)[:,1], np.array(gate)[:,0], linewidth = 4,color = 'red',transform=ccrs.PlateCarree())

#ax2.text(-16,81.7,'G',transform=ccrs.PlateCarree(),fontsize=20,weight='bold')

plt.savefig('MT_plot/validation/FSMB.png')
