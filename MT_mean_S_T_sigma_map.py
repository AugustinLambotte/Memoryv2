import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import xarray as xr
from scipy import stats
import numpy.ma as ma
from datetime import date
import seawater as sw
from ORAS5.useful_func_oras5 import plot
"""
    This script is used to plot on a map the value of the trend over all the time span for the EGC energy and the fresh water flux and for each pixel
"""
earth_radius = 6370*1e3
lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

"""
    return maps of the trend of sim and EGC energy
"""    
EGC_path = np.loadtxt('ORAS5/Data/geo/Current_path.txt')

S0               = []
T0               = []
sigma0               = []
S30               = []
T30               = []
sigma30               = []
S60               = []
T60               = []
sigma60               = []

day_btwn_2011_start_and_2018_end = []
for file in os.listdir('ORAS5/Data/ke'):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        
        
        Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
        Tds = xr.open_dataset(f'ORAS5/Data/votemper/{file[:-3]}nc')

        S0.append(Sds['vosaline'].sel(deptht = 0, method = 'nearest')[:])
        T0.append(Tds['votemper'].sel(deptht = 0, method = 'nearest')[:])
        sigma0.append(sw.pden(S0[-1],T0[-1],sw.pres(   0,lat))-1000)
        S30.append(Sds['vosaline'].sel(deptht = 30, method = 'nearest')[:])
        T30.append(Tds['votemper'].sel(deptht = 30, method = 'nearest')[:])
        sigma30.append(sw.pden(S30[-1],T30[-1],sw.pres( 30,lat))-1000)
        S60.append(Sds['vosaline'].sel(deptht = 60, method = 'nearest')[:])
        T60.append(Tds['votemper'].sel(deptht = 60, method = 'nearest')[:])
        sigma60.append(sw.pden(S60[-1],T60[-1],sw.pres( 60,lat))-1000)
        Sds.close()
        Tds.close()

S0               = np.nan_to_num(np.array(S0))
T0               = np.nan_to_num(np.array(T0))
sigma0           = np.nan_to_num(np.array(sigma0))
S30              = np.nan_to_num(np.array(S30))
T30              = np.nan_to_num(np.array(T30))
sigma30          = np.nan_to_num(np.array(sigma30))
S60              = np.nan_to_num(np.array(S60))
T60              = np.nan_to_num(np.array(T60))
sigma60          = np.nan_to_num(np.array(sigma60))

S0               = np.mean(np.array(S0)     , axis = 0)
T0               = np.mean(np.array(T0)     , axis = 0)
sigma0           = np.mean(np.array(sigma0) , axis = 0)
S30              = np.mean(np.array(S30)    , axis = 0)
T30              = np.mean(np.array(T30)    , axis = 0)
sigma30          = np.mean(np.array(sigma30), axis = 0)
S60              = np.mean(np.array(S60)    , axis = 0)
T60              = np.mean(np.array(T60)    , axis = 0)
sigma60          = np.mean(np.array(sigma60), axis = 0)


#####################################
#### - Plotting S, T and sigma - ####
#####################################

titles_list = ['S\nSurface', 'T\nSurface',r'$\sigma_{\theta}$'+'\nSurface','30m', '30m','30m','60m','60m','60m']

label_list = [r'$[psu]$',r'$[°C]$',r'$[kg*m^{-3}]$',r'$[psu]$',r'$[°C]$',r'$[kg*m^{-3}]$',r'$[psu]$',r'$[°C]$',r'$[kg*m^{-3}]$']
v_list = [[32,35],[-2.5,4],[26,28],[32,35],[-2.5,4],[26,28],[32,35],[-2.5,4],[26,28]]
cmap = ['cmo.haline','cmo.thermal','cmo.dense','cmo.haline','cmo.thermal','cmo.dense','cmo.haline','cmo.thermal','cmo.dense']
titlesize = 25
labelsize = 17
fig,axes = plt.subplots(ncols = 3, nrows = 3, figsize = (25,18), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -6)})
plt.subplots_adjust(left=0.01,
                    bottom=0.05, 
                    right=0.9, 
                    top=0.93,)


for trend_map,i in zip([S0,T0,sigma0,S30,T30,sigma30,S60,T60,sigma60],range(9)):
    
    xlim = [-33, 20]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axes.flatten()[i]) - axes.flatten()[i].transData
    rect_in_target = proj_to_data.transform_path(rect)
    axes.flatten()[i].set_boundary(rect_in_target)
    axes.flatten()[i].set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

    axes.flatten()[i].coastlines()
    axes.flatten()[i].gridlines()        

    cs = axes.flatten()[i].pcolormesh(lon, lat, trend_map,vmin = v_list[i][0], vmax = v_list[i][1],cmap = cmap[i], transform=ccrs.PlateCarree())
    contour = axes.flatten()[i].contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.01,axes.flatten()[i].get_position().y0 - 0.02,0.015,axes.flatten()[i].get_position().height])
    cb = plt.colorbar(cs, cax = cax)
    cb.ax.tick_params(labelsize=labelsize)
    cb.set_label(label = label_list[i],size = labelsize)

    
    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    

plt.savefig(f'MT_plot/mean_S_T_sigma.png')
plt.close()