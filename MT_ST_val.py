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
from scipy import interpolate


earth_radius = 6370*1e3
lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

"""
    return maps of the trend of sim and EGC energy
"""    
EGC_path = np.loadtxt('ORAS5/Data/geo/Current_path.txt')

##########################
##### - Copernicus - #####
##########################

ds = xr.open_dataset('obs_data/Data/Sal_dens/cmems_obs-mob_glo_phy-sss_my_multi_P1M_1712234106978.nc')
print(ds)
S_obsMarch     = np.array(ds['sos'].where((ds['sos'].time.dt.month == 3),drop = True).sel(depth = 0))
sigma_obsMarch = np.array(ds['dos'].where((ds['dos'].time.dt.month == 3),drop = True).sel(depth = 0))-1000

S_obsSept     = np.array(ds['sos'].where((ds['sos'].time.dt.month == 9),drop = True).sel(depth = 0))
sigma_obsSept = np.array(ds['dos'].where((ds['dos'].time.dt.month == 9),drop = True).sel(depth = 0))-1000
latObs = np.array(ds['latitude'])
lonObs = np.array(ds['longitude'])

lonObs,latObs = np.meshgrid(lonObs,latObs)
ds.close()


S_obsMarch     = np.mean(np.array(S_obsMarch)     , axis = 0)
sigma_obsMarch = np.mean(np.array(sigma_obsMarch) , axis = 0)

S_obsSept      = np.mean(np.array(S_obsSept)      , axis = 0)
sigma_obsSept  = np.mean(np.array(sigma_obsSept)  , axis = 0)



#####################
##### - ORAS5 - #####
#####################

S_oras5March     = []
sigma_oras5March = []

S_oras5Sept      = []
sigma_oras5Sept  = []

for file in os.listdir('ORAS5/Data/ke'):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        
        if int(file[5:7]) == 3:
            Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
            Tds = xr.open_dataset(f'ORAS5/Data/votemper/{file[:-3]}nc')

            S_oras5March.append(Sds['vosaline'].sel(deptht = 0, method = 'nearest')[:])
            T = Tds['votemper'].sel(deptht = 0, method = 'nearest')[:]
            sigma_oras5March.append(sw.dens0(S_oras5March[-1],T)-1000)
            Sds.close()
            Tds.close()

        if int(file[5:7]) == 9:
            Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
            Tds = xr.open_dataset(f'ORAS5/Data/votemper/{file[:-3]}nc')

            S_oras5Sept.append(Sds['vosaline'].sel(deptht = 0, method = 'nearest')[:])
            T = Tds['votemper'].sel(deptht = 0, method = 'nearest')[:]
            sigma_oras5Sept.append(sw.dens0(S_oras5Sept[-1],T)-1000)
            Sds.close()
            Tds.close()


S_oras5March     = np.mean(np.array(S_oras5March)     , axis = 0)
sigma_oras5March = np.mean(np.array(sigma_oras5March) , axis = 0)

S_oras5Sept      = np.mean(np.array(S_oras5Sept)      , axis = 0)
sigma_oras5Sept  = np.mean(np.array(sigma_oras5Sept)  , axis = 0)

############################################
#####- Interpolation over ORAS5 grid - #####
############################################


points = [] # list of length NxM containing all the coordinates [lat,lon] of all points from si drift map
values_Smarch = []
values_Ssept = []
values_sigmmarch = []
values_sigmsept = []

#data[data == 0] = np.nan
for i in range(len(latObs)):
    for j in range(len(lonObs[0])):
        #if (data[i,j] !=0)  and (not np.isnan(data[i,j])):
        points.append([latObs[i,j],lonObs[i,j]])
        values_Smarch.append(S_obsMarch[i,j])
        values_Ssept.append(S_obsSept[i,j])
        values_sigmmarch.append(sigma_obsMarch[i,j])
        values_sigmsept.append(sigma_obsSept[i,j])
points = np.array(points)
values_Smarch    = np.array(values_Smarch)
values_Ssept     = np.array(values_Ssept)
values_sigmmarch = np.array(values_sigmmarch)
values_sigmsept  = np.array(values_sigmsept)

S_marchInterp    = interpolate.griddata(points, values_Smarch,    (lat, lon), method='nearest')
S_septInterp     = interpolate.griddata(points, values_Ssept,     (lat, lon), method='nearest')
sigm_marchInterp = interpolate.griddata(points, values_sigmmarch, (lat, lon), method='nearest')
sigm_septInterp  = interpolate.griddata(points, values_sigmsept,  (lat, lon), method='nearest')

diff_S_march    = S_oras5March     - S_marchInterp
diff_S_sept     = S_oras5Sept      - S_septInterp
diff_sigm_march = sigma_oras5March - sigm_marchInterp
diff_sigm_sept  = sigma_oras5Sept  - sigm_septInterp

##################################
#### - Plotting S and sigma - ####
##################################


titles_list = ['March\nCopernicus', 'September\nCopernicus','ORAS5','ORAS5', 'Difference','Difference']

min = 31
v_list = [[min,35],[min,35],[min,35],[min,35],[-1.5,1.5],[-1.5,1.5]]
cmap = ['cmo.haline','cmo.haline','cmo.haline','cmo.haline','cmo.balance','cmo.balance',]
titlesize = 25
labelsize = 17
fig,axes = plt.subplots(ncols = 2, nrows = 3, figsize = (15,18), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -8)})
plt.subplots_adjust(left=0.01,
                    bottom=0.05, 
                    right=0.9, 
                    top=0.93,
                    wspace = 0)


for map,i in zip([S_obsMarch,S_obsSept,S_oras5March,S_oras5Sept, diff_S_march, diff_S_sept],range(6)):
    
    xlim = [-33, 16]
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
    if i < 2:
        cs = axes.flatten()[i].pcolormesh(lonObs, latObs, map,vmin = v_list[i][0], vmax = v_list[i][1],cmap = cmap[i], transform=ccrs.PlateCarree())
    else:
        cs = axes.flatten()[i].pcolormesh(lon, lat, map,vmin = v_list[i][0], vmax = v_list[i][1],cmap = cmap[i], transform=ccrs.PlateCarree())
    
    if i == 3:
        cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.01,axes.flatten()[i].get_position().y0+0.03,0.03,axes.flatten()[i].get_position().height*2])
        cb = plt.colorbar(cs, cax = cax)
        cb.ax.tick_params(labelsize=labelsize)
        cb.set_label(label = '[psu]',size = labelsize)

    
    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)

cax = fig.add_axes([axes.flatten()[-1].get_position().x1+0.01,axes.flatten()[-1].get_position().y0,0.03,axes.flatten()[-1].get_position().height])
cb = plt.colorbar(cs, cax = cax)
cb.ax.tick_params(labelsize=labelsize)
cb.set_label(label = '[psu]',size = labelsize)
plt.savefig(f'MT_plot/validation/salinity.png')
plt.close()

###########################
####### - Density - #######
###########################

titles_list = ['March\nCopernicus', 'September\nCopernicus','ORAS5','ORAS5', 'Difference','Difference']

min = 23
max = 29
balance = 1.5
v_list = [[min,max],[min,max],[min,max],[min,max],[-balance,balance],[-balance,balance]]
cmap = ['cmo.dense','cmo.dense','cmo.dense','cmo.dense','cmo.balance','cmo.balance']
titlesize = 25
labelsize = 17

fig,axes = plt.subplots(ncols = 2, nrows = 3, figsize = (15,18), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -8)})
plt.subplots_adjust(left=0.01,
                    bottom=0.05, 
                    right=0.9, 
                    top=0.93,
                    wspace = 0)

for map,i in zip([sigma_obsMarch,sigma_obsSept,sigma_oras5March,sigma_oras5Sept, diff_sigm_march, diff_sigm_sept],range(6)):
    
    xlim = [-33, 16]
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
    if i < 2:
        cs = axes.flatten()[i].pcolormesh(lonObs, latObs, map,vmin = v_list[i][0], vmax = v_list[i][1],cmap = cmap[i], transform=ccrs.PlateCarree())
    else:
        cs = axes.flatten()[i].pcolormesh(lon, lat, map,vmin = v_list[i][0], vmax = v_list[i][1],cmap = cmap[i], transform=ccrs.PlateCarree())
    
    if i == 3:
        cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.01,axes.flatten()[i].get_position().y0+0.03,0.03,axes.flatten()[i].get_position().height*2])
        cb = plt.colorbar(cs, cax = cax)
        cb.ax.tick_params(labelsize=labelsize)
        cb.set_label(label = r'$[kg*m^{-3}]$',size = labelsize)

    
    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)

cax = fig.add_axes([axes.flatten()[-1].get_position().x1+0.01,axes.flatten()[-1].get_position().y0,0.03,axes.flatten()[-1].get_position().height])
cb = plt.colorbar(cs, cax = cax)
cb.ax.tick_params(labelsize=labelsize)
cb.set_label(label = r'$[kg*m^{-3}]$',size = labelsize)
plt.savefig(f'MT_plot/validation/density.png')
plt.close()