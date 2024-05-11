import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import os
from ORAS5.useful_func_oras5 import plot, vector_plot, interpol_obs_to_oras5
import seawater as sw

lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")


sit_summer   = []
sal_summer   = []
temp_summer  = []
dens_summer  = []

sit_winter   = []
sal_winter   = []
temp_winter  = []
dens_winter  = []

for file in os.listdir('ORAS5/Data/iicethic'):
    if int(file[:4]) < 2011 or int(file[:4]) > 2018:
        continue
    if int(file[5:7]) in [3,4,5]:
        sit_winter.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
        Tds = xr.open_dataset(f'ORAS5/Data/votemper/{file[:-3]}nc')

        sal_winter.append( Sds['vosaline'].sel(deptht = 0, method = 'nearest')[:])
        temp_winter.append(Tds['votemper'].sel(deptht = 0, method = 'nearest')[:])
        dens_winter.append(sw.pden(sal_winter[-1],temp_winter[-1],sw.pres(0,lat))-1000)
        Sds.close()
        Tds.close()

    if int(file[5:7]) in [7,8,9]:
        sit_summer.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))

        Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
        Tds = xr.open_dataset(f'ORAS5/Data/votemper/{file[:-3]}nc')

        sal_summer.append( Sds['vosaline'].sel(deptht = 0, method = 'nearest')[:])
        temp_summer.append(Tds['votemper'].sel(deptht = 0, method = 'nearest')[:])
        dens_summer.append(sw.pden(sal_summer[-1],temp_summer[-1],sw.pres(0,lat))-1000)
        Sds.close()
        Tds.close()


sit_winter  = np.mean(sit_winter, axis = 0)
sal_winter  = np.mean(sal_winter, axis = 0)
temp_winter = np.mean(temp_winter,axis = 0)
dens_winter = np.mean(dens_winter,axis = 0)

sit_summer  = np.mean(sit_summer, axis = 0)
sal_summer  = np.mean(sal_summer, axis = 0)
temp_summer = np.mean(temp_summer,axis = 0)
dens_summer = np.mean(dens_summer,axis = 0)

sit_winter[sit_winter <0.1] = np.nan
sit_summer[sit_summer <0.1] = np.nan


################
### - PLOT - ###
################

#titles_list = ['SIT\nWinter', 'Salinity\nWinter',r'$\sigma_\theta$'+'\nWinter','Summer','Summer','Summer']
titles_list = ['Mean sea ice thickness\nWinter', 'Mean surface salinity\nWinter','Mean surface density anomaly'+'\nWinter','Summer','Summer','Summer']
label_list  = [r'$[m]$', r'$[psu]$', r'$[kg\cdot m^{-3}]$']
cmap_list =['cmo.ice','cmo.haline','cmo.dense','cmo.ice','cmo.haline','cmo.dense']
v = [[0,3],[30,35],[24,29],[0,3],[30,35],[24,29]]
titlesize = 30
labelsize = 30
fig,axes = plt.subplots(ncols = 3, nrows = 2, figsize = (23,13), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -7)})
plt.subplots_adjust(left=0.002,
                    bottom=0.05,
                    wspace=0.2,  
                    top=0.93)


for map,i in zip([sit_winter,sal_winter,dens_winter,sit_summer,sal_summer,dens_summer],range(6)):
    
    
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

    cs = axes.flatten()[i].pcolormesh(lon, lat, map,vmin = v[i][0], vmax = v[i][1],cmap = cmap_list[i], transform=ccrs.PlateCarree())

    #axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    if i > 2:
        cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.0,axes.flatten()[i].get_position().y0 + 0.05,0.01,axes.flatten()[i].get_position().height * 2])
        cb = plt.colorbar(cs, cax = cax,)
        cb.set_label(label = label_list[i-3],size = labelsize)
        cb.ax.tick_params(labelsize=labelsize)

plt.savefig(f'MT_plot/sit_s_dens.png')
plt.close()

