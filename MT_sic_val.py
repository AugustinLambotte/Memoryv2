import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import os
from ORAS5.useful_func_oras5 import plot, vector_plot, interpol_oras5_to_obs


lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

lonCS2 = np.loadtxt('C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/obs_data/Data/geo/longitude.txt')
latCS2 = np.loadtxt('C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/obs_data/Data/geo/latitude.txt')


ORAS5sic_march = []
CS2sic_march = []
ORAS5sic_sept  = [] 
CS2sic_sept  = [] 
for file in os.listdir('ORAS5/Data/iicethic'):
    if int(file[:4]) < 2011 or int(file[:4]) > 2018:
        continue
    if int(file[5:7]) == 3:
        ORAS5sic_march.append(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))
        CS2sic_march.append(np.loadtxt(f'obs_data/Data/SI/sic/month_mean/{file}'))
    if int(file[5:7]) == 9:
        ORAS5sic_sept.append(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))
        CS2sic_sept.append(np.loadtxt(f'obs_data/Data/SI/sic/month_mean/{file}'))

ORAS5sic_march = np.mean(ORAS5sic_march,axis = 0)
ORAS5sic_sept  = np.mean(ORAS5sic_sept,axis = 0)

CS2sic_march = np.mean(CS2sic_march,axis = 0)
CS2sic_sept  = np.mean(CS2sic_sept,axis = 0)

ORAS5sic_march[ORAS5sic_march <0.1] = np.nan
ORAS5sic_sept[ORAS5sic_sept <0.1] = np.nan
CS2sic_march[CS2sic_march <0.1] = np.nan
CS2sic_sept[CS2sic_sept <0.1] = np.nan

####################################################
##### - Interpolation of ORAS5 over CS2 grid - #####
####################################################

ORAS5sic_march_interp = interpol_oras5_to_obs(ORAS5sic_march)
ORAS5sic_sept_interp  = interpol_oras5_to_obs(ORAS5sic_sept)

diff_march = ORAS5sic_march_interp - CS2sic_march
diff_sept  = ORAS5sic_sept_interp  - CS2sic_sept

################
### - PLOT - ###
################

titles_list = ['ORAS5\nMarch', 'OSI-SAF\nMarch','Difference\nMarch', 'September', 'September','September']
titlesize = 32
labelsize = 28
fig,axes = plt.subplots(ncols = 3, nrows = 2, figsize = (23,13), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -7)})
plt.subplots_adjust(left=0.002,
                    bottom=0.02,
                    wspace=0.2,  
                    top=0.93)


for map,i in zip([ORAS5sic_march*1e+2,CS2sic_march*1e+2,diff_march*1e+2,ORAS5sic_sept*1e+2,CS2sic_sept*1e+2,diff_sept*1e+2],range(6)):
    
    
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

    if i in [0,3]:
        cs = axes.flatten()[i].pcolormesh(lon, lat, map,vmin = 0, vmax = 100,cmap = "cmo.ice", transform=ccrs.PlateCarree())
    elif i in [1,4]:
        cs = axes.flatten()[i].pcolormesh(lonCS2, latCS2, map,vmin = 0, vmax = 100,cmap = "cmo.ice", transform=ccrs.PlateCarree())
    else:
        cs = axes.flatten()[i].pcolormesh(lonCS2, latCS2, map,vmin = -10, vmax = 10,cmap = "cmo.balance", transform=ccrs.PlateCarree())


    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    
    if i == 4:
        cax = fig.add_axes([axes.flatten()[4].get_position().x1,axes.flatten()[4].get_position().y0+0.05 ,0.015,axes.flatten()[4].get_position().height * 2])
        cb = plt.colorbar(cs, cax = cax,)
        cb.set_label(label = '[%]',size = labelsize)
        cb.ax.tick_params(labelsize=labelsize)

cax = fig.add_axes([axes.flatten()[-1].get_position().x1,axes.flatten()[-1].get_position().y0+0.05 ,0.015,axes.flatten()[-1].get_position().height*2])
cb = plt.colorbar(cs, cax = cax,)
cb.set_label(label = '[%]',size = labelsize)
cb.ax.tick_params(labelsize=labelsize)
plt.savefig(f'MT_plot/validation/SIC_val')
plt.close()

