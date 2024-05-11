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


ORAS5sit_march = []
CS2sit_march = []
ORAS5sit_sept  = [] 
CS2sit_sept  = [] 
for file in os.listdir('ORAS5/Data/iicethic'):
    if int(file[:4]) < 2011 or int(file[:4]) > 2018:
        continue
    if int(file[5:7]) == 3:
        ORAS5sit_march.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        CS2sit_march.append(np.loadtxt(f'obs_data/Data/SI/sit/month_mean/{file}'))
    if int(file[5:7]) == 9:
        ORAS5sit_sept.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        CS2sit_sept.append(np.loadtxt(f'obs_data/Data/SI/sit/month_mean/{file}'))

ORAS5sit_march = np.mean(ORAS5sit_march,axis = 0)
ORAS5sit_sept  = np.mean(ORAS5sit_sept,axis = 0)

CS2sit_march = np.mean(CS2sit_march,axis = 0)
CS2sit_sept  = np.mean(CS2sit_sept,axis = 0)

ORAS5sit_march[ORAS5sit_march <0.1] = np.nan
ORAS5sit_sept[ORAS5sit_sept <0.1] = np.nan
CS2sit_march[CS2sit_march <0.1] = np.nan
CS2sit_sept[CS2sit_sept <0.1] = np.nan

####################################################
##### - Interpolation of ORAS5 over CS2 grid - #####
####################################################

ORAS5sit_march_interp = interpol_oras5_to_obs(ORAS5sit_march)
ORAS5sit_sept_interp  = interpol_oras5_to_obs(ORAS5sit_sept)

diff_march = ORAS5sit_march_interp - CS2sit_march
diff_sept  = ORAS5sit_sept_interp  - CS2sit_sept

################
### - PLOT - ###
################

titles_list = ['ORAS5\nMarch', 'CS2\nMarch','Difference\nMarch', 'September', 'September','September']
titlesize = 32
labelsize = 28
fig,axes = plt.subplots(ncols = 3, nrows = 2, figsize = (23,13), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -7)})
plt.subplots_adjust(left=0.002,
                    bottom=0.02,
                    wspace=0.2,  
                    top=0.93)


for map,i in zip([ORAS5sit_march,CS2sit_march,diff_march,ORAS5sit_sept,CS2sit_sept,diff_sept],range(6)):
    
    
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
        cs = axes.flatten()[i].pcolormesh(lon, lat, map,vmin = 0, vmax = 3,cmap = "cmo.ice", transform=ccrs.PlateCarree())
    elif i in [1,4]:
        cs = axes.flatten()[i].pcolormesh(lonCS2, latCS2, map,vmin = 0, vmax = 3,cmap = "cmo.ice", transform=ccrs.PlateCarree())
    else:
        cs = axes.flatten()[i].pcolormesh(lonCS2, latCS2, map,vmin = -1.5, vmax = 1.5,cmap = "cmo.balance", transform=ccrs.PlateCarree())


    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    
    if i == 4:
        cax = fig.add_axes([axes.flatten()[4].get_position().x1,axes.flatten()[4].get_position().y0+0.05 ,0.015,axes.flatten()[4].get_position().height * 2])
        cb = plt.colorbar(cs, cax = cax,)
        cb.set_label(label = '[m]',size = labelsize)
        cb.ax.tick_params(labelsize=labelsize)

cax = fig.add_axes([axes.flatten()[-1].get_position().x1,axes.flatten()[-1].get_position().y0+0.05 ,0.015,axes.flatten()[-1].get_position().height*2])
cb = plt.colorbar(cs, cax = cax,)
cb.set_label(label = '[m]',size = labelsize)
cb.ax.tick_params(labelsize=labelsize)
plt.savefig(f'MT_plot/validation/SIT_val')
plt.close()


####################################
####### - Plot pré-défense - #######
####################################

titles_list = ['March', 'September']
titlesize = 25
labelsize = 25
fig,axes = plt.subplots(ncols = 2, figsize = (16,8), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -7)})
plt.subplots_adjust(left=0.01,
                    right = 0.85,
                    wspace = 0.05)


for map,i in zip([ORAS5sit_march,ORAS5sit_sept],range(2)):
    
    
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

    cs = axes.flatten()[i].pcolormesh(lon, lat, map,vmin = 0, vmax = 3,cmap = "cmo.ice", transform=ccrs.PlateCarree())


    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    
    
cax = fig.add_axes([axes.flatten()[-1].get_position().x1+0.05,axes.flatten()[-1].get_position().y0 ,0.015,axes.flatten()[-1].get_position().height])
cb = plt.colorbar(cs, cax = cax,)
cb.set_label(label = '[m]',size = labelsize)
cb.ax.tick_params(labelsize=labelsize)
plt.savefig(f'MT_plot/pre_def/SIT_cycle')
plt.close()

