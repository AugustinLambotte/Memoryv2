import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import xarray as xr
from scipy import stats
from ORAS5.useful_func_oras5 import plot, show, vector_plot
from datetime import date, timedelta

earth_radius = 6370 * 1e3

""" ########################
#### - Copernicus - ####
########################

lat_min = 60
lat_max = 83
lon_min = -40
lon_max = 20

uCopernicus = []
vCopernicus = []

uCopernicus_winter = []
vCopernicus_winter = []

uCopernicus_summer = []
vCopernicus_summer = []


v_file = "C:/Users/Augustin/Downloads/Northward_sea_water_velocity"
u_file = "C:/Users/Augustin/Downloads/Eastward_sea_water_velocity"

v_ds = xr.open_dataset(v_file, decode_times = False)
u_ds = xr.open_dataset(u_file, decode_times = False)
    
u_gos= u_ds['ugos'].where((u_ds.longitude > lon_min) & (u_ds.longitude < lon_max) & (u_ds.latitude > lat_min) & (u_ds.latitude > 65.4 + (76.5-65.4)/(9+17) * (u_ds.longitude + 17)) & (u_ds.latitude < lat_max),drop = True)
v_gos= v_ds['vgos'].where((v_ds.longitude > lon_min) & (v_ds.longitude < lon_max) & (v_ds.latitude > lat_min) & (v_ds.latitude > 65.4 + (76.5-65.4)/(9+17) * (v_ds.longitude + 17)) & (v_ds.latitude < lat_max),drop = True)

for time in u_gos.time:
    time = int(time)
    corresponding_date = date(2010,10,1) + timedelta(days = time - int(u_gos.time[0]))
    if corresponding_date.year >= 2011 and corresponding_date.year <= 2018:
        print(corresponding_date)
        uCopernicus.append(np.array(u_gos.sel(time =time)))
        vCopernicus.append(np.array(v_gos.sel(time =time)))

        if corresponding_date.month in [1,2,3]:
            uCopernicus_winter.append(np.array(u_gos.sel(time =time)))
            vCopernicus_winter.append(np.array(v_gos.sel(time =time)))
        if corresponding_date.month in [7,8,9]:
            uCopernicus_summer.append(np.array(u_gos.sel(time =time)))
            vCopernicus_summer.append(np.array(v_gos.sel(time =time)))

lonCopernicus = u_ds.coords['longitude']
latCopernicus = u_ds.coords['latitude']

lonCopernicus, latCopernicus = np.meshgrid(lonCopernicus[:-1],latCopernicus[:-1])

u_ds.close()
v_ds.close()


uCopernicus = np.mean(np.nan_to_num(uCopernicus),axis = 0)
vCopernicus = np.mean(np.nan_to_num(vCopernicus),axis = 0)

uCopernicus_winter = np.mean(np.nan_to_num(uCopernicus_winter),axis = 0)
vCopernicus_winter = np.mean(np.nan_to_num(vCopernicus_winter),axis = 0)

uCopernicus_summer = np.mean(np.nan_to_num(uCopernicus_summer),axis = 0)
vCopernicus_summer = np.mean(np.nan_to_num(vCopernicus_summer),axis = 0)

uCopernicus[uCopernicus == 0] = np.nan
vCopernicus[vCopernicus == 0] = np.nan
uCopernicus_winter[uCopernicus_winter == 0] = np.nan
vCopernicus_winter[vCopernicus_winter == 0] = np.nan
uCopernicus_summer[uCopernicus_summer == 0] = np.nan
vCopernicus_summer[vCopernicus_summer == 0] = np.nan
 """
###################
#### - ORAS5 - ####
###################
uORAS5 = []
vORAS5 = []

uORAS5_winter = []
vORAS5_winter = []

uORAS5_summer = []
vORAS5_summer = []
for file in os.listdir(f'ORAS5/Data/vozocrte'):
    year = int(file[:4])
    month = int(file[5:7])
    
    if year in [i for i in range(2011,2019)]:
        uORAS5.append(np.loadtxt(f'ORAS5/Data/vozocrte/{file}'))
        vORAS5.append(np.loadtxt(f'ORAS5/Data/vomecrtn/{file}'))
         
        if month == 3:
            uORAS5_winter.append(np.loadtxt(f'ORAS5/Data/vozocrte/{file}'))
            vORAS5_winter.append(np.loadtxt(f'ORAS5/Data/vomecrtn/{file}'))
        
        if month == 9:
            uORAS5_summer.append(np.loadtxt(f'ORAS5/Data/vozocrte/{file}'))
            vORAS5_summer.append(np.loadtxt(f'ORAS5/Data/vomecrtn/{file}'))

uORAS5 = np.mean(uORAS5,axis = 0)
vORAS5 = np.mean(vORAS5,axis = 0)

uORAS5_winter = np.mean(uORAS5_winter,axis = 0)
vORAS5_winter = np.mean(vORAS5_winter,axis = 0)

uORAS5_summer = np.mean(uORAS5_summer,axis = 0)
vORAS5_summer = np.mean(vORAS5_summer,axis = 0)

lonORAS5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
latORAS5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

##################
#### - PLOT - ####
##################

titles_list = ['March', 'September']
titlesize = 30
labelsize = 25
fig,axes = plt.subplots(ncols = 2, figsize = (20,10), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -7)})
plt.subplots_adjust(left=0.01,
                    bottom=0.05,
                    right = 0.8,
                    wspace =0,  
                    top=0.93,)


for u,v,i in zip([uORAS5_winter,uORAS5_summer],[vORAS5_winter,vORAS5_summer],range(2)):
    
    u = np.array(u)
    v = np.array(v)

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

    dlat = (v * 360)/(2*np.pi * earth_radius)
    dlon = (u * 360)/(2*np.pi * earth_radius * np.cos(latORAS5/360 * 2*np.pi))
    #Convert from m/s to degree/s
    current_intensity = np.sqrt(v**2 + u**2)
    current_intensity[current_intensity == 0] = np.nan
    cs = axes.flatten()[i].pcolormesh(lonORAS5, latORAS5, current_intensity,vmin = 0, vmax = 0.2,cmap = "cmo.speed", transform=ccrs.PlateCarree())

    #Vector plot
    current_intensity_vector_plot = np.sqrt(dlon**2 + dlat**2)
    dlon /= (current_intensity_vector_plot*18)
    dlat /= (current_intensity_vector_plot*18)
    axes.flatten()[i].quiver(np.array(lonORAS5),np.array(latORAS5),np.array(dlon),np.array(dlat),transform = ccrs.PlateCarree())


    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    
cax = fig.add_axes([axes.flatten()[-1].get_position().x1+0.05,axes.flatten()[-1].get_position().y0 + 0.05,0.03,axes.flatten()[-1].get_position().height])
cb = plt.colorbar(cs, cax = cax,)
cb.set_label(label = '[m/s]',size = labelsize)
cb.ax.tick_params(labelsize=labelsize)

plt.savefig(f'MT_plot/validation/EGC_validation.png')
plt.close()
