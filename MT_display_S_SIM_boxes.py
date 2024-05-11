import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import cartopy.crs as ccrs
import gsw
import seawater as sw
import os
import matplotlib.path as mpath
import matplotlib.patches as mpatches

lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")
earth_radius = 6370 * 1e3
### - Mean EGC Path - ###
u = []
v = []

for file in os.listdir(f'ORAS5/Data/vozocrte'):
    year = int(file[:4])
    month = int(file[5:7])

    u.append(np.loadtxt(f'ORAS5/Data/vozocrte/{file}'))
    v.append(np.loadtxt(f'ORAS5/Data/vomecrtn/{file}'))


u = np.mean(u,axis = 0)
v = np.mean(v,axis = 0)


areaS   = [76.5,78.5,-13,-4]
areaSIM = [78.5,80.5,-11,-3]

areaSIM = np.where((lat >= areaSIM[0]) & (lat <= areaSIM[1]) & (lon >= areaSIM[2]) & (lon <= areaSIM[3]),1,0)
areaS   = np.where((lat >= areaS[0]) & (lat <= areaS[1]) & (lon >= areaS[2]) & (lon <= areaS[3]),1,0)

titlesize = 17
labelsize = 17
fig,axes = plt.subplots( figsize = (16,10), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -6)})

xlim = [-33, 18]
ylim = [73, 82]
lower_space = 3 
rect = mpath.Path([[xlim[0], ylim[0]],
                [xlim[1], ylim[0]],
                [xlim[1], ylim[1]],
                [xlim[0], ylim[1]],
                [xlim[0], ylim[0]],
                ]).interpolated(20)
proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axes) - axes.transData
rect_in_target = proj_to_data.transform_path(rect)
axes.set_boundary(rect_in_target)
axes.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

axes.coastlines()
axes.gridlines()        

### - Vector plot - ###

dlat = (v * 360)/(2*np.pi * earth_radius)
dlon = (u * 360)/(2*np.pi * earth_radius * np.cos(lat/360 * 2*np.pi))
#Convert from m/s to degree/s
current_intensity = np.sqrt(v**2 + u**2)
current_intensity[current_intensity == 0] = np.nan
cs = axes.pcolormesh(lon, lat, current_intensity,vmin = 0, vmax = 0.3,cmap = "cmo.speed", transform=ccrs.PlateCarree())
current_intensity_vector_plot = np.sqrt(dlon**2 + dlat**2)
dlon /= (current_intensity_vector_plot*18)
dlat /= (current_intensity_vector_plot*18)
axes.quiver(np.array(lon),np.array(lat),np.array(dlon),np.array(dlat),transform = ccrs.PlateCarree())

### - boxes - ###
axes.add_patch(mpatches.Rectangle(xy=[-13.5, 76], width=9, height=3,
                                    facecolor='red',
                                    alpha = 0.5,
                                    transform=ccrs.Geodetic()))
axes.add_patch(mpatches.Rectangle(xy=[-11.5, 78], width=8, height=3,
                                    facecolor='blue',
                                    alpha = 0.5,
                                    transform=ccrs.Geodetic()))
#axes.contour(lon,lat,areaS, levels = [0,1], colors = 'blue',interpolation = 'none',linewidths = 3,transform=ccrs.PlateCarree())
#axes.contour(lon,lat,areaSIM, levels = [0,1],  colors = 'red',interpolation = 'none',linewidths = 3,transform=ccrs.PlateCarree())

plt.savefig(f'MT_plot/S_SIM_boxes.png')
plt.close()
