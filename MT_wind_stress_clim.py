import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import xarray as xr
from scipy import stats
import numpy.ma as ma
from ORAS5.useful_func_oras5 import plot, show, vector_plot
from datetime import date
earth_radius = 6370 * 1e3

u = []
v = []
for file in os.listdir(f'ORAS5/Data/vozocrte'):
    year = int(file[:4])
    month = int(file[5:7])

    u.append(np.loadtxt(f'ORAS5/Data/sozotaux/{file}'))
    v.append(np.loadtxt(f'ORAS5/Data/sometauy/{file}'))

u = np.mean(u,axis = 0)
v = np.mean(v,axis = 0)

lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

fig,axs = plt.subplots( figsize = (10,7), subplot_kw= {'projection' : ccrs.LambertConformal(central_longitude = -7)})
xlim = [-33, 20]
ylim = [65, 81]
plt.subplots_adjust(left=0.01,
                    bottom=0.05, 
                    top=0.93,)
lower_space = 3 
rect = mpath.Path([[xlim[0], ylim[0]],
                [xlim[1], ylim[0]],
                [xlim[1], ylim[1]],
                [xlim[0], ylim[1]],
                [xlim[0], ylim[0]],
                ]).interpolated(20)

proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axs) - axs.transData
rect_in_target = proj_to_data.transform_path(rect)
axs.set_boundary(rect_in_target)
axs.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

axs.coastlines()
axs.gridlines()



dlat = (v * 360)/(2*np.pi * earth_radius)
dlon = (u * 360)/(2*np.pi * earth_radius * np.cos(lat/360 * 2*np.pi))
#Convert from m/s to degree/s
current_intensity = np.sqrt(v**2 + u**2)
current_intensity[current_intensity == 0] = np.nan
axs.set_title('Wind stress', fontsize = 21)
cs = axs.pcolormesh(lon, lat, current_intensity,vmin = 0, vmax = 0.2,cmap = "cmo.speed", transform=ccrs.PlateCarree())

#Vector plot
current_intensity_vector_plot = np.sqrt(dlon**2 + dlat**2)
dlon /= (current_intensity_vector_plot*18)
dlat /= (current_intensity_vector_plot*18)
axs.quiver(np.array(lon),np.array(lat),np.array(dlon),np.array(dlat),transform = ccrs.PlateCarree())
cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.025,axs.get_position().height])
cb = plt.colorbar(cs, cax = cax)
cb.ax.tick_params(labelsize=15)
cb.set_label(label = r'$[N/m^2]$',size = 15)

plt.savefig('MT_plot/Wind_stress.png')
plt.close()
