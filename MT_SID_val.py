import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import os
from scipy import interpolate

lonORAS5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
latORAS5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")
earth_radius = 6370 * 1e3
# ORAS5
u_april  = []
v_april  = []

u_august = []
v_august = []

for file in os.listdir(f'ORAS5/Data/iicevelu'):
    year = int(file[:4])
    month = int(file[5:7])
    if year >= 2011 and year <= 2018:
        if month == 3:
            u_april.append(np.loadtxt(f'ORAS5/Data/iicevelu/{file}'))
            v_april.append(np.loadtxt(f'ORAS5/Data/iicevelv/{file}'))
        if month == 9:
            u_august.append(np.loadtxt(f'ORAS5/Data/iicevelu/{file}'))
            v_august.append(np.loadtxt(f'ORAS5/Data/iicevelv/{file}'))

u_april = np.array(u_april)
v_april = np.array(v_april)
u_august = np.array(u_august)
v_august = np.array(v_august)

u_april  = np.nan_to_num(u_april)
v_april  = np.nan_to_num(v_april)
u_august = np.nan_to_num(u_august)
v_august = np.nan_to_num(v_august)

u_aprilOras5  = np.mean(u_april,axis = 0)
v_aprilOras5  = np.mean(v_april,axis = 0)

u_augustOras5 = np.mean(u_august,axis = 0)
v_augustOras5 = np.mean(v_august,axis = 0)

u_aprilOras5[u_aprilOras5 == 0]   = np.nan
v_aprilOras5[v_aprilOras5 == 0]   = np.nan
u_augustOras5[u_augustOras5 == 0] = np.nan
v_augustOras5[v_augustOras5 == 0] = np.nan

# OSI-SAF
u_april  = []
v_april  = []

u_august = []
v_august = []
lat_min = 60
lat_max = 83
lon_min = -40
lon_max = 20
lonlat_already_acquired = False
for year in range(2011,2019):
    print(year)
    for file in os.listdir(f"C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/{year}/03"):
        if 'ice_drift_nh_ease2-750_cdr-v1p0_24h-' in file:
            fileApril = f"C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/{year}/03/{file}"
            ds = xr.open_dataset(fileApril, decode_times = False)

            if ~lonlat_already_acquired:
                lonOSISAF = np.nan_to_num(np.array(ds['lon'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)))
                latOSISAF = np.nan_to_num(np.array(ds['lat'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)))
                

                
                lonlat_already_acquired = True

            dlat = np.array((ds['lat1'] - ds['lat']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True))[0]
            dlon = np.array((ds['lon1'] - ds['lon']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True))[0]
            v_april.append((dlat/360 * 2*np.pi) * earth_radius /(24*60*60))                                  #[m/s]
            u_april.append((dlon/360 * 2*np.pi) * np.cos(latOSISAF/360 * 2*np.pi)* earth_radius /(24*60*60)) #[m/s]
            
            ds.close()
    for file in os.listdir(f"C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/{year}/09"):
        if 'ice_drift_nh_ease2-750_cdr-v1p0_24h-' in file:
            fileAugust = f"C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/{year}/09/{file}"
        
            ds = xr.open_dataset(fileAugust, decode_times = False)
        
            dlat = np.array((ds['lat1'] - ds['lat']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True))[0]
            dlon = np.array((ds['lon1'] - ds['lon']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True))[0]

            v_august.append((dlat/360 * 2*np.pi) * earth_radius /(24*60*60))                                  #[m/s]
            u_august.append((dlon/360 * 2*np.pi) * np.cos(latOSISAF/360 * 2*np.pi)* earth_radius /(24*60*60)) #[m/s]

            ds.close()


u_april = np.array(u_april)
v_april = np.array(v_april)
u_august = np.array(u_august)
v_august = np.array(v_august)

u_april  = np.nan_to_num(u_april)
v_april  = np.nan_to_num(v_april)
u_august = np.nan_to_num(u_august)
v_august = np.nan_to_num(v_august)

u_aprilCS2  = np.mean(u_april,axis = 0)
v_aprilCS2  = np.mean(v_april,axis = 0)

u_augustCS2 = np.mean(u_august,axis = 0)
v_augustCS2 = np.mean(v_august,axis = 0)

u_aprilCS2[u_aprilCS2 == 0]   = np.nan
v_aprilCS2[v_aprilCS2 == 0]   = np.nan
u_augustCS2[u_augustCS2 == 0] = np.nan
v_augustCS2[v_augustCS2 == 0] = np.nan



u_aprilCS2[2,:]    = np.nan
v_aprilCS2[2,:]    = np.nan
u_augustCS2[2,:]   = np.nan
v_augustCS2[2,:]   = np.nan
u_aprilCS2[3,:]    = np.nan
v_aprilCS2[3,:]    = np.nan
u_augustCS2[3,:]   = np.nan
v_augustCS2[3,:]   = np.nan
u_aprilCS2[4,:]    = np.nan
v_aprilCS2[4,:]    = np.nan
u_augustCS2[4,:]   = np.nan
v_augustCS2[4,:]   = np.nan
u_aprilCS2[11,34]  = np.nan
v_aprilCS2[11,34]  = np.nan
u_augustCS2[11,34] = np.nan
v_augustCS2[11,34] = np.nan


##############################################
#####- Interpolation over OSI-SAF grid - #####
##############################################

absVelocity_OSISAF_august = np.sqrt(u_augustCS2**2   + v_augustCS2**2  )
absVelocity_OSISAF_april  = np.sqrt(u_aprilCS2**2    + v_aprilCS2**2   )
absVelocity_ORAS5_august  = np.sqrt(u_augustOras5**2 + v_augustOras5**2)
absVelocity_ORAS5_april   = np.sqrt(u_aprilOras5**2  + v_aprilOras5**2 )

points = [] # list of length NxM containing all the coordinates [lat,lon] of all points from si drift map
values_april  = []
values_august = []

#data[data == 0] = np.nan
for i in range(len(latORAS5)):
    for j in range(len(lonORAS5[0])):
        #if (data[i,j] !=0)  and (not np.isnan(data[i,j])):
        points.append([latORAS5[i,j],lonORAS5[i,j]])
        values_april.append(absVelocity_ORAS5_april[i,j])
        values_august.append(absVelocity_ORAS5_august[i,j])

points = np.array(points)
values_april  = np.array(values_april)
values_august = np.array(values_august)

april_interp  = interpolate.griddata(points, values_april,  (latOSISAF, lonOSISAF), method='nearest')
august_interp = interpolate.griddata(points, values_august, (latOSISAF, lonOSISAF), method='nearest')

diff_april  = april_interp - absVelocity_OSISAF_april
diff_august = august_interp - absVelocity_OSISAF_august


title_list = ['March\nORAS5','September\nORAS5','OSI-SAF','OSI-SAF','Difference','Difference']



fig,axes = plt.subplots(nrows = 3, ncols = 2,figsize = (16,14), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -6)})
titlesize = 25
labelsize = 25
plt.subplots_adjust(left=0.01,
                    bottom=0.05,
                    wspace=-0.3,  
                    top=0.93,)

for u,v,i in zip([u_aprilOras5, u_augustOras5, u_aprilCS2,u_augustCS2], [v_aprilOras5, v_augustOras5, v_aprilCS2,v_augustCS2], range(4)):
    if i in [0,1]:
        dlat = (v * 360)/(2*np.pi * earth_radius)
        dlon = (u * 360)/(2*np.pi * earth_radius * np.cos(latORAS5/360 * 2*np.pi))
    if i in [2,3]:
        dlat = (v * 360)/(2*np.pi * earth_radius)
        dlon = (u * 360)/(2*np.pi * earth_radius * np.cos(latOSISAF/360 * 2*np.pi))
           
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

    #Convert from m/s to degree/s
    current_intensity = np.sqrt(v**2 + u**2)
    current_intensity[current_intensity == 0] = np.nan
    
    if i in [0,1]:
        cs = axes.flatten()[i].pcolormesh(lonORAS5, latORAS5,   current_intensity, vmin = 0, vmax = 0.25,cmap = "cmo.speed", transform=ccrs.PlateCarree())
    if i in [2,3]:
        cs = axes.flatten()[i].pcolormesh(lonOSISAF, latOSISAF, current_intensity, vmin = 0, vmax = 0.25,cmap = "cmo.speed", transform=ccrs.PlateCarree())
    
    #Vector plot
    current_intensity_vector_plot = np.sqrt(dlon**2 + dlat**2)
    dlon /= (current_intensity_vector_plot*18)
    dlat /= (current_intensity_vector_plot*18)
    if i in [0,1]:
        axes.flatten()[i].quiver(np.array(lonORAS5),np.array(latORAS5),np.array(dlon),np.array(dlat),transform = ccrs.PlateCarree())
    if i in [2,3]:
        
        axes.flatten()[i].quiver(np.array(lonOSISAF),np.array(latOSISAF),np.array(dlon),np.array(dlat),transform = ccrs.PlateCarree())

    axes.flatten()[i].set_title(title_list[i], fontsize = titlesize)
    
cax = fig.add_axes([axes.flatten()[3].get_position().x1+0.02,axes.flatten()[3].get_position().y0,0.02,axes.flatten()[3].get_position().height * 2])
cb = plt.colorbar(cs, cax = cax,)
cb.set_label(label = '[m/s]',size = labelsize)
cb.ax.tick_params(labelsize=labelsize)


##############################
##### - Plot difference - ####
##############################
for map,i in zip([diff_april,diff_august], range(4,6)):
    
           
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
   
    cs = axes.flatten()[i].pcolormesh(lonOSISAF, latOSISAF, map, vmin = -0.1, vmax = 0.1,cmap = "cmo.balance", transform=ccrs.PlateCarree())
        
    axes.flatten()[i].set_title(title_list[i], fontsize = titlesize)
    
cax = fig.add_axes([axes.flatten()[-1].get_position().x1+0.02,axes.flatten()[-1].get_position().y0 ,0.02,axes.flatten()[-1].get_position().height])
cb = plt.colorbar(cs, cax = cax,)
cb.set_label(label = '[m/s]',size = labelsize)
cb.ax.tick_params(labelsize=labelsize)
plt.savefig('MT_plot/validation/SID_val.png')