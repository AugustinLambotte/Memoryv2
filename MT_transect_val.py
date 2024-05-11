import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
from ORAS5.useful_func_oras5 import plot, plot_bathy
import cartopy.crs as ccrs
import gsw
import seawater as sw
import matplotlib.path as mpath

import os



month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July','August','September','October','November','December']
grid_transect1 = [[77.5,-8 ],[77.5,-7 ],[77.5,-6 ],[77.5,-5 ],[77.5,-4 ],[77.5,-3 ],[77.5,-2],[77.5,-1]]
grid_transect2  = [[75.5,-13],[75.5,-12],[75.5,-11],[75.5,-10],[75.5,-9 ],[75.5,-8 ],[75.5,-7 ]]
grid_transect3  = [[73.5,-18],[73.5,-17],[73.5,-16],[73.5,-15],[73.5,-14],[73.5,-13]]
grid_transect4  = [[71.5,-21],[71.5,-20],[71.5,-19],[71.5,-18],[71.5,-17]]
grid_transect5  = [[69.5,-21],[69.5,-20],[69.5,-19],[69.5,-18],[69.5,-17]]
grid_transect6  = [[68.5,-24],[68.5,-23],[68.5,-22]]
grid_transect7  = [[67.5,-26],[67.5,-25],[67.5,-24],[67.5,-23],[67.5,-22]]

transect = [grid_transect1,grid_transect2,grid_transect3,grid_transect4,grid_transect5,grid_transect6,grid_transect7]
transectName = ['Transect 1', 'Transect 2','Transect 3','Transect 4','Transect 5','Transect 6','Transect 7']

###########################################
## - Plotting bathyemtry with transect - ##
###########################################
    
lon_oras5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat_oras5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

""" fig = plt.figure(figsize=(12,8))
axs = plt.axes(projection = ccrs.LambertConformal(central_longitude = -6))            
xlim = [-33, 20]
ylim = [65, 80]
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
gate_map = np.zeros(np.shape(lon_oras5))
for k in range(len(transect)):
    # Adding line
    geodetic = ccrs.Geodetic()
    ad_lon_t, ad_lat_t = ccrs.PlateCarree().transform_point(transect[k][0][1], transect[k][0][0], geodetic)
    liv_lon_t, liv_lat_t = ccrs.PlateCarree().transform_point(transect[k][-1][1], transect[k][-1][0], geodetic)

    axs.plot([ad_lon_t, liv_lon_t], [ad_lat_t, liv_lat_t],
            color='green', linewidth=2, marker='o', markersize=2,
            transform=ccrs.PlateCarree())
    
    for i in range(len(lat_oras5)):
        for j in range(len(lat_oras5[0])):
            if [lat_oras5[i,j],lon_oras5[i,j]] in transect[k]:
                #Colouring grid cell
                gate_map[i,j] = 1e6

                
gate_map[gate_map == 0] = np.nan
axs.coastlines()
axs.gridlines()
ds_bathy = xr.open_dataset('C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/ETOPO_2022_v1_60s_N90W180_surface.nc')
lat_bathy = np.array(ds_bathy['lat'].where((ds_bathy.lat < 81) & (ds_bathy.lat > 60) & (ds_bathy.lon> -40) & (ds_bathy.lon < 20)))
lon_bathy = np.array(ds_bathy['lon'].where((ds_bathy.lat < 81) & (ds_bathy.lat > 60) & (ds_bathy.lon> -40) & (ds_bathy.lon < 20))).T
z         = np.array(ds_bathy['z'].where(  (ds_bathy.lat < 81) & (ds_bathy.lat > 60) & (ds_bathy.lon> -40) & (ds_bathy.lon < 20)))
z         = z[        9000:10200,8430:11990]
lat_bathy = lat_bathy[9000:10200,8430:11990]
lon_bathy = lon_bathy[9000:10200,8430:11990]

lat_bathy = np.nan_to_num(lat_bathy)
lon_bathy = np.nan_to_num(lon_bathy)
z = np.nan_to_num(z)
z[z>0] = 0
z *= -1

csBathy = axs.pcolormesh(lon_bathy, lat_bathy, z,cmap = "cmo.dense", transform=ccrs.PlateCarree())
csGate  = axs.pcolormesh(lon_oras5, lat_oras5, gate_map,cmap = "cmo.balance",vmax = 1,vmin = 0, transform=ccrs.PlateCarree())
cax = fig.add_axes([axs.get_position().x1+0.01,axs.get_position().y0 - 0.02,0.04,axs.get_position().height])
cb = plt.colorbar(csBathy, cax = cax)
cb.set_label('Depth [m]', fontsize = 20)
ds_bathy.close()
      

#axs.set_title(title, fontsize = titlesize)
cb.ax.tick_params(labelsize=17)

plt.savefig('MT_plot/bathy_transect.png')
plt.close() """

#######################################################
#### - Plotting winter and summer mean vert prof - ####
#######################################################

# Plot S,T,sigma mean over the 12 months
SProfile_list_winter     = [[],[],[],[],[],[],[]]
TProfile_list_winter     = [[],[],[],[],[],[],[]]
sigmaProfile_list_winter = [[],[],[],[],[],[],[]]

SProfile_list_summer     = [[],[],[],[],[],[],[]]
TProfile_list_summer     = [[],[],[],[],[],[],[]]
sigmaProfile_list_summer = [[],[],[],[],[],[],[]]
for file in os.listdir('ORAS5/Data/votemper'):
    if int(file[5:7]) in [11,12,1,2,3]: #Winter
        dsS = xr.open_dataset(f"ORAS5/Data/vosaline/{file}")
        dsT = xr.open_dataset(f"ORAS5/Data/votemper/{file}")
        depth = np.array(dsS['deptht'])

        for transect_, transectName_,i in zip(transect,transectName,range(8)):
            SProfile = []
            TProfile = []
            sigmaProfile = []
            for point in range(len(transect_)):
                SProfile.append(np.array(dsS['vosaline'].sel(lon = transect_[point][1], lat = transect_[point][0], method = 'nearest')))
                TProfile.append(np.array(dsT['votemper'].sel(lon = transect_[point][1], lat = transect_[point][0], method = 'nearest')))
                sigmaProfile.append(sw.pden(SProfile[-1],TProfile[-1],sw.pres(depth,transect_[point][0]))-1000)
            SProfile_list_winter [i].append(np.array(SProfile))
            TProfile_list_winter [i].append(np.array(TProfile))
            sigmaProfile_list_winter [i].append(np.array(sigmaProfile))
        dsS.close()
        dsT.close()
    if int(file[5:7]) in [5,6,7,8,9]: #Winter
        dsS = xr.open_dataset(f"ORAS5/Data/vosaline/{file}")
        dsT = xr.open_dataset(f"ORAS5/Data/votemper/{file}")
        depth = np.array(dsS['deptht'])

        for transect_, transectName_,i in zip(transect,transectName,range(8)):
            SProfile = []
            TProfile = []
            sigmaProfile = []
            for point in range(len(transect_)):
                SProfile.append(np.array(dsS['vosaline'].sel(lon = transect_[point][1], lat = transect_[point][0], method = 'nearest')))
                TProfile.append(np.array(dsT['votemper'].sel(lon = transect_[point][1], lat = transect_[point][0], method = 'nearest')))
                sigmaProfile.append(sw.pden(SProfile[-1],TProfile[-1],sw.pres(depth,transect_[point][0]))-1000)
            SProfile_list_summer[i].append(np.array(SProfile))
            TProfile_list_summer[i].append(np.array(TProfile))
            sigmaProfile_list_summer[i].append(np.array(sigmaProfile))
        dsS.close()
        dsT.close()


        

##############################
########## - Plot - ##########
##############################
for i in range(len(SProfile_list_summer)):
    lonList = np.array(transect[i])[:,1]
    SProfile_winter          = np.mean(SProfile_list_winter[i]    ,axis = 0)
    TProfile_winter          = np.mean(TProfile_list_winter[i]    ,axis = 0)
    sigmaProfile_winter      = np.mean(sigmaProfile_list_winter[i],axis = 0)

    SProfile_summer          = np.mean(SProfile_list_summer[i]    ,axis = 0)
    TProfile_summer          = np.mean(TProfile_list_summer[i]    ,axis = 0)
    sigmaProfile_summer      = np.mean(sigmaProfile_list_summer[i],axis = 0)

    fig, axes = plt.subplots(nrows = 3, ncols = 2, figsize = (15,9), sharex = True, sharey = True)
    plt.subplots_adjust(left=0.1, 
                    right=0.99, 
                    wspace=0.1)    
    fig.suptitle(f'{transectName[i]}',fontsize = 20)

    labelsize = 18
    titlesize = 18
    ticklabelsize = 15
    legendsize = 15 
    cs = axes[0,0].pcolormesh(lonList,depth,SProfile_winter.T, vmin = 33.5, vmax = 35.1, cmap = 'cmo.haline')
    axes[0,0].set_title('Winter \nSalinity [psu]',fontsize = titlesize)
    #axes[0,0].set_xlabel('Longitude [°]',size = labelsize)
    axes[0,0].set_ylabel('Depth [m]',size = labelsize)
    axes[0,0].tick_params(labelsize = ticklabelsize)
    axes[0,0].set_ylim(300,0)
    axes[0,0].grid()
    cb = plt.colorbar(cs)
    cb.ax.tick_params(labelsize=ticklabelsize)

    cs = axes[1,0].pcolormesh(lonList,depth,TProfile_winter.T, vmin = -2.5, vmax = 3, cmap = 'cmo.thermal')
    axes[1,0].set_title('Potential temperature [°C]',fontsize = titlesize)
    #axes[1,0].set_xlabel('Longitude [°]',size = labelsize)
    axes[1,0].set_ylabel('Depth [m]',size = labelsize)
    axes[1,0].tick_params(labelsize = ticklabelsize)
    axes[1,0].set_ylim(300,0)
    axes[1,0].grid()

    cb = plt.colorbar(cs)
    cb.ax.tick_params(labelsize=ticklabelsize)

    cs = axes[2,0].pcolormesh(lonList,depth,sigmaProfile_winter.T, vmin = 26, vmax = 28.5, cmap = 'cmo.dense')
    axes[2,0].set_title(r'$\sigma_{\theta} [kg/m^3]$' ,fontsize = titlesize)
    axes[2,0].set_xlabel('Longitude [°]',size = labelsize)
    axes[2,0].set_ylabel('Depth [m]',size = labelsize)
    axes[2,0].tick_params(labelsize = ticklabelsize)
    axes[2,0].set_ylim(300,0)
    axes[2,0].grid()
    cb = plt.colorbar(cs)
    cb.ax.tick_params(labelsize=ticklabelsize)

    cs = axes[0,1].pcolormesh(lonList,depth,SProfile_summer.T, vmin = 33.5, vmax = 35.1, cmap = 'cmo.haline')
    axes[0,1].set_title('Summer \nSalinity [psu]',fontsize = titlesize)
    #axes[0,1].set_xlabel('Longitude [°]',size = labelsize)
    #axes[0,1].set_ylabel('Depth [m]',size = labelsize)
    axes[0,1].tick_params(labelsize = ticklabelsize)
    axes[0,1].set_ylim(300,0)
    axes[0,1].grid()
    cb = plt.colorbar(cs)
    cb.ax.tick_params(labelsize=ticklabelsize)

    cs = axes[1,1].pcolormesh(lonList,depth,TProfile_summer.T, vmin = -2.5, vmax = 3, cmap = 'cmo.thermal')
    axes[1,1].set_title('Potential temperature [°C]',fontsize = titlesize)
    #axes[1,1].set_xlabel('Longitude [°]',size = labelsize)
    axes[1,1].tick_params(labelsize = ticklabelsize)
    axes[1,1].set_ylim(300,0)
    axes[1,1].grid()

    cb = plt.colorbar(cs)
    cb.ax.tick_params(labelsize=ticklabelsize)

    cs = axes[2,1].pcolormesh(lonList,depth,sigmaProfile_summer.T, vmin = 26, vmax = 28.5, cmap = 'cmo.dense')
    axes[2,1].set_title(r'$\sigma_{\theta} [kg/m^3]$' ,fontsize = titlesize)
    axes[2,1].set_xlabel('Longitude [°]',size = labelsize)
    axes[2,1].tick_params(labelsize = ticklabelsize)
    axes[2,1].set_ylim(300,0)
    axes[2,1].grid()
    cb = plt.colorbar(cs)
    cb.ax.tick_params(labelsize=ticklabelsize)
    
    plt.savefig(f'MT_plot/Profile/{transectName[i]}.png')

