import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import os
from ORAS5.useful_func_oras5 import plot, vector_plot, interpol_obs_to_oras5

earth_radius = 6370 * 1e3
month_name = ['Ja','Fe','Mar','Ap','May','Ju','July','Au','Se','Oc','No','De']

lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

lonCS2 = np.loadtxt('C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/obs_data/Data/geo/longitude.txt')
latCS2 = np.loadtxt('C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/obs_data/Data/geo/latitude.txt')

#######################
##### - SIT MAP - #####
#######################
ORAS5sit_march = []
CS2sit_march = []
ORAS5sit_sept  = [] 
CS2sit_sept  = [] 
for file in os.listdir('ORAS5/Data/iicethic'):
    if int(file[:4]) < 2011 or int(file[:4]) > 2018:
        continue
    if int(file[5:7]) == 4:
        ORAS5sit_march.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        CS2sit_march.append(np.loadtxt(f'obs_data/Data/SI/sit/month_mean/{file}'))
    if int(file[5:7]) == 8:
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


#######################################
############# - SIV SIE - #############
#######################################


#################################
########### - ORAS5 - ###########
#################################
sit_oras5 = []
sic_oras5 = []

EGC_path_oras5 = np.nan_to_num(np.loadtxt('ORAS5/Data/geo/EGC_path.txt'))
surface_oras5  = np.loadtxt('ORAS5/Data/geo/surface.txt')
for file in os.listdir('ORAS5/Data/ke'):
    
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        sit_oras5.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        sic_oras5.append(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))

sit_oras5 = np.array(sit_oras5)
sic_oras5 = np.array(sic_oras5)
sit_oras5[:,EGC_path_oras5 == 0] = np.nan
sic_oras5[:,EGC_path_oras5 == 0] = np.nan

#sit_oras5[(sit_oras5 == 0) | (sic_oras5 < 0.15)] = np.nan
#sic_oras5[(sic_oras5 == 0) | (sic_oras5 < 0.15)] = np.nan

#sit_oras5[(sic_oras5 < 0.15)] = 0
#sic_oras5[(sic_oras5 < 0.15)] = 0

siv_oras5 = []
sie_oras5 = []
for i in range(len(sit_oras5)):
    siv_oras5.append(np.nansum(sit_oras5[i] * sic_oras5[i] * surface_oras5) *1e-9)
    sie_oras5.append(np.nansum(np.where(sic_oras5[i]>0.15,surface_oras5,np.nan))*1e-6)
siv_oras5 = np.array(siv_oras5)
sie_oras5 = np.array(sie_oras5)



sit_oras5 = np.nanmean(sit_oras5,axis = (1,2))
sic_oras5 = np.nanmean(sic_oras5,axis = (1,2))*1e2

siv_oras5_ra = np.convolve(siv_oras5, np.ones(12)/12, mode = 'valid')
sie_oras5_ra = np.convolve(sie_oras5, np.ones(12)/12, mode = 'valid')




###############################
########### - OBS - ###########
###############################
sit_obs = []
sic_obs = []

EGC_path_obs = np.nan_to_num(np.loadtxt('obs_data/Data/geo/EGC_path.txt'))
surface_obs  = np.loadtxt('obs_data/Data/geo/surface.txt')
print(f'surface ORAS 5 - CS2 = {(np.nansum(np.where(EGC_path_oras5 == 1, surface_oras5,0)) - np.nansum(np.where(EGC_path_obs == 1, surface_obs,0)))*1e-6} km^2')
for file in os.listdir('obs_data/Data/gos/ke/month_mean'):
    
    if int(file[:4]) >= 2011 and int(file[:4]) < 2019:
        sit_obs.append(np.loadtxt(f'obs_data/Data/SI/sit/month_mean/{file}' ))
        sic_obs.append(np.loadtxt(f'obs_data/Data/SI/sic/month_mean/{file}' ))

sit_obs = np.array(sit_obs)
sic_obs = np.array(sic_obs)

sit_obs[:,EGC_path_obs == 0] = np.nan
sic_obs[:,EGC_path_obs == 0] = np.nan

#sit_obs[(sit_obs == 0)|(sic_obs <0.15)] = np.nan
#sic_obs[(sic_obs == 0)|(sic_obs <0.15)] = np.nan

#sit_obs[(sic_obs <0.15)] = 0
#sic_obs[(sic_obs <0.15)] = 0

siv_obs = []
sie_obs = []
for i in range(len(sit_obs)):
    siv_obs.append(np.nansum(sit_obs[i] * sic_obs[i] * surface_obs) *1e-9)
    sie_obs.append(np.nansum(np.where(sic_obs[i]>0.15,surface_obs,np.nan))*1e-6)
siv_obs = np.array(siv_obs)
sie_obs = np.array(sie_obs)

siv_obs_ra = np.convolve(siv_obs, np.ones(12)/12, mode = 'valid')
sie_obs_ra = np.convolve(sie_obs, np.ones(12)/12, mode = 'valid')

#############################
########## - SID - ##########
#############################

# ORAS5
u_april  = []
v_april  = []

u_august = []
v_august = []

for file in os.listdir(f'ORAS5/Data/iicevelu'):
    year = int(file[:4])
    month = int(file[5:7])
    if year >= 2011 and year <= 2018:
        if month == 4:
            u_april.append(np.loadtxt(f'ORAS5/Data/iicevelu/{file}'))
            v_april.append(np.loadtxt(f'ORAS5/Data/iicevelv/{file}'))
        if month == 8:
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
    for file in os.listdir(f"C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/{year}/04"):
        if 'ice_drift_nh_ease2-750_cdr-v1p0_24h-' in file:
            fileApril = f"C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/{year}/04/{file}"
            ds = xr.open_dataset(fileApril, decode_times = False)

            if ~lonlat_already_acquired:
                lonOSISAF = np.nan_to_num(np.array(ds['lon'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)))
                latOSISAF = np.nan_to_num(np.array(ds['lat'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)))
                """ plt.subplot(221)
                plt.imshow(lonOSISAF)
                plt.subplot(222)
                plt.imshow(latOSISAF)
                plt.subplot(223)
                plt.imshow(np.array(ds['lon'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)))
                plt.subplot(224)
                plt.imshow(np.array(ds['lat'].where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True)))
                plt.show() """

                
                lonlat_already_acquired = True

            dlat = np.array((ds['lat1'] - ds['lat']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True))[0]
            dlon = np.array((ds['lon1'] - ds['lon']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True))[0]
            v_april.append((dlat/360 * 2*np.pi) * earth_radius /(24*60*60))                                  #[m/s]
            u_april.append((dlon/360 * 2*np.pi) * np.cos(latOSISAF/360 * 2*np.pi)* earth_radius /(24*60*60)) #[m/s]
            
            ds.close
    for file in os.listdir(f"C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/{year}/08"):
        if 'ice_drift_nh_ease2-750_cdr-v1p0_24h-' in file:
            fileAugust = f"C:/Users/Augustin/Downloads/osisaf.met.no/reprocessed/ice/drift_lr/v1/merged/{year}/08/{file}"
        
            ds = xr.open_dataset(fileAugust, decode_times = False)
        
            dlat = np.array((ds['lat1'] - ds['lat']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True))[0]
            dlon = np.array((ds['lon1'] - ds['lon']).where((ds.lon > lon_min) & (ds.lon < lon_max) & (ds.lat > lat_min) & (ds.lat < lat_max), drop = True))[0]

            v_august.append((dlat/360 * 2*np.pi) * earth_radius /(24*60*60))                                  #[m/s]
            u_august.append((dlon/360 * 2*np.pi) * np.cos(latOSISAF/360 * 2*np.pi)* earth_radius /(24*60*60)) #[m/s]

            ds.close


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



u_aprilCS2[2,:]  = np.nan
v_aprilCS2[2,:]  = np.nan
u_augustCS2[2,:] = np.nan
v_augustCS2[2,:] = np.nan
u_aprilCS2[3,:]  = np.nan
v_aprilCS2[3,:]  = np.nan
u_augustCS2[3,:] = np.nan
v_augustCS2[3,:] = np.nan
u_aprilCS2[4,:]  = np.nan
v_aprilCS2[4,:]  = np.nan
u_augustCS2[4,:] = np.nan
v_augustCS2[4,:] = np.nan
u_aprilCS2[11,34]  = np.nan
v_aprilCS2[11,34]  = np.nan
u_augustCS2[11,34] = np.nan
v_augustCS2[11,34] = np.nan

##########################
######## - PLOT - ########
##########################


titles_list = ['ORAS5 - April', 'ORAS5 - August','CS2 - April', 'CS2 - August']
titlesize = 30
labelsize = 25
ticklabelsize = 21
legendsize = 15   
plt.subplots_adjust(left=0.002,
                    bottom=0.05,
                    wspace=-0.3,  
                    top=0.93)
fig = plt.figure(figsize = (26,30))
# add grid specifications
gs = fig.add_gridspec(5,2)
axs_map  = [fig.add_subplot(gs[0,0], projection = ccrs.LambertConformal(central_longitude = -7)),
            fig.add_subplot(gs[0,1], projection = ccrs.LambertConformal(central_longitude = -7)),
            fig.add_subplot(gs[1,0], projection = ccrs.LambertConformal(central_longitude = -7)),
            fig.add_subplot(gs[1,1], projection = ccrs.LambertConformal(central_longitude = -7)),
            fig.add_subplot(gs[3,0], projection = ccrs.LambertConformal(central_longitude = -7)),
            fig.add_subplot(gs[3,1], projection = ccrs.LambertConformal(central_longitude = -7)),
            fig.add_subplot(gs[4,0], projection = ccrs.LambertConformal(central_longitude = -7)),
            fig.add_subplot(gs[4,1], projection = ccrs.LambertConformal(central_longitude = -7)),]
axs_plot = [fig.add_subplot(gs[2,0]),
            fig.add_subplot(gs[2,1])]

plt.subplots_adjust(top = 0.95, bottom = 0.1, hspace = 0.3, left = 0.01,wspace=0.1,right =0.95)

##########################
###### - SIT maps - ######
##########################
    
    
for map,i in zip([ORAS5sit_march,ORAS5sit_sept,CS2sit_march,CS2sit_sept],range(4)):
    
    
    xlim = [-33, 20]
    ylim = [65, 81]
    lower_space = 3 
    rect = mpath.Path([[xlim[0], ylim[0]],
                    [xlim[1], ylim[0]],
                    [xlim[1], ylim[1]],
                    [xlim[0], ylim[1]],
                    [xlim[0], ylim[0]],
                    ]).interpolated(20)
    proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axs_map.flatten()[i]) - axs_map.flatten()[i].transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs_map.flatten()[i].set_boundary(rect_in_target)
    axs_map.flatten()[i].set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

    axs_map.flatten()[i].coastlines()
    axs_map.flatten()[i].gridlines()        

    if i <2:
        cs = axs_map.flatten()[i].pcolormesh(lon, lat, map,vmin = 0, vmax = 3,cmap = "cmo.ice", transform=ccrs.PlateCarree())
    else:
        cs = axs_map.flatten()[i].pcolormesh(lonCS2, latCS2, map,vmin = 0, vmax = 3,cmap = "cmo.ice", transform=ccrs.PlateCarree())


    axs_map.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    
cax = fig.add_axs_map([axs_map.flatten()[3].get_position().x1+0.05,axs_map.flatten()[3].get_position().y0 ,0.03,axs_map.flatten()[3].get_position().height * 2])
cb = plt.colorbar(cs, cax = cax,)
cb.set_label(label = '[m]',size = labelsize)
cb.ax.tick_params(labelsize=labelsize)

###############################
######### - SIV SIE - #########
###############################



sie_obs_monthly = sie_obs.reshape((8,12)).T
siv_obs_monthly = siv_obs.reshape((8,12)).T

sie_oras5_monthly = sie_oras5.reshape((8,12)).T
siv_oras5_monthly = siv_oras5.reshape((8,12)).T



axs_plot.flatten()[0].plot(month_name,sie_obs_monthly*1e-5,color = 'blue',alpha = 0.3)
axs_plot.flatten()[0].plot(month_name,np.nanmean(sie_obs_monthly*1e-5,axis=1),color = 'blue',linewidth = 3,label = 'CryoSat-2')
axs_plot.flatten()[0].plot(month_name,sie_oras5_monthly*1e-5,color = 'red',alpha = 0.3)
axs_plot.flatten()[0].plot(month_name,np.nanmean(sie_oras5_monthly*1e-5,axis=1),color = 'red',linewidth = 3,label = 'ORAS5')
axs_plot.flatten()[0].set_title(f'SIE',fontsize = titlesize)
axs_plot.flatten()[0].grid()
axs_plot.flatten()[0].legend(fontsize = legendsize)
axs_plot.flatten()[0].tick_params(labelsize = ticklabelsize)
axs_plot.flatten()[0].set_ylabel(r'$10^5km^2$',size = labelsize)

axs_plot.flatten()[1].plot(month_name,siv_obs_monthly,color = 'blue',alpha = 0.3)
axs_plot.flatten()[1].plot(month_name,np.nanmean(siv_obs_monthly,axis=1),color = 'blue',linewidth = 3,label = 'CryoSat-2')
axs_plot.flatten()[1].plot(month_name,siv_oras5_monthly,color = 'red',alpha = 0.3)
axs_plot.flatten()[1].plot(month_name,np.nanmean(siv_oras5_monthly,axis=1),color = 'red',linewidth = 3,label = 'ORAS5')
axs_plot.flatten()[1].set_title(f'SIV',fontsize = titlesize)
axs_plot.flatten()[1].legend(fontsize = legendsize)
axs_plot.flatten()[1].grid()
axs_plot.flatten()[1].tick_params(labelsize = ticklabelsize)
axs_plot.flatten()[1].set_ylabel(r'$km^3$',size = labelsize)

#######################
####### - SID - #######
#######################

title_list = ['ORAS5 - April','ORAS5 - August','OSI-SAF - April','OSI-SAF - August']
lonORAS5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
latORAS5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

for u,v,i in zip([u_aprilOras5, u_augustOras5, u_aprilCS2,u_augustCS2], [v_aprilOras5, v_augustOras5, v_aprilCS2,v_augustCS2], range(6,10)):
    if i in [6,7]:
        dlat = (v * 360)/(2*np.pi * earth_radius)
        dlon = (u * 360)/(2*np.pi * earth_radius * np.cos(latORAS5/360 * 2*np.pi))
    if i in [8,9]:
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
    
    if i in [6,7]:
        cs = axes.flatten()[i].pcolormesh(lonORAS5, latORAS5,   current_intensity, vmin = 0, vmax = 0.25,cmap = "cmo.speed", transform=ccrs.PlateCarree())
    if i in [8,9]:
        cs = axes.flatten()[i].pcolormesh(lonOSISAF, latOSISAF, current_intensity, vmin = 0, vmax = 0.25,cmap = "cmo.speed", transform=ccrs.PlateCarree())
    
    #Vector plot
    current_intensity_vector_plot = np.sqrt(dlon**2 + dlat**2)
    dlon /= (current_intensity_vector_plot*18)
    dlat /= (current_intensity_vector_plot*18)
    if i in [6,7]:
        axes.flatten()[i].quiver(np.array(lonORAS5),np.array(latORAS5),np.array(dlon),np.array(dlat),transform = ccrs.PlateCarree())
    if i in [8,9]:
        
        axes.flatten()[i].quiver(np.array(lonOSISAF),np.array(latOSISAF),np.array(dlon),np.array(dlat),transform = ccrs.PlateCarree())

    axes.flatten()[i].set_title(title_list[i-6], fontsize = titlesize)
    
cax = fig.add_axes([axes.flatten()[-1].get_position().x1+0.02,axes.flatten()[-1].get_position().y0 ,0.03,axes.flatten()[-1].get_position().height * 2])
cb = plt.colorbar(cs, cax = cax,)
cb.set_label(label = '[m/s]',size = labelsize)
cb.ax.tick_params(labelsize=labelsize)


plt.savefig(f'MT_plot/validation/SI.png')
plt.close()

