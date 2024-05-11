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



Vert_profile = [[78.5,-4 ]]
    
lon_oras5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat_oras5 = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

#######################################################
#### - Plotting winter and summer mean vert prof - ####
#######################################################

# Plot S,T,sigma mean over the 12 months
SProfile_list_winter     = []
TProfile_list_winter     = []
sigmaProfile_list_winter = []

SProfile_list_summer     = []
TProfile_list_summer     = []
sigmaProfile_list_summer = []
for file in os.listdir('ORAS5/Data/votemper'):
    if int(file[5:7]) in [11,12,1,2,3]: #Winter
        dsS = xr.open_dataset(f"ORAS5/Data/vosaline/{file}")
        dsT = xr.open_dataset(f"ORAS5/Data/votemper/{file}")
        depth = np.array(dsS['deptht'])

        SProfile = []
        TProfile = []
        sigmaProfile = []
        for cell_coord in Vert_profile:
            SProfile.append(np.array(dsS['vosaline'].sel(lon = cell_coord[1], lat = cell_coord[0], method = 'nearest')))
            TProfile.append(np.array(dsT['votemper'].sel(lon = cell_coord[1], lat = cell_coord[0], method = 'nearest')))
            sigmaProfile.append(sw.pden(SProfile[-1],TProfile[-1],sw.pres(depth,cell_coord[0]))-1000)
        
        SProfile_list_winter.append(np.mean(SProfile,axis = 0))
        TProfile_list_winter.append(np.mean(TProfile,axis = 0))
        sigmaProfile_list_winter.append(np.mean(sigmaProfile,axis = 0))
        dsS.close()
        dsT.close()
    if int(file[5:7]) in [5,6,7,8,9]: #Winter
        dsS = xr.open_dataset(f"ORAS5/Data/vosaline/{file}")
        dsT = xr.open_dataset(f"ORAS5/Data/votemper/{file}")
        depth = np.array(dsS['deptht'])

        SProfile = []
        TProfile = []
        sigmaProfile = []
        for cell_coord in Vert_profile:
            SProfile.append(np.array(dsS['vosaline'].sel(lon = cell_coord[1], lat = cell_coord[0], method = 'nearest')))
            TProfile.append(np.array(dsT['votemper'].sel(lon = cell_coord[1], lat = cell_coord[0], method = 'nearest')))
            sigmaProfile.append(sw.pden(SProfile[-1],TProfile[-1],sw.pres(depth,cell_coord[0]))-1000)
        SProfile_list_summer.append(np.mean(SProfile,axis = 0))
        TProfile_list_summer.append(np.mean(TProfile,axis = 0))
        sigmaProfile_list_summer.append(np.mean(sigmaProfile,axis = 0))
        dsS.close()
        dsT.close()


##############################
########## - Plot - ##########
##############################
fig, axs = plt.subplots(ncols = 2, sharey = True, figsize = (10,7))
labelsize = 18
titlesize = 18
ticklabelsize = 15
legendsize = 15 
for i in range(len(SProfile_list_summer)):
    axs[0].plot(SProfile_list_summer[i],depth, alpha = 0.2, color = 'red')
    axs[0].plot(SProfile_list_winter[i],depth, alpha = 0.2, color = 'blue')

axs[0].plot(np.mean(SProfile_list_summer,axis = 0),depth, linewidth = 2, color = 'red', label = 'Summer')
axs[0].plot(np.mean(SProfile_list_winter,axis = 0),depth, linewidth = 2, color = 'blue', label = 'Winter')
axs[0].set_xlabel('[psu]',size = labelsize)
axs[0].set_title('Salinity',fontsize = titlesize)
axs[0].set_ylim(900,0)
axs[0].set_xlim(29,35.5)
axs[0].tick_params(labelsize = ticklabelsize)
axs[0].set_ylabel('Depth [m]',size = labelsize)
axs[0].legend(fontsize = legendsize)
axs[0].grid()

for i in range(len(SProfile_list_summer)):
    axs[1].plot(TProfile_list_summer[i],depth, alpha = 0.2, color = 'red')
    axs[1].plot(TProfile_list_winter[i],depth, alpha = 0.2, color = 'blue')

axs[1].plot(np.mean(TProfile_list_summer,axis = 0),depth, linewidth = 2, color = 'red', label = 'Summer')
axs[1].plot(np.mean(TProfile_list_winter,axis = 0),depth, linewidth = 2, color = 'blue', label = 'Winter')
axs[1].set_xlabel('[°C]',size = labelsize)
axs[1].set_title('Potential temperature',fontsize = titlesize)
axs[1].tick_params(labelsize = ticklabelsize)
axs[1].set_ylim(900,0)
axs[1].set_xlim(-2,2.5)
axs[1].legend(fontsize = legendsize)
axs[1].grid()
plt.savefig('MT_plot/profile/vert_prof.png')

