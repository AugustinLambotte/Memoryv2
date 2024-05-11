import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import xarray as xr
from scipy import stats
import numpy.ma as ma
from ORAS5.useful_func_oras5 import plot as plot_oras
from ORAS5.useful_func_oras5 import show as show_ORAS
from obs_data.useful_func import plot as plot_obs
from obs_data.useful_func import show as show_obs
from datetime import date
import seawater as sw

month_name = ['Ja','Fe','Mar','Ap','May','Ju','July','Au','Se','Oc','No','De']

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

print(np.shape(siv_oras5))
print(np.shape(sie_oras5))

siv_oras5_ra = np.convolve(siv_oras5, np.ones(12)/12, mode = 'valid')
sie_oras5_ra = np.convolve(sie_oras5, np.ones(12)/12, mode = 'valid')




###############################
########### - OBS - ###########
###############################
sit_obs = []
sic_obs = []

EGC_path_obs = np.nan_to_num(np.loadtxt('obs_data/Data/geo/EGC_path.txt'))
surface_obs  = np.loadtxt('obs_data/Data/geo/surface.txt')

for file in os.listdir('obs_data/Data/gos/ke/month_mean'):
    
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2019:
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



##################################
###### - PLOT ANNUALY AV. - ######
##################################
labelsize = 25
titlesize = 25
ticklabelsize = 21
legendsize = 15    

fig,axes = plt.subplots(nrows = 1, ncols = 2,figsize = (18,7),sharex = True)

axes[0].plot(np.linspace(2011,2019,len(sie_oras5))   ,sie_oras5   ,color = 'red',alpha = 0.4)
axes[0].plot(np.linspace(2012,2019,len(sie_oras5_ra)),sie_oras5_ra,color = 'red',linewidth = 3,label = 'ORAS5')
axes[0].set_title(f'SIE',fontsize = titlesize)
axes[0].grid()
axes[0].tick_params(labelsize = ticklabelsize)
axes[0].set_ylabel(r'$km^2$',size = labelsize)
axes[0].set_xlim(2011,2020)

axes[0].plot(np.linspace(2011,2020,len(sie_obs))   ,sie_obs   ,color = 'blue',alpha = 0.4)
axes[0].plot(np.linspace(2012,2020,len(sie_obs_ra)),sie_obs_ra,color = 'blue',linewidth = 3,label = 'observations')
axes[0].legend(fontsize = legendsize)
################
axes[1].plot(np.linspace(2011,2019,len(siv_oras5))   ,siv_oras5   ,color = 'red',alpha = 0.4)
axes[1].plot(np.linspace(2012,2019,len(siv_oras5_ra)),siv_oras5_ra,color = 'red',linewidth = 3,label = 'ORAS5')
axes[1].set_title(f'SIV',fontsize = titlesize)
axes[1].grid()
axes[1].tick_params(labelsize = ticklabelsize)
axes[1].set_ylabel(r'$km^3$',size = labelsize)
axes[1].set_xlim(2011,2020)

axes[1].plot(np.linspace(2011,2020,len(siv_obs))   ,siv_obs   ,color = 'blue',alpha = 0.4)
axes[1].plot(np.linspace(2012,2020,len(siv_obs_ra)),siv_obs_ra,color = 'blue',linewidth = 3,label = 'observations')
axes[1].legend(fontsize = legendsize)


plt.savefig('validation/MT_sie_siv_averaged.png')
plt.close()

##################################
######## - Monthly plot - ########
##################################
sie_obs_monthly = sie_obs[:-12].reshape((8,12)).T
siv_obs_monthly = siv_obs[:-12].reshape((8,12)).T

sie_oras5_monthly = sie_oras5.reshape((8,12)).T
siv_oras5_monthly = siv_oras5.reshape((8,12)).T


fig,axes = plt.subplots(nrows = 1, ncols = 2,figsize = (18,7),sharex = True)


axes[0].plot(month_name,sie_obs_monthly,color = 'blue',alpha = 0.3)
axes[0].plot(month_name,np.nanmean(sie_obs_monthly,axis=1),color = 'blue',linewidth = 3,label = 'observations')
axes[0].plot(month_name,sie_oras5_monthly,color = 'red',alpha = 0.3)
axes[0].plot(month_name,np.nanmean(sie_oras5_monthly,axis=1),color = 'red',linewidth = 3,label = 'ORAS5')
axes[0].set_title(f'SIE',fontsize = titlesize)
axes[0].grid()
axes[0].legend(fontsize = legendsize)
axes[0].tick_params(labelsize = ticklabelsize)
axes[0].set_ylabel(r'$km^3$',size = labelsize)

axes[1].plot(month_name,siv_obs_monthly,color = 'blue',alpha = 0.3)
axes[1].plot(month_name,np.nanmean(siv_obs_monthly,axis=1),color = 'blue',linewidth = 3,label = 'observations')
axes[1].plot(month_name,siv_oras5_monthly,color = 'red',alpha = 0.3)
axes[1].plot(month_name,np.nanmean(siv_oras5_monthly,axis=1),color = 'red',linewidth = 3,label = 'ORAS5')
axes[1].set_title(f'SIV',fontsize = titlesize)
axes[1].legend(fontsize = legendsize)
axes[1].grid()
axes[1].tick_params(labelsize = ticklabelsize)
axes[1].set_ylabel(r'$km^3$',size = labelsize)


plt.savefig('validation/MT_sie_siv_month.png')
