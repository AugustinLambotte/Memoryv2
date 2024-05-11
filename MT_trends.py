import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import matplotlib.path as mpath
import cmocean
import xarray as xr
from scipy import stats
import numpy.ma as ma
from datetime import date
import seawater as sw
from ORAS5.useful_func_oras5 import plot
"""
    This script is used to plot on a map the value of the trend over all the time span for the EGC energy and the fresh water flux and for each pixel
"""
earth_radius = 6370*1e3
lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

"""
    return maps of the trend of sim and EGC energy
"""    
EGC_path = np.loadtxt('ORAS5/Data/geo/Current_path.txt')

sim             = []
ke              = []
u               = []
v               = []
sit             = []
sic             = []
drift_magnitude = []
eke             = []
wind_stressU    = []
wind_stressV    = []
wind_stress     = []
S0              = []
T0              = []
sigma0          = []
S30             = []
T30             = []
sigma30         = []
S60             = []
T60             = []
sigma60         = []
sohtc300        = []
sowaflup        = []

day_btwn_2011_start_and_2018_end = []
for file in os.listdir('ORAS5/Data/ke'):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        if file != '2018-12.txt':
            sim.append(np.loadtxt(f'ORAS5/Data/sim/{file}'))
        ke.append(np.loadtxt(f'ORAS5/Data/ke/{file}'))
        u.append(np.loadtxt(f'ORAS5/Data/vozocrte/{file}'))
        v.append(np.loadtxt(f'ORAS5/Data/vomecrtn/{file}'))
        sit.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        sic.append(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))
        drift_magnitude.append(np.sqrt(np.loadtxt(f'ORAS5/Data/iicevelu/{file}')**2 + np.loadtxt(f'ORAS5/Data/iicevelv/{file}')**2))
        eke.append(np.loadtxt(f'ORAS5/Data/eke/{file}'))
        wind_stress.append(np.sqrt(np.loadtxt(f'ORAS5/Data/sozotaux/{file}')**2 + np.loadtxt(f'ORAS5/Data/sometauy/{file}')**2))
        wind_stressV.append(np.loadtxt(f'ORAS5/Data/sometauy/{file}'))
        wind_stressU.append(np.loadtxt(f'ORAS5/Data/sozotaux/{file}'))
        sohtc300.append(np.loadtxt(f'ORAS5/Data/sohtc300/{file}'))
        sowaflup.append(np.loadtxt(f'ORAS5/Data/sowaflup/{file}'))
        
        Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
        Tds = xr.open_dataset(f'ORAS5/Data/votemper/{file[:-3]}nc')

        S0.append(Sds['vosaline'].sel(deptht = 0, method = 'nearest')[:])
        T0.append(Tds['votemper'].sel(deptht = 0, method = 'nearest')[:])
        sigma0.append(sw.pden(S0[-1],T0[-1],sw.pres(   0,lat))-1000)
        S30.append(Sds['vosaline'].sel(deptht = 50, method = 'nearest')[:])
        T30.append(Tds['votemper'].sel(deptht = 50, method = 'nearest')[:])
        sigma30.append(sw.pden(S30[-1],T30[-1],sw.pres( 50,lat))-1000)
        S60.append(Sds['vosaline'].sel(deptht = 260, method = 'nearest')[:])
        T60.append(Tds['votemper'].sel(deptht = 260, method = 'nearest')[:])
        sigma60.append(sw.pden(S60[-1],T60[-1],sw.pres( 260,lat))-1000)
        Sds.close()
        Tds.close()
        day_btwn_2011_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)

ke               = np.nan_to_num(np.array(ke))
u                = np.nan_to_num(np.array(u))
v                = np.nan_to_num(np.array(v))
sic              = np.nan_to_num(np.array(sic))
sim              = np.nan_to_num(np.array(sim))
sit              = np.nan_to_num(np.array(sit))
eke              = np.nan_to_num(np.array(eke))
S0               = np.nan_to_num(np.array(S0))
T0               = np.nan_to_num(np.array(T0))
sigma0           = np.nan_to_num(np.array(sigma0))
S30              = np.nan_to_num(np.array(S30))
T30              = np.nan_to_num(np.array(T30))
sigma30          = np.nan_to_num(np.array(sigma30))
S60              = np.nan_to_num(np.array(S60))
T60              = np.nan_to_num(np.array(T60))
sigma60          = np.nan_to_num(np.array(sigma60))
wind_stressU     = np.nan_to_num(np.array(wind_stressU))
wind_stressV     = np.nan_to_num(np.array(wind_stressV))
wind_stress      = np.nan_to_num(np.array(wind_stress))
drift_magnitude  = np.nan_to_num(np.array(drift_magnitude))
sohtc300         = np.nan_to_num(np.array(sohtc300))
sowaflup         = np.nan_to_num(np.array(sowaflup))

complex_wind_stress = wind_stressU +1j*wind_stressV
angle_wind_stress   = np.angle(complex_wind_stress,deg = True)

#ke[:,EGC_path != 1]               = np.nan
#u[:,EGC_path != 1]               = np.nan
#v[:,EGC_path != 1]               = np.nan
#sic[:,EGC_path != 1]              = np.nan
#sim[:,EGC_path != 1]              = np.nan
#sit[:,EGC_path != 1]              = np.nan
#eke[:,EGC_path != 1]              = np.nan
#wind_stress[:,EGC_path != 1]      = np.nan
drift_magnitude[:,EGC_path != 1]  = np.nan
#S[:,EGC_path != 1]                = np.nan
#T[:,EGC_path != 1]                = np.nan

trend_sim               = np.zeros((np.shape(sim)[1:]))
trend_ke                = np.zeros((np.shape(sim)[1:]))
trend_u                 = np.zeros((np.shape(sim)[1:]))
trend_v                 = np.zeros((np.shape(sim)[1:]))
trend_sit               = np.zeros((np.shape(sim)[1:]))
trend_sic               = np.zeros((np.shape(sim)[1:]))
trend_drift             = np.zeros((np.shape(sim)[1:]))
trend_eke               = np.zeros((np.shape(sim)[1:]))
trend_wind_stressU      = np.zeros((np.shape(sim)[1:]))
trend_wind_stressV      = np.zeros((np.shape(sim)[1:]))
trend_wind_stress       = np.zeros((np.shape(sim)[1:]))
trend_S0                = np.zeros((np.shape(sim)[1:]))
trend_T0                = np.zeros((np.shape(sim)[1:]))
trend_sigma0            = np.zeros((np.shape(sim)[1:]))
trend_S30               = np.zeros((np.shape(sim)[1:]))
trend_T30               = np.zeros((np.shape(sim)[1:]))
trend_sigma30           = np.zeros((np.shape(sim)[1:]))
trend_S60               = np.zeros((np.shape(sim)[1:]))
trend_T60               = np.zeros((np.shape(sim)[1:]))
trend_sigma60           = np.zeros((np.shape(sim)[1:]))
trend_angle_wind_stress = np.zeros((np.shape(sim)[1:]))
trend_sohtc300          = np.zeros((np.shape(sim)[1:]))
trend_sowaflup          = np.zeros((np.shape(sim)[1:]))

p_sim                    = np.zeros((np.shape(sim)[1:]))
p_ke                     = np.zeros((np.shape(sim)[1:]))
p_u                      = np.zeros((np.shape(sim)[1:]))
p_v                      = np.zeros((np.shape(sim)[1:]))
p_sit                    = np.zeros((np.shape(sim)[1:]))
p_sic                    = np.zeros((np.shape(sim)[1:]))
p_drift                  = np.zeros((np.shape(sim)[1:]))
p_eke                    = np.zeros((np.shape(sim)[1:]))
p_wind_stressU           = np.zeros((np.shape(sim)[1:]))
p_wind_stressV           = np.zeros((np.shape(sim)[1:]))
p_wind_stress            = np.zeros((np.shape(sim)[1:]))
p_S0                     = np.zeros((np.shape(sim)[1:]))
p_T0                     = np.zeros((np.shape(sim)[1:]))
p_sigma0                 = np.zeros((np.shape(sim)[1:]))
p_S30                    = np.zeros((np.shape(sim)[1:]))
p_T30                    = np.zeros((np.shape(sim)[1:]))
p_sigma30                = np.zeros((np.shape(sim)[1:]))
p_S60                    = np.zeros((np.shape(sim)[1:]))
p_T60                    = np.zeros((np.shape(sim)[1:]))
p_sigma60                = np.zeros((np.shape(sim)[1:]))
p_angle_wind_stress      = np.zeros((np.shape(sim)[1:]))
p_sohtc300               = np.zeros((np.shape(sim)[1:]))
p_sowaflup              = np.zeros((np.shape(sim)[1:]))

for l in range(np.shape(trend_sim)[0]):
    for c in range(np.shape(trend_sim)[1]):
        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:-1],np.convolve(sim[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sim[l,c]             = stats.linregress(day_btwn_2011_start_and_2018_end[11:-1],np.convolve(sim[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sim[l,c]         = slope * (365)

        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(ke[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_ke[l,c]              = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(ke[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_ke[l,c]          = slope * (365)

        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(u[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_u[l,c]               = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(u[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_u[l,c]           = slope * (365)

        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(v[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_v[l,c]               = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(v[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_v[l,c]           = slope * (365)
        
        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sit[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sit[l,c]             = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sit[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sit[l,c]         = slope * (365)
        
        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sic[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sic[l,c]             = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sic[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sic[l,c]         = slope * (365)

        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(drift_magnitude[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_drift[l,c]           = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(drift_magnitude[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_drift[l,c]       = slope * (365)

        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(eke[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_eke[l,c]             = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(eke[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_eke[l,c]         = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(wind_stressU[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_wind_stressU[l,c]          = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(wind_stressU[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_wind_stressU[l,c]      = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(wind_stressV[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_wind_stressV[l,c]          = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(wind_stressV[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_wind_stressV[l,c]      = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(wind_stress[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_wind_stress[l,c]           = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(wind_stress[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_wind_stress[l,c]       = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S0[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_S0[l,c]                    = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S0[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_S0[l,c]                = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(T0[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_T0[l,c]                    = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(T0[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_T0[l,c]                = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S30[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_S30[l,c]                   = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S30[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_S30[l,c]               = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(T30[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_T30[l,c]                   = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(T30[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_T30[l,c]               = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S60[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_S60[l,c]                   = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S60[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_S60[l,c]               = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(T60[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_T60[l,c]                   = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(T60[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_T60[l,c]               = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sigma0[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sigma0[l,c]                = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sigma0[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sigma0[l,c]            = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sigma30[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sigma30[l,c]               = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sigma30[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sigma30[l,c]           = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sigma60[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sigma60[l,c]               = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sigma60[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sigma60[l,c]           = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(angle_wind_stress[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_angle_wind_stress[l,c]     = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(angle_wind_stress[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_angle_wind_stress[l,c] = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sohtc300[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sohtc300[l,c]              = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sohtc300[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sohtc300[l,c]          = slope * (365)

        slope                        = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sowaflup[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sowaflup[l,c]              = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sowaflup[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sowaflup[l,c]          = slope * (365)



trend_ke[trend_ke == 0]                               = np.nan
trend_u[trend_u == 0]                                 = np.nan
trend_v[trend_v == 0]                                 = np.nan
trend_sim[trend_sim == 0]                             = np.nan
trend_sit[trend_sit == 0]                             = np.nan
trend_sic[trend_sic == 0]                             = np.nan
trend_drift[trend_drift == 0]                         = np.nan
trend_eke[trend_eke == 0]                             = np.nan
trend_wind_stressU[trend_wind_stressU == 0]           = np.nan
trend_wind_stressV[trend_wind_stressV == 0]           = np.nan
trend_wind_stress[trend_wind_stress == 0]             = np.nan
trend_S0[trend_S0 == 0]                               = np.nan
trend_T0[trend_T0 == 0]                               = np.nan
trend_sigma0[trend_sigma0 == 0]                       = np.nan
trend_S30[trend_S30 == 0]                             = np.nan
trend_T30[trend_T30 == 0]                             = np.nan
trend_sigma30[trend_sigma30 == 0]                     = np.nan
trend_S60[trend_S60 == 0]                             = np.nan
trend_T60[trend_T60 == 0]                             = np.nan
trend_sigma60[trend_sigma60 == 0]                     = np.nan
trend_angle_wind_stress[trend_angle_wind_stress == 0] = np.nan
trend_sohtc300[trend_sohtc300 == 0]                   = np.nan
trend_sowaflup[trend_sowaflup == 0]                   = np.nan


p_sim[p_sim == 1]                             = np.nan
p_ke[p_ke == 1]                               = np.nan
p_v[p_v == 1]                                 = np.nan
p_u[p_u == 1]                                 = np.nan
p_sit[p_sit == 1]                             = np.nan
p_sic[p_sic == 1]                             = np.nan
p_drift[p_drift == 1]                         = np.nan
p_eke[p_eke == 1]                             = np.nan
p_wind_stressU[p_wind_stressU == 1]           = np.nan
p_wind_stressV[p_wind_stressV == 1]           = np.nan
p_wind_stress[p_wind_stress == 1]             = np.nan
p_S0[p_S0 == 1]                               = np.nan
p_T0[p_T0 == 1]                               = np.nan
p_sigma0[p_sigma0 == 1]                       = np.nan
p_S30[p_S30 == 1]                             = np.nan
p_T30[p_T30 == 1]                             = np.nan
p_sigma30[p_sigma30 == 1]                     = np.nan
p_S60[p_S60 == 1]                             = np.nan
p_T60[p_T60 == 1]                             = np.nan
p_sigma60[p_sigma60 == 1]                     = np.nan
p_angle_wind_stress[p_angle_wind_stress == 1] = np.nan
p_sohtc300[p_sohtc300 == 1]                   = np.nan
p_sowaflup[p_sowaflup == 1]                   = np.nan



####################################
#### - Plotting SIT SIC KE EKE- ####
####################################

""" titles_list = ['SIT', 'SIC']#,'KE', 'EKE']
label_list = [r'$[m*y^{-1}]$', r'$[\%*y^{-1}]$']#, r'$[mJ*kg^{-1}*y^{-1}]$',r'$[mJ*kg^{-1}*y^{-1}]$']
v_list = [0.15,5,1,0.5]
titlesize = 30
labelsize = 23
fig,axes = plt.subplots(ncols = 2, figsize = (16,8), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -11)})
plt.subplots_adjust(left=0.001,
                    bottom=0.05,
                    top=0.93, 
                    hspace=0.1,
                    wspace = 0.25)



for trend_map,p_map,i in zip([trend_sit,trend_sic*1e2],[p_sit,p_sic],range(2)):
    
    xlim = [-33, 15]
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

    cs = axes.flatten()[i].pcolormesh(lon, lat, trend_map,vmin = -v_list[i], vmax = v_list[i],cmap = 'cmo.balance', transform=ccrs.PlateCarree())
    contour = axes.flatten()[i].contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.01,axes.flatten()[i].get_position().y0 - 0.02,0.015,axes.flatten()[i].get_position().height])
    cb = plt.colorbar(cs, cax = cax)
    cb.ax.tick_params(labelsize=labelsize)
    cb.set_label(label = label_list[i],size = labelsize)

    coords_p_signi = []
    for k in range(len(p_map)):
        for l in range(len(p_map[0])):
            if p_map[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axes.flatten()[i].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())

    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    

plt.savefig(f'MT_plot/Trend_sit_sic.png')
plt.close() """
"""
#################################
#### - Plotting KE EKE u v - ####
#################################

titles_list = ['KE', 'EKE', 'U', 'V']
label_list = [ r'$[mJ*kg^{-1}*y^{-1}]$',r'$[mJ*kg^{-1}*y^{-1}]$',r'$[cm*s^{-1}*y^{-1}]$',r'$[cm*s^{-1}*y^{-1}]$']
v_list = [1,0.5,1,1]

titlesize = 25
labelsize = 20
fig,axes = plt.subplots(ncols = 2,nrows = 2, figsize = (16,11), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -6)})
plt.subplots_adjust(left=0.01,
                    bottom=0.05,
                    top=0.93, 
                    hspace=0.1,
                    wspace = 0)


for trend_map,p_map,i in zip([trend_ke*1e3,trend_eke*1e3,trend_u*1e2,trend_v*1e2],[p_ke,p_eke,p_u,p_v],range(4)):
    
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

    cs = axes.flatten()[i].pcolormesh(lon, lat, trend_map,vmin = -v_list[i], vmax = v_list[i],cmap = 'cmo.balance', transform=ccrs.PlateCarree())
    contour = axes.flatten()[i].contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.01,axes.flatten()[i].get_position().y0 - 0.02,0.015,axes.flatten()[i].get_position().height])
    cb = plt.colorbar(cs, cax = cax)
    cb.ax.tick_params(labelsize=labelsize)
    cb.set_label(label = label_list[i],size = labelsize)

    coords_p_signi = []
    for k in range(len(p_map)):
        for l in range(len(p_map[0])):
            if p_map[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axes.flatten()[i].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())

    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    

plt.savefig(f'MT_plot/Trend_ke_eke_u_v.png')
plt.close()

#############################
#### - Plotting KE EKE - ####
#############################

titles_list = ['KE', 'EKE']
label_list = [ r'$[mJ*kg^{-1}*y^{-1}]$',r'$[mJ*kg^{-1}*y^{-1}]$']
v_list = [1,1]

titlesize = 25
labelsize = 20
fig,axes = plt.subplots(ncols = 2, figsize = (18,8), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -6)})
plt.subplots_adjust(left=0.01,
                    bottom=0.05,
                    top=0.93, 
                    hspace=0.1,
                    wspace = 0)


for trend_map,p_map,i in zip([trend_ke*1e3,trend_eke*1e3],[p_ke,p_eke],range(2)):
    
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

    cs = axes.flatten()[i].pcolormesh(lon, lat, trend_map,vmin = -v_list[i], vmax = v_list[i],cmap = 'cmo.balance', transform=ccrs.PlateCarree())
    contour = axes.flatten()[i].contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    coords_p_signi = []
    for k in range(len(p_map)):
        for l in range(len(p_map[0])):
            if p_map[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axes.flatten()[i].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())

    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    
cax = fig.add_axes([axes.flatten()[-1].get_position().x1+0.01,axes.flatten()[-1].get_position().y0 - 0.02,0.015,axes.flatten()[-1].get_position().height])
cb = plt.colorbar(cs, cax = cax,ticks = [-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
cb.ax.tick_params(labelsize=labelsize)
cb.set_label(label = label_list[i],size = labelsize)
plt.savefig(f'MT_plot/Trend_ke_eke.png')
plt.close()
"""
#####################################
#### - Plotting S, T and sigma - ####
#####################################

titles_list = ['Salinity\nSurface', 'Potential temperature\nSurface','Density anomaly'+'\nSurface','50m depth', '50m depth','50m depth']

label_list = [r'$[psu\cdot y^{-1}]$',r'$[°C\cdot y^{-1}]$',r'$[kg\cdot m^{-3}\cdot y^{-1}]$']
ticks_list = [[-0.1,-0.05,0,0.05,0.1],[-0.4,-0.2,0,0.2,0.4],[-0.1,-0.05,0,0.05,0.1]]
v_list = [0.1,0.5,0.1,0.1,0.5,0.1]
titlesize = 40
labelsize = 35
fig,axes = plt.subplots(ncols = 3, nrows = 2, figsize = (30,16), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -8)})
plt.subplots_adjust(left=0.01,
                    bottom=0.05, 
                    right=0.9, 
                    top=0.91,
                    wspace = 0.1,
                    hspace = 0.15)

mean_ke = np.nanmean(ke,axis = 0)
for trend_map,p_map,i in zip([trend_S0,trend_T0,trend_sigma0,trend_S30,trend_T30,trend_sigma30],[p_S0,p_T0,p_sigma0,p_S30,p_T30,p_sigma30],range(6)):
    trend_map[trend_map<-1] = np.nan
    if i == 2 or i ==5:
        trend_map[trend_map<-0.1] = np.nan
    xlim = [-33, 18]
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

    cs = axes.flatten()[i].pcolormesh(lon, lat, trend_map,vmin = -v_list[i], vmax = v_list[i],cmap = 'cmo.balance', transform=ccrs.PlateCarree())
    contour = axes.flatten()[i].contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())
    if i > 2:
        cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.01,axes.flatten()[i].get_position().y0 + 0.03,0.015,axes.flatten()[i].get_position().height*2])
        cb = plt.colorbar(cs, cax = cax, ticks = ticks_list[i-3])
        cb.ax.tick_params(labelsize=labelsize)
        cb.set_label(label = label_list[i-3],size = labelsize)

    coords_p_signi = []
    for k in range(len(p_map)):
        for l in range(len(p_map[0])):
            if p_map[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axes.flatten()[i].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())
    
    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    

plt.savefig(f'MT_plot/Trend_S_T_depth.png')
plt.close()
"""
###########################
##### - Wind Stress - #####
###########################

titles_list = ['Magnitude','Angle']

label_list = [r'$[mN\cdot m^{-2}\cdot y^{-1}]$',r'[$\circ \cdot y^{-1}$]']
v_list = [10,15]
titlesize = 30
labelsize = 22
fig,axes = plt.subplots( nrows = 1, ncols = 2,figsize = (16,7), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -9)})

plt.subplots_adjust(left=0.01,
                    bottom=0.05,  
                    top=0.93,
                    wspace = 0.17)

for wind_stress_comp,p_wind_stress_comp,i in zip([trend_wind_stress*1e+3,trend_angle_wind_stress],[p_wind_stress,p_angle_wind_stress],range(2)):


    lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
    lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

    xlim = [-33, 17]
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

    
    cs = axes.flatten()[i].pcolormesh(lon, lat,wind_stress_comp,vmin = -v_list[i], vmax = v_list[i],cmap = "cmo.balance", transform=ccrs.PlateCarree())
    contour = axes.flatten()[i].contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    
    cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.005,axes.flatten()[i].get_position().y0,0.015,axes.flatten()[i].get_position().height])
    if i == 0:
        cb = plt.colorbar(cs, cax = cax,ticks = [-10,-5,0,5,10])
    else:
        cb = plt.colorbar(cs, cax = cax)

    cb.ax.tick_params(labelsize=labelsize)
    cb.set_label(label = label_list[i],size = labelsize) 

    coords_p_signi = []
    for k in range(len(p_u)):
        for l in range(len(p_u[0])):
            if p_wind_stress_comp[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axes.flatten()[i].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())
    
    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    

plt.savefig(f'MT_plot/Trend_wind_stress.png')
plt.close() 

""" ############################
####### - sohtc300 - #######
############################

titles_list = ['sohtc300']

label_list = [r'$[GJ*m^{-2}*y^{-1}]$']
v_list = [0.5]
titlesize = 25
labelsize = 17
fig,axes = plt.subplots(figsize = (10,7), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -9)})

plt.subplots_adjust(left=0.01,
                    bottom=0.05,  
                    top=0.93,
                    wspace = 0.13)

for map,p,i in zip([trend_sohtc300*1e-9],[p_sohtc300],range(1)):


    lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
    lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

    xlim = [-33, 17]
    ylim = [65, 81]
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

    
    cs = axes.pcolormesh(lon, lat,map,vmin = -v_list[i], vmax = v_list[i],cmap = "cmo.balance", transform=ccrs.PlateCarree())
    contour = axes.contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    
    cax = fig.add_axes([axes.get_position().x1+0.01,axes.get_position().y0 - 0.02,0.015,axes.get_position().height])
    cb = plt.colorbar(cs, cax = cax)
    cb.ax.tick_params(labelsize=15)
    cb.set_label(label = label_list[i],size = 15) 

    coords_p_signi = []
    for k in range(len(p_u)):
        for l in range(len(p_u[0])):
            if p[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axes.scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())
    
    axes.set_title(titles_list[i],fontsize = titlesize)
    

plt.savefig(f'MT_plot/Trend_sohtc300.png')
plt.close()

############################
####### - sowaflup - #######
############################

titles_list = ['sowaflup']

label_list = [r'$[GJ*m^{-2}*y^{-1}]$']
v_list = [0.000000000000003]
titlesize = 25
labelsize = 17
fig,axes = plt.subplots(figsize = (10,7), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -9)})

plt.subplots_adjust(left=0.01,
                    bottom=0.05,  
                    top=0.93,
                    wspace = 0.13)

for map,p,i in zip([trend_sowaflup*1e-9],[p_sowaflup],range(1)):


    lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
    lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

    xlim = [-33, 17]
    ylim = [65, 81]
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

    
    cs = axes.pcolormesh(lon, lat,map,vmin = -v_list[i], vmax = v_list[i],cmap = "cmo.balance", transform=ccrs.PlateCarree())
    contour = axes.contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    
    cax = fig.add_axes([axes.get_position().x1+0.01,axes.get_position().y0 - 0.02,0.015,axes.get_position().height])
    cb = plt.colorbar(cs, cax = cax)
    cb.ax.tick_params(labelsize=15)
    cb.set_label(label = label_list[i],size = 15) 

    coords_p_signi = []
    for k in range(len(p_u)):
        for l in range(len(p_u[0])):
            if p[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axes.scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())
    
    axes.set_title(titles_list[i],fontsize = titlesize)
    

plt.savefig(f'MT_plot/Trend_sowaflup.png')
plt.close()

####################################
####### - Plot pré-défense - #######
####################################

titles_list = ['Sea ice thickness', 'Sea ice concentration']
label_list = [r'$[m\cdot y^{-1}]$', r'$[\%\cdot y^{-1}]$']
v_list = [0.15,5,1,0.5]
titlesize = 30
labelsize = 25
fig,axes = plt.subplots(ncols = 2, figsize = (17,7), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -7)})
plt.subplots_adjust(left=0.001,
                    bottom=0.05,
                    top=0.93, 
                    wspace = 0.4)



for trend_map,p_map,i in zip([trend_sit,trend_sic*1e2],[p_sit,p_sic],range(2)):
    
    xlim = [-30, 17.5]
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

    cs = axes.flatten()[i].pcolormesh(lon, lat, trend_map,vmin = -v_list[i], vmax = v_list[i],cmap = 'cmo.balance', transform=ccrs.PlateCarree())
    contour = axes.flatten()[i].contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.01,axes.flatten()[i].get_position().y0 - 0.02,0.015,axes.flatten()[i].get_position().height])
    cb = plt.colorbar(cs, cax = cax)
    cb.ax.tick_params(labelsize=labelsize)
    cb.set_label(label = label_list[i],size = labelsize)

    coords_p_signi = []
    for k in range(len(p_map)):
        for l in range(len(p_map[0])):
            if p_map[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axes.flatten()[i].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())

    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    

plt.savefig(f'MT_plot/pre_def/Trend_sit_sic.png')
plt.close()