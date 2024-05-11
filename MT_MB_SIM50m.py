import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import cartopy.crs as ccrs
import gsw
import os
import matplotlib.path as mpath
from ORAS5.useful_func_oras5 import plot
from scipy import stats
from datetime import date
import seawater as sw



lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")

S0 = []

day_btwn_2011_start_and_2018_end = []
for file in os.listdir('ORAS5/Data/ke'):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        
        Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
        Tds = xr.open_dataset(f'ORAS5/Data/votemper/{file[:-3]}nc')
        day_btwn_2011_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)
        S0.append(Sds['vosaline'].sel(deptht = 50, method = 'nearest')[:])
S0       = np.nan_to_num(np.array(S0))
trend_S0 = np.zeros((np.shape(S0)[1:]))
p_S0     = np.zeros((np.shape(S0)[1:]))

for l in range(np.shape(trend_S0)[0]):
    for c in range(np.shape(trend_S0)[1]):
        slope         = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S0[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_S0[l,c]     = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S0[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_S0[l,c] = slope * (365)


trend_S0[trend_S0 == 0] = np.nan
p_S0[p_S0 == 1]         = np.nan
trend_S0[trend_S0<-1] = np.nan

surface   = np.loadtxt('ORAS5/Data/geo/Surface.txt')
X_dist    = np.loadtxt('ORAS5/Data/geo/X_dist.txt')
Y_dist    = np.loadtxt('ORAS5/Data/geo/Y_dist.txt')
earth_radius = 6370 * 1e3
EGC_path = np.loadtxt('ORAS5/Data/geo/Current_path.txt')

def MB_G1G2(file):
    iicevelv = np.loadtxt(f'ORAS5/Data/iicevelv/{file}')
    iicevelu = np.loadtxt(f'ORAS5/Data/iicevelu/{file}')
    sit      = np.loadtxt(f'ORAS5/Data/iicethic/{file}')
    sic      = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')

    ##### - GATE 1 - #####

    fromNorth = 0
    for lon_ in range(-12,18):
        sit_interp = (np.nan_to_num(float(sit[     (lat == 82.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == 81.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == 82.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == 81.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == 82.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == 81.5) & (lon == lon_)])))/2

        fromNorth -= sit_interp * sic_interp * v_interp * float(X_dist[(lat == 82.5) & (lon == lon_)])
    
    fromEast = 0
    for lat_ in [81.5,80.5]:
        sit_interp = (np.nan_to_num(float(sit[     (lat == lat_) & (lon == 17)])) + np.nan_to_num(float(sit[     (lat == lat_) & (lon == 18)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == lat_) & (lon == 17)])) + np.nan_to_num(float(sic[     (lat == lat_) & (lon == 18)])))/2
        u_interp   = (np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == 17)])) + np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == 18)])))/2

        fromEast -= sit_interp * sic_interp * u_interp * float(Y_dist[(lat == lat_) & (lon == 17)])
    
    
    MBG1 = (fromNorth + fromEast)*60*60*24*30 #m^3/month

    ##### - GATE 2 - #####

    fromNorth = 0
    for lon_ in range(-19,12):
        sit_interp = (np.nan_to_num(float(sit[     (lat == 78.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == 79.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == 78.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == 79.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == 78.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == 79.5) & (lon == lon_)])))/2

        fromNorth -= sit_interp * sic_interp * v_interp * float(X_dist[(lat == 79.5) & (lon == lon_)])
   

    MBG2 = (fromNorth)*60*60*24*30 #m^3/month
    
    

    return MBG1,MBG2

def SIM_A1A2(file):
    iicevelv = np.loadtxt(f'ORAS5/Data/iicevelv/{file}')
    iicevelu = np.loadtxt(f'ORAS5/Data/iicevelu/{file}')
    sit      = np.loadtxt(f'ORAS5/Data/iicethic/{file}')
    sic      = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')

    
    ########################    
    ###### - ZONE 1 - ######
    ########################
    
    ##### - GATE 1 - #####

    fromNorth = 0
    for lon_ in range(-12,18):
        sit_interp = (np.nan_to_num(float(sit[     (lat == 82.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == 81.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == 82.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == 81.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == 82.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == 81.5) & (lon == lon_)])))/2

        fromNorth += sit_interp * sic_interp * v_interp * float(X_dist[(lat == 82.5) & (lon == lon_)])
    
    fromEast = 0
    for lat_ in [81.5,80.5]:
        sit_interp = (np.nan_to_num(float(sit[     (lat == lat_) & (lon == 17)])) + np.nan_to_num(float(sit[     (lat == lat_) & (lon == 18)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == lat_) & (lon == 17)])) + np.nan_to_num(float(sic[     (lat == lat_) & (lon == 18)])))/2
        u_interp   = (np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == 17)])) + np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == 18)])))/2

        fromEast += sit_interp * sic_interp * u_interp * float(Y_dist[(lat == lat_) & (lon == 17)])
    
    
    MBG1 = (fromNorth + fromEast)*60*60*24*30 #m^3/month

    ##### - GATE 2 - #####

    fromNorth = 0
    for lon_ in range(-19,12):
        sit_interp = (np.nan_to_num(float(sit[     (lat == 78.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == 79.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == 78.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == 79.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == 78.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == 79.5) & (lon == lon_)])))/2

        fromNorth += sit_interp * sic_interp * v_interp * float(X_dist[(lat == 79.5) & (lon == lon_)])
   

    MBG2 = (fromNorth)*60*60*24*30 #m^3/month

    MB = MBG2 - MBG1 #Positive when SI coming in A1

    year = int(file[:4])
    month = int(file[5:7])
    if year == 2011 and month == 1:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/{file}')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/2011-02.txt')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/2011-02.txt')
    elif year == 2018 and month == 12:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/2018-11.txt')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/{file}')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/2018-11.txt')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')
    elif month == 1:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/{year-1}-12.txt')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/{year  }-02.txt')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/{year-1}-12.txt')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/{year  }-02.txt')
    elif month == 12:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/{year  }-11.txt')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/{year+1}-01.txt')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/{year  }-11.txt')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/{year+1}-01.txt')
    else:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/{year}-{month-1:02d}.txt')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/{year}-{month+1:02d}.txt')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/{year}-{month-1:02d}.txt')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/{year}-{month+1:02d}.txt')
    
    previousSIT = np.where((lat <= 81.5) & (lat >= 79.5) & (lon < 17) & (lon > -20),previousSIT,0)
    previousSIC = np.where((lat <= 81.5) & (lat >= 79.5) & (lon < 17) & (lon > -20),previousSIC,0)
    nextSIT     = np.where((lat <= 81.5) & (lat >= 79.5) & (lon < 17) & (lon > -20),nextSIT,    0)
    nextSIC     = np.where((lat <= 81.5) & (lat >= 79.5) & (lon < 17) & (lon > -20),nextSIC,    0)

    startMonthSIV = (np.nansum(previousSIT*previousSIC*surface) + np.nansum(np.where((lat <= 81.5) & (lat >= 79.5) & (lon < 17) & (lon > -20),sit*sic*surface,0)))/2
    endMonthSIV   = (np.nansum(nextSIT*nextSIC*surface)         + np.nansum(np.where((lat <= 81.5) & (lat >= 79.5) & (lon < 17) & (lon > -20),sit*sic*surface,0)))/2

    bilanSIV = endMonthSIV - startMonthSIV #Positive when SIV increase over area

    SIMA1 = MB - bilanSIV

    ########################    
    ###### - ZONE 2 - ######
    ########################

    fromNorth = 0
    for lon_ in range(-13,-2):
        sit_interp = (np.nan_to_num(float(sit[     (lat == 78.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == 79.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == 78.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == 79.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == 78.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == 79.5) & (lon == lon_)])))/2

        fromNorth -= sit_interp * sic_interp * v_interp * float(X_dist[(lat == 79.5) & (lon == lon_)])
   
    fromSouth = 0
    for lon_ in range(-13,-2):
        sit_interp = (np.nan_to_num(float(sit[     (lat == 77.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == 76.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == 77.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == 76.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == 77.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == 76.5) & (lon == lon_)])))/2

        fromSouth += sit_interp * sic_interp * v_interp * float(X_dist[(lat == 76.5) & (lon == lon_)])

    fromWest = 0
    for lat_ in [77.5,78.5]:
        sit_interp = (np.nan_to_num(float(sit[     (lat == lat_) & (lon == -14)])) + np.nan_to_num(float(sit[     (lat == lat_) & (lon == -13)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == lat_) & (lon == -14)])) + np.nan_to_num(float(sic[     (lat == lat_) & (lon == -13)])))/2
        v_interp   = (np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == -14)])) + np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == -13)])))/2

        fromWest += sit_interp * sic_interp * v_interp * float(Y_dist[(lat == lat_) & (lon == lon_)])

    fromEast = 0
    for lat_ in [77.5,78.5]:
        sit_interp = (np.nan_to_num(float(sit[     (lat == lat_) & (lon == -3)])) + np.nan_to_num(float(sit[     (lat == lat_) & (lon == -2)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == lat_) & (lon == -3)])) + np.nan_to_num(float(sic[     (lat == lat_) & (lon == -2)])))/2
        v_interp   = (np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == -3)])) + np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == -2)])))/2

        fromEast -= sit_interp * sic_interp * v_interp * float(Y_dist[(lat == lat_) & (lon == lon_)])
    
    
    MB = (fromNorth + fromSouth + fromEast + fromWest)*60*60*24*30 #m^3/month

    year = int(file[:4])
    month = int(file[5:7])
    if year == 2011 and month == 1:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/{file}')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/2011-02.txt')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/2011-02.txt')
    elif year == 2018 and month == 12:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/2018-11.txt')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/{file}')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/2018-11.txt')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')
    elif month == 1:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/{year-1}-12.txt')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/{year  }-02.txt')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/{year-1}-12.txt')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/{year  }-02.txt')
    elif month == 12:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/{year  }-11.txt')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/{year+1}-01.txt')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/{year  }-11.txt')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/{year+1}-01.txt')
    else:
        previousSIT = np.loadtxt(f'ORAS5/Data/iicethic/{year}-{month-1:02d}.txt')
        nextSIT     = np.loadtxt(f'ORAS5/Data/iicethic/{year}-{month+1:02d}.txt')

        previousSIC = np.loadtxt(f'ORAS5/Data/ileadfra/{year}-{month-1:02d}.txt')
        nextSIC     = np.loadtxt(f'ORAS5/Data/ileadfra/{year}-{month+1:02d}.txt')
    
    previousSIT = np.where((lat >= 77.5) & (lat <= 78.5) & (lon <= -3) & (lon >= -13),previousSIT,0)
    previousSIC = np.where((lat >= 77.5) & (lat <= 78.5) & (lon <= -3) & (lon >= -13),previousSIC,0)
    nextSIT     = np.where((lat >= 77.5) & (lat <= 78.5) & (lon <= -3) & (lon >= -13),nextSIT,    0)
    nextSIC     = np.where((lat >= 77.5) & (lat <= 78.5) & (lon <= -3) & (lon >= -13),nextSIC,    0)


    startMonthSIV = (np.nansum(previousSIT*previousSIC*surface) + np.nansum(np.where((lat >= 77.5) & (lat <= 78.5) & (lon <= -3) & (lon >= -13),sit*sic*surface,0)))/2
    endMonthSIV   = (np.nansum(nextSIT*nextSIC*surface)         + np.nansum(np.where((lat >= 77.5) & (lat <= 78.5) & (lon <= -3) & (lon >= -13),sit*sic*surface,0)))/2

    bilanSIV = endMonthSIV - startMonthSIV #Positive when SIV increase over area

    SIMA2 = MB - bilanSIV

    return SIMA1, SIMA2

def S_A2(file,depth):
    Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
    Smap = Sds['vosaline'].sel(deptht = depth, method = 'nearest')[:]
    Sds.close()
    S = np.nanmean(np.where((lat >= 77.5) & (lat <= 78.5) & (lon < -2) & (lon >= -15),Smap,np.nan))
    return S

gate1 = np.concatenate(([[82,i] for i in range(-12,18)],[[81,17],[80,17]]))
gate2 = [[79,i] for i in range(-19,12)]
A2 = np.concatenate(([[77,i] for i in range(-13,-2)],[[79,-3]],[[79,i] for i in range(-3,-14,-1)],[[77,-13]]))
labelsize = 20
titlesize = 27
ticklabelsize = 20
legendsize = 20
    
fig = plt.figure(figsize = (27,16))
# add grid specifications
gs = fig.add_gridspec(3,2)
axs_map = fig.add_subplot(gs[:,0], projection = ccrs.LambertConformal(central_longitude = -3))
axs_plot = []
for i in range(3):
    axs_plot.append(fig.add_subplot(gs[i,1]))
plt.subplots_adjust(top = 0.95, bottom = 0.1, hspace = 0.3, left = 0.01,wspace=0.1,right =0.95)

##############################
#### - plotting the map - ####
##############################
xlim = [-25, 20]
ylim = [72, 83]
lower_space = 3 
rect = mpath.Path([[xlim[0], ylim[0]],
                [xlim[1], ylim[0]],
                [xlim[1], ylim[1]],
                [xlim[0], ylim[1]],
                [xlim[0], ylim[0]],
                ]).interpolated(20)
proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axs_map) - axs_map.transData
rect_in_target = proj_to_data.transform_path(rect)
axs_map.set_boundary(rect_in_target)
axs_map.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

axs_map.coastlines()
axs_map.gridlines()        

cs = axs_map.pcolormesh(lon, lat, trend_S0,vmin = -0.1, vmax = 0.1,cmap = "cmo.balance", transform=ccrs.PlateCarree())
cax = fig.add_axes([axs_map.get_position().x1 - 0.4,axs_map.get_position().y0 - 0.02,axs_map.get_position().width*0.8,axs_map.get_position().height*0.04])
cb = plt.colorbar(cs, cax = cax, orientation = 'horizontal', ticks = [-0.1,-0.05,0,0.05,0.1])
cb.ax.tick_params(labelsize=labelsize)
cb.set_label(label = r'$[psu*y^{-1}]$',size = labelsize)
### - Transect - ###

axs_map.plot(np.array(gate1)[:,1], np.array(gate1)[:,0], linewidth = 4,color = 'black',transform=ccrs.PlateCarree())
axs_map.plot(np.array(gate2)[:,1], np.array(gate2)[:,0], linewidth = 4,color = 'black',transform=ccrs.PlateCarree())
axs_map.plot(np.array(A2)[:,1], np.array(A2)[:,0], linewidth = 4,color = 'red', linestyle = 'dashed',transform=ccrs.PlateCarree())

axs_map.text(np.array(gate1)[0,1]-3,np.array(gate1)[0,0],'G1',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
axs_map.text(np.array(gate1)[0,1]+10,np.array(gate1)[0,0]-1,'A1',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
axs_map.text(np.array(gate2)[0,1]-3,np.array(gate2)[0,0],'G2',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
axs_map.text(np.array(A2)[0,1]+3,np.array(A2)[0,0]+0.2,'A2',color = 'red',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')

#################################
########### - Plots - ###########
#################################
MBG1 = []
MBG2 = []
MBG3 = []

SIMA1 = []
SIMA2 = []

SA2 = []

for file in os.listdir('ORAS5/Data/iicethic'):
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        a,b = MB_G1G2(file)
        MBG1.append(a*1e-9)
        MBG2.append(b*1e-9)
        a,b = SIM_A1A2(file)
        SIMA1.append(a*1e-9)
        SIMA2.append(b*1e-9)
        SA2.append(S_A2(file,50))

MBG1 = np.convolve(MBG1,np.ones(12)/12, mode = 'valid')
MBG2 = np.convolve(MBG2,np.ones(12)/12, mode = 'valid')

SIMA1 = np.convolve(SIMA1,np.ones(12)/12, mode = 'valid')
SIMA2 = np.convolve(SIMA2,np.ones(12)/12, mode = 'valid')

SA2 = np.convolve(SA2,np.ones(12)/12, mode = 'valid')

MBG1linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],MBG1)
MBG2linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],MBG2)

SIMA1linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],SIMA1)
SIMA2linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],SIMA2)

SA2linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],SA2)

axs_plot[0].plot(np.linspace(2012,2019,len(MBG1)),MBG1,linewidth = 2,label = 'G1', color = 'blue')
axs_plot[0].plot(np.linspace(2012,2019,len(MBG2)),MBG2,linewidth = 2,label = 'G2', color = 'red')
#axs_plot[0].plot(np.linspace(2012,2019,len(MBG3)),MBG3,label = 'G3', color = 'green')

axs_plot[0].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),MBG1linregress.intercept + MBG1linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'blue')
axs_plot[0].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),MBG2linregress.intercept + MBG2linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'red')
#axs_plot[0].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),MBG3linregress.intercept + MBG3linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'green')
axs_plot[0].set_title(f'Sea ice transport through G1 and G2',fontsize = titlesize)
axs_plot[0].legend(fontsize = legendsize)
axs_plot[0].grid()
axs_plot[0].set_ylabel(r'$[km^3*month^{-1}]$',size = labelsize)
axs_plot[0].tick_params(labelsize = ticklabelsize)

#axs_plot[1].plot(np.linspace(2012,2019,len(SIMA1)),SIMA1,linewidth = 2,label = 'A1', color = 'blue')
axs_plot[1].plot(np.linspace(2012,2019,len(SIMA2)),SIMA2,linewidth = 2,label = 'A2', color = 'red')

#axs_plot[1].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),SIMA1linregress.intercept + SIMA1linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'blue')
axs_plot[1].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),SIMA2linregress.intercept + SIMA2linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'red')
axs_plot[1].set_title(f'Sea ice melt over A2',fontsize = titlesize)
axs_plot[1].legend(fontsize = legendsize)
axs_plot[1].grid()
axs_plot[1].set_ylabel(r'$[km^3*month^{-1}]$',size = labelsize)
axs_plot[1].tick_params(labelsize = ticklabelsize)

axs_plot[2].plot(np.linspace(2012,2019,len(SA2)),SA2,linewidth = 2,label = 'A2', color = 'red')

axs_plot[2].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),SA2linregress.intercept + SA2linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'red')
axs_plot[2].set_title(f'Salinity at 50m depth averaged over A2',fontsize = titlesize)
#axs_plot[2].legend(fontsize = legendsize)
axs_plot[2].grid()
axs_plot[2].set_ylabel(r'$[psu]$',size = labelsize)
axs_plot[2].tick_params(labelsize = ticklabelsize)

plt.savefig('MT_plot/MB_SIM_S50m.png')
plt.close()
