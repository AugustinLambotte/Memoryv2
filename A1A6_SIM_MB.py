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
import matplotlib.ticker as ticker


lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")
surface   = np.loadtxt('ORAS5/Data/geo/Surface.txt')
X_dist    = np.loadtxt('ORAS5/Data/geo/X_dist.txt')
Y_dist    = np.loadtxt('ORAS5/Data/geo/Y_dist.txt')
earth_radius = 6370 * 1e3
EGC_path = np.loadtxt('ORAS5/Data/geo/Current_path.txt')

sit             = []
sic             = []
day_btwn_2011_start_and_2018_end = []

for file in os.listdir('ORAS5/Data/ke'):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        sit.append(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        sic.append(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))
        day_btwn_2011_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)

sic = np.nan_to_num(np.array(sic))
sit = np.nan_to_num(np.array(sit))


trend_sit = np.zeros((np.shape(sic)[1:]))
trend_sic = np.zeros((np.shape(sic)[1:]))
p_sit     = np.zeros((np.shape(sic)[1:]))
p_sic     = np.zeros((np.shape(sic)[1:]))

for l in range(np.shape(trend_sit)[0]):
    for c in range(np.shape(trend_sit)[1]):

        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sit[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sit[l,c]             = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sit[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sit[l,c]         = slope * (365)
        
        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sic[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sic[l,c]             = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(sic[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sic[l,c]         = slope * (365)

trend_sit[trend_sit == 0] = np.nan
trend_sic[trend_sic == 0] = np.nan
p_sit[p_sit == 1]         = np.nan
p_sic[p_sic == 1]         = np.nan





def MB(coords, file):
    iicevelv = np.loadtxt(f'ORAS5/Data/iicevelv/{file}')
    iicevelu = np.loadtxt(f'ORAS5/Data/iicevelu/{file}')
    sit      = np.loadtxt(f'ORAS5/Data/iicethic/{file}')
    sic      = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')

    ##### - GATE 1 - #####

    fromNorth = 0
    for lon_ in range(int(coords[0,1]),int(coords[1,1])):
        sit_interp = (np.nan_to_num(float(sit[     (lat == coords[0,0]+0.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == coords[0,0]-0.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == coords[0,0]+0.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == coords[0,0]-0.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == coords[0,0]+0.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == coords[0,0]-0.5) & (lon == lon_)])))/2

        fromNorth -= sit_interp * sic_interp * v_interp * float(X_dist[(lat == coords[0,0]+0.5) & (lon == lon_)])
   
    MB = fromNorth*60*60*24*30 #m^3/month
 
    return MB

def SIM(coords, file):
    iicevelv = np.loadtxt(f'ORAS5/Data/iicevelv/{file}')
    iicevelu = np.loadtxt(f'ORAS5/Data/iicevelu/{file}')
    sit      = np.loadtxt(f'ORAS5/Data/iicethic/{file}')
    sic      = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')

    fromNorth = 0
    fromSouth = 0
    for lon_ in range(int(coords[0,1]),int(coords[1,1])):
        
        sit_interp = (np.nan_to_num(float(sit[     (lat == coords[1,0]-0.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == coords[1,0]+0.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == coords[1,0]-0.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == coords[1,0]+0.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == coords[1,0]-0.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == coords[1,0]+0.5) & (lon == lon_)])))/2

        fromNorth -= sit_interp * sic_interp * v_interp * float(X_dist[(lat == coords[1,0]-0.5) & (lon == lon_)])
   
        sit_interp = (np.nan_to_num(float(sit[     (lat == coords[0,0]-0.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == coords[0,0]+0.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == coords[0,0]-0.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == coords[0,0]+0.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == coords[0,0]-0.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == coords[0,0]+0.5) & (lon == lon_)])))/2

        fromSouth += sit_interp * sic_interp * v_interp * float(X_dist[(lat == coords[0,0]-0.5) & (lon == lon_)])

    fromWest = 0
    fromEast = 0
    for lat_ in [coords[0,0]+0.5,coords[1,0]-0.5]:
        sit_interp = (np.nan_to_num(float(sit[     (lat == lat_) & (lon == int(coords[0,1]))])) + np.nan_to_num(float(sit[     (lat == lat_) & (lon == int(coords[0,1])-1)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == lat_) & (lon == int(coords[0,1]))])) + np.nan_to_num(float(sic[     (lat == lat_) & (lon == int(coords[0,1])-1)])))/2
        u_interp   = (np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == int(coords[0,1]))])) + np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == int(coords[0,1])-1)])))/2
        
        fromWest += sit_interp * sic_interp * u_interp * float(Y_dist[(lat == lat_) & (lon == lon_)])

        sit_interp = (np.nan_to_num(float(sit[     (lat == lat_) & (lon == int(coords[1,1]))])) + np.nan_to_num(float(sit[     (lat == lat_) & (lon == int(coords[1,1])-1)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == lat_) & (lon == int(coords[1,1]))])) + np.nan_to_num(float(sic[     (lat == lat_) & (lon == int(coords[1,1])-1)])))/2
        u_interp   = (np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == int(coords[1,1]))])) + np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == int(coords[1,1])-1)])))/2

        fromEast -= sit_interp * sic_interp * u_interp * float(Y_dist[(lat == lat_) & (lon == lon_)])
    
    
    MB = (fromNorth + fromSouth + fromEast + fromWest)*60*60*24*30 #m^3/month

    year = int(file[:4])
    month = int(file[5:7])
    if year == 2011 and month == 1:
        previousSIT = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
        nextSIT     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/2011-02.txt'))

        previousSIC = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))
        nextSIC     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/2011-02.txt'))
    elif year == 2018 and month == 12:
        previousSIT = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/2018-11.txt'))
        nextSIT     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))

        previousSIC = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/2018-11.txt'))
        nextSIC     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))
    elif month == 1:
        previousSIT = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{year-1}-12.txt'))
        nextSIT     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{year  }-02.txt'))

        previousSIC = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{year-1}-12.txt'))
        nextSIC     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{year  }-02.txt'))
    elif month == 12:
        previousSIT = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{year  }-11.txt'))
        nextSIT     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{year+1}-01.txt'))

        previousSIC = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{year  }-11.txt'))
        nextSIC     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{year+1}-01.txt'))
    else:
        previousSIT = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{year}-{month-1:02d}.txt'))
        nextSIT     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{year}-{month+1:02d}.txt'))

        previousSIC = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{year}-{month-1:02d}.txt'))
        nextSIC     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{year}-{month+1:02d}.txt'))
    
    previousSIT = np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),previousSIT,0)
    previousSIC = np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),previousSIC,0)
    nextSIT     = np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),nextSIT,    0)
    nextSIC     = np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),nextSIC,    0)


    startMonthSIV = (np.nansum(previousSIT*previousSIC*surface) + np.nansum(np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),sit*sic*surface,0)))/2
    endMonthSIV   = (np.nansum(nextSIT*nextSIC*surface)         + np.nansum(np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),sit*sic*surface,0)))/2

    bilanSIV = endMonthSIV - startMonthSIV #Positive when SIV increase over area

    SIM = MB - bilanSIV

    return SIM


#################
####- GATE - ####
#################

gate1_coord = np.array([[79,-19.5],[79,10.5]])
gate2_coord = np.array([[76,-19.5],[76,2.5]])
gate3_coord = np.array([[73,-20.5],[73,-3.5]])
gate4_coord = np.array([[70,-21.5],[70,-9.5]])

gate1 = [[gate1_coord[0,0],i] for i in np.linspace(gate1_coord[0,1],gate1_coord[1,1],10)]
gate2 = [[gate2_coord[0,0],i] for i in np.linspace(gate2_coord[0,1],gate2_coord[1,1],10)]
gate3 = [[gate3_coord[0,0],i] for i in np.linspace(gate3_coord[0,1],gate3_coord[1,1],10)]
gate4 = [[gate4_coord[0,0],i] for i in np.linspace(gate4_coord[0,1],gate4_coord[1,1],10)]

MB_list       = [[],[],[],[]]

for coords,i in zip([gate1_coord,gate2_coord,gate3_coord,gate4_coord],range(4)):
    for file in os.listdir('ORAS5/Data/iicethic'):
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
            MB_list[i].append(MB(coords,file)*1e-9)
    MB_list[i] = np.convolve(MB_list[i], np.ones(12)/12, mode = 'valid')



##################
#### - ZONE - ####
##################
A1_corners = np.array([[77,-16.5],[79,-0.5]])
A2_corners = np.array([[75,-17.5],[77,-2.5]])
A3_corners = np.array([[73,-18.5],[75,-10.5]])
A4_corners = np.array([[71,-20.5],[73,-12.5]])
A5_corners = np.array([[69,-21.5],[71,-14.5]])
A6_corners = np.array([[67,-27.5],[69,-19.5]])


A1 = np.concatenate(([[A1_corners[0,0],i] for i in np.linspace(A1_corners[0,1],A1_corners[1,1],10)],[[A1_corners[1,0],A1_corners[1,1]]],[[A1_corners[1,0],i] for i in np.linspace(A1_corners[1,1],A1_corners[0,1],10)],[[A1_corners[0,0],A1_corners[0,1]]]))
A2 = np.concatenate(([[A2_corners[0,0],i] for i in np.linspace(A2_corners[0,1],A2_corners[1,1],10)],[[A2_corners[1,0],A2_corners[1,1]]],[[A2_corners[1,0],i] for i in np.linspace(A2_corners[1,1],A2_corners[0,1],10)],[[A2_corners[0,0],A2_corners[0,1]]]))
A3 = np.concatenate(([[A3_corners[0,0],i] for i in np.linspace(A3_corners[0,1],A3_corners[1,1],10)],[[A3_corners[1,0],A3_corners[1,1]]],[[A3_corners[1,0],i] for i in np.linspace(A3_corners[1,1],A3_corners[0,1],10)],[[A3_corners[0,0],A3_corners[0,1]]]))
A4 = np.concatenate(([[A4_corners[0,0],i] for i in np.linspace(A4_corners[0,1],A4_corners[1,1],10)],[[A4_corners[1,0],A4_corners[1,1]]],[[A4_corners[1,0],i] for i in np.linspace(A4_corners[1,1],A4_corners[0,1],10)],[[A4_corners[0,0],A4_corners[0,1]]]))
A5 = np.concatenate(([[A5_corners[0,0],i] for i in np.linspace(A5_corners[0,1],A5_corners[1,1],10)],[[A5_corners[1,0],A5_corners[1,1]]],[[A5_corners[1,0],i] for i in np.linspace(A5_corners[1,1],A5_corners[0,1],10)],[[A5_corners[0,0],A5_corners[0,1]]]))
A6 = np.concatenate(([[A6_corners[0,0],i] for i in np.linspace(A6_corners[0,1],A6_corners[1,1],10)],[[A6_corners[1,0],A6_corners[1,1]]],[[A6_corners[1,0],i] for i in np.linspace(A6_corners[1,1],A6_corners[0,1],10)],[[A6_corners[0,0],A6_corners[0,1]]]))


SIM_list      = [[],[],[],[],[],[]]

for coords,i in zip([A1_corners,A2_corners,A3_corners,A4_corners,A5_corners,A6_corners],range(6)):
    for file in os.listdir('ORAS5/Data/iicethic'):
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
            SIM_list[i].append(SIM(coords,file)*1e-9)
    SIM_list[i]      = np.convolve(SIM_list[i], np.ones(12)/12, mode = 'valid')

    
####################
#### - Figure - ####
####################

labelsize = 28
titlesize = 27
ticklabelsize = 25
legendsize = 20

fig = plt.figure(figsize = (22,19))
# add grid specifications
gs = fig.add_gridspec(4,2)
axs_map1 = fig.add_subplot(gs[:2,1], projection = ccrs.LambertConformal(central_longitude = -9))
axs_plot = []
for i in range(3):
    axs_plot.append(fig.add_subplot(gs[i,0]))
axs_plot.append(fig.add_subplot(gs[2,1]))
axs_plot.append(fig.add_subplot(gs[3,0]))
axs_plot.append(fig.add_subplot(gs[3,1]))
plt.subplots_adjust(top = 0.93, bottom = 0.05,left = 0.07,right = 0.9)#, hspace = 0.3, left = 0.01,wspace=0.1,right =0.95)

##############################
#### - plotting the map - ####
##############################
xlim = [-30, 13]
ylim = [65, 80]
lower_space = 3 
rect = mpath.Path([[xlim[0], ylim[0]],
                [xlim[1], ylim[0]],
                [xlim[1], ylim[1]],
                [xlim[0], ylim[1]],
                [xlim[0], ylim[0]],
                ]).interpolated(20)
proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axs_map1) - axs_map1.transData
rect_in_target = proj_to_data.transform_path(rect)
axs_map1.set_boundary(rect_in_target)
axs_map1.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

axs_map1.coastlines()
axs_map1.gridlines()        
axs_map1.set_title('Sea ice thickness\ntrend',fontsize = 35)        

cs = axs_map1.pcolormesh(lon, lat, trend_sit,vmin = -0.15, vmax = 0.15,cmap = "cmo.balance", transform=ccrs.PlateCarree())

coords_p_signi = []
for k in range(len(p_sit)):
    for l in range(len(p_sit[0])):
        if p_sit[k,l] <= 0.05:
            coords_p_signi.append([lon[k,l],lat[k,l]])
coords_p_signi = np.array(coords_p_signi)
if len(coords_p_signi) != 0:
    cs_non_signi = axs_map1.scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())

cax = fig.add_axes([axs_map1.get_position().x1+0.01,axs_map1.get_position().y0+0.05,0.015,axs_map1.get_position().height*0.95])
cb = plt.colorbar(cs, cax = cax)
cb.ax.tick_params(labelsize=labelsize)
cb.set_label(label = r'[$m\cdot y^{-1}$]',size = labelsize)
### - Transect - ###


axs_map1.plot(np.array(gate1)[:,1], np.array(gate1)[:,0], linewidth = 3,color = 'purple',transform=ccrs.PlateCarree())

axs_map1.text(np.array(gate1)[-1,1]+3,np.array(gate1)[0,0]+0.2,'G1',color = 'purple',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')

axs_map1.plot(np.array(A1)[:,1], np.array(A1)[:,0], linewidth = 2,color = 'black', linestyle = 'dashed',transform=ccrs.PlateCarree())
axs_map1.plot(np.array(A2)[:,1], np.array(A2)[:,0], linewidth = 2,color = 'black', linestyle = 'dashed',transform=ccrs.PlateCarree())
axs_map1.plot(np.array(A3)[:,1], np.array(A3)[:,0], linewidth = 2,color = 'black', linestyle = 'dashed',transform=ccrs.PlateCarree())
axs_map1.plot(np.array(A4)[:,1], np.array(A4)[:,0], linewidth = 2,color = 'black', linestyle = 'dashed',transform=ccrs.PlateCarree())
axs_map1.plot(np.array(A5)[:,1], np.array(A5)[:,0], linewidth = 2,color = 'black', linestyle = 'dashed',transform=ccrs.PlateCarree())
axs_map1.plot(np.array(A6)[:,1], np.array(A6)[:,0], linewidth = 2,color = 'black', linestyle = 'dashed',transform=ccrs.PlateCarree())

axs_map1.text(np.array(A1)[0,1]+3,np.array(A1)[0,0]+0.2,'A1',color = 'black',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
axs_map1.text(np.array(A2)[0,1]+1,np.array(A2)[0,0]+0.2,'A2',color = 'black',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
axs_map1.text(np.array(A3)[0,1]+1,np.array(A3)[0,0]+0.2,'A3',color = 'black',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
axs_map1.text(np.array(A4)[0,1]+1,np.array(A4)[0,0]+0.2,'A4',color = 'black',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
axs_map1.text(np.array(A5)[0,1]+1,np.array(A5)[0,0]+0.2,'A5',color = 'black',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
axs_map1.text(np.array(A6)[0,1]+1,np.array(A6)[0,0]+0.2,'A6',color = 'black',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')



#################################
########### - Plots - ###########
#################################

lag_time = [0,1,2,3,4,5]
MB_list_memory = MB_list[0]
SIM_list_memory = np.copy(SIM_list)
for i in range(6):
    """ if lag_time[i] != 0:
        SIM_list[i] = SIM_list[i][lag_time[i]:]
        MB_list[0] = MB_list_memory[:-lag_time[i]] """

    axs_plot[i].plot(np.linspace(2012 + lag_time[i]/12,2019,len(MB_list[0])),MB_list[0],linewidth = 3, color = 'blue')
    if i != 0:
        print(len(SIM_list[i][lag_time[i]:]))
        print(len(MB_list[0][:-lag_time[i]]))

        r = stats.pearsonr(SIM_list[i][lag_time[i]:],MB_list[0][:-lag_time[i]])[0]
        p = stats.pearsonr(SIM_list[i][lag_time[i]:],MB_list[0][:-lag_time[i]])[1]   
    else:
        r = stats.pearsonr(SIM_list[i],MB_list[0])[0]
        p = stats.pearsonr(SIM_list[i],MB_list[0])[1]   
    if i == 0:    
        axs_plot[i].set_title(f'Correlation between sea ice transport\nthrough G1 and sea ice melt over\nA{i+1}, r: {round(r,2)},'+r' $p < 10^{-10}$',fontsize = titlesize)
    elif i == 1:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)},'+r' $p < 10^{-7}$',fontsize = titlesize)
    elif i == 2:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)},'+r' $p < 10^{-3}$',fontsize = titlesize)
    elif i == 3:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)},'+r' $p < 10^{-16}$',fontsize = titlesize)
    elif i == 4:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)},'+r' $p < 10^{-11}$',fontsize = titlesize)
    else:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)}, p: {round(p,3)}',fontsize = titlesize)
    
    if lag_time[i] == 1:
        axs_plot[i].text(2016.5,120,f'lag time: {lag_time[i]} month',fontsize = legendsize,weight='bold')
    elif lag_time[i] != 0 and i != 5:
        axs_plot[i].text(2016.5,120,f'lag time: {lag_time[i]} months',fontsize = legendsize,weight='bold')
    elif i == 5:
        axs_plot[i].text(2015.5,60,f'lag time: {lag_time[i]} months',fontsize = legendsize,weight='bold')
    else:
        axs_plot[i].text(2016.5,120,f'no lag time',fontsize = legendsize,weight='bold')

    ax_sal = plt.twinx(axs_plot[i])
    ax_sal.plot(np.linspace(2012,2019-lag_time[i]/12,len(SIM_list[i])),SIM_list[i],linewidth = 2,label = f'A1', color = 'red')

    

    ax_sal.tick_params(labelsize = ticklabelsize, color = 'red')
    ax_sal.tick_params(axis='y', colors='red')
    axs_plot[i].grid()
    if i != 3 and i != 5:
        axs_plot[i].set_ylabel(r'$[km^3\cdot month^{-1}]$',size = labelsize, color = 'blue')
    if i != 2 and i != 4:
        ax_sal.set_ylabel(r'$[km^3\cdot month^{-1}]$',size = labelsize, color = 'red')

    axs_plot[i].tick_params(labelsize = ticklabelsize, color = 'blue')
    axs_plot[i].tick_params(axis='y', colors='blue')
    if i < 4:
        axs_plot[i].set_xticklabels([])

    linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],SIM_list_memory[i])
    print(f'A{i+1}: slope:{linregress.slope *365}km^3/month/year, p_lin:{linregress.pvalue}, p_corr:{p}')

plt.savefig('MT_plot/correlation/A1A6_SIM_MB.png')
plt.close()

