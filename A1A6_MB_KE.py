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

ke             = []
day_btwn_2011_start_and_2018_end = []

for file in os.listdir('ORAS5/Data/ke'):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        ke.append(np.where(EGC_path == 1,np.loadtxt(f'ORAS5/Data/ke/{file}')*1e+3,np.nan))
        day_btwn_2011_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)

ke = np.nan_to_num(np.array(ke))


trend_ke = np.zeros((np.shape(ke)[1:]))
p_ke     = np.zeros((np.shape(ke)[1:]))

for l in range(np.shape(trend_ke)[0]):
    for c in range(np.shape(trend_ke)[1]):

        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(ke[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_ke[l,c]             = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(ke[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_ke[l,c]         = slope * (365)
        

trend_ke[trend_ke == 0] = np.nan
p_ke[p_ke == 1]         = np.nan





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

def KE(coords, file):
    ke = np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]) &
                  (EGC_path == 1),
                  np.loadtxt(f'ORAS5/Data/ke/{file}'),
                  np.nan)


    return np.nanmean(ke)


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

KE_list      = [[],[],[],[],[],[]]

for coords,i in zip([A1_corners,A2_corners,A3_corners,A4_corners,A5_corners,A6_corners],range(6)):
    for file in os.listdir('ORAS5/Data/iicethic'):
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
            KE_list[i].append(KE(coords,file)*1e+3)
    KE_list[i]      = np.convolve(KE_list[i], np.ones(12)/12, mode = 'valid')

    
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
axs_map1.set_title('East Greenland Current\n kinetic energy trend',fontsize = 35)        

cs = axs_map1.pcolormesh(lon, lat, trend_ke,vmin = -1.4, vmax = 1.4,cmap = "cmo.balance", transform=ccrs.PlateCarree())

coords_p_signi = []
for k in range(len(p_ke)):
    for l in range(len(p_ke[0])):
        if p_ke[k,l] <= 0.05:
            coords_p_signi.append([lon[k,l],lat[k,l]])
coords_p_signi = np.array(coords_p_signi)
if len(coords_p_signi) != 0:
    cs_non_signi = axs_map1.scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())

cax = fig.add_axes([axs_map1.get_position().x1+0.01,axs_map1.get_position().y0+0.05,0.015,axs_map1.get_position().height*0.95])
cb = plt.colorbar(cs, cax = cax)
cb.ax.tick_params(labelsize=labelsize)
cb.set_label(label = r'[$mJ\cdot kg^{-1} \cdot y^{-1}$]',size = labelsize)
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
KE_list_memory = np.copy(KE_list)
for i in range(6):
    if lag_time[i] != 0:
        KE_list[i] = KE_list[i][lag_time[i]:]
        MB_list[0] = MB_list_memory[:-lag_time[i]]

    axs_plot[i].plot(np.linspace(2012 + lag_time[i]/12,2019,len(MB_list[0])),MB_list[0],linewidth = 3, color = 'blue')
    #linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],MB_list[0])
    #axs_plot[i].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),linregress.intercept + linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'blue')

    r = stats.pearsonr(KE_list[i],MB_list[0])[0]
    p = stats.pearsonr(KE_list[i],MB_list[0])[1] 
    if i == 0:    
        axs_plot[i].set_title(f'Correlation between sea ice transport through G1\nand EGC kinetic energy over\nA{i+1}, r: {round(r,2)}, p: {round(p,3)}',fontsize = titlesize)
    elif i == 1:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)},'+r' $p < 10^{-7}$',fontsize = titlesize)
    elif i == 2:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)},'+r' $p < 10^{-5}$',fontsize = titlesize)
    elif i == 4:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)},'+r' $p < 10^{-5}$',fontsize = titlesize)
    elif i == 5:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)},'+r' $p < 10^{-6}$',fontsize = titlesize)
    else:
        axs_plot[i].set_title(f'A{i+1}, r: {round(r,2)}, p: {round(p,3)}',fontsize = titlesize)
    

    if lag_time[i] == 1:
        axs_plot[i].text(2016.6,120,f'lag time: {lag_time[i]} month',fontsize = legendsize,weight='bold')
    elif lag_time[i] != 0 and i != 5 and i != 4:
        axs_plot[i].text(2016.5,120,f'lag time: {lag_time[i]} months',fontsize = legendsize,weight='bold')
    elif i == 5 or i == 4:
        axs_plot[i].text(2015.9,50,f'lag time: {lag_time[i]} months',fontsize = legendsize,weight='bold')
    else:
        axs_plot[i].text(2016.6,120,f'no lag time',fontsize = legendsize,weight='bold')


    ax_sal = plt.twinx(axs_plot[i])
    ax_sal.plot(np.linspace(2012,2019-lag_time[i]/12,len(KE_list[i])),KE_list[i],linewidth = 2,label = f'A1', color = 'red')

    ax_sal.tick_params(labelsize = ticklabelsize, color = 'red')
    ax_sal.tick_params(axis='y', colors='red')
    axs_plot[i].grid()
    if i != 3 and i != 5:
        axs_plot[i].set_ylabel(r'$[km^3\cdot month^{-1}]$',size = labelsize, color = 'blue')
    if i != 2 and i != 4:
        ax_sal.set_ylabel(r'$[mJ \cdot kg^{-1}]$',size = labelsize, color = 'red')

    axs_plot[i].tick_params(labelsize = ticklabelsize, color = 'blue')
    axs_plot[i].tick_params(axis='y', colors='blue')
    if i < 4:
        axs_plot[i].set_xticklabels([])

    linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],KE_list_memory[i])
    print(f'A{i+1}: slope:{linregress.slope *365}mJ/kg/year, p_lin:{linregress.pvalue}, p_corr:{p}')

plt.savefig('MT_plot/correlation/A1A6_MB_KE.png')
plt.close()

