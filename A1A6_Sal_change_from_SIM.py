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

earth_radius = 6370 * 1e3

lon          = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat          = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")
surface      = np.loadtxt('ORAS5/Data/geo/Surface.txt')
X_dist       = np.loadtxt('ORAS5/Data/geo/X_dist.txt')
Y_dist       = np.loadtxt('ORAS5/Data/geo/Y_dist.txt')
EGC_path     = np.loadtxt('ORAS5/Data/geo/Current_path.txt')

depth = 100
S_    = []

day_btwn_2011_start_and_2018_end = []
for file in os.listdir('ORAS5/Data/ke'):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        
        
        day_btwn_2011_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)
        Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
        Smap = Sds['vosaline'].where(Sds.deptht < depth, drop = True)[:]
        if file == '2011-01.txt':
            # height list = list og the height in meters of each cells
            depth_list = np.array(Sds.deptht)
            height_list = [2*depth_list[0]]
            for j in range(1,len(depth_list)):
                if depth_list[j] > depth:
                    break

                height_list.append(2*(depth_list[j] - np.sum(height_list)))
        S_.append(np.apply_along_axis(lambda x : np.nansum(x*height_list)/np.sum(height_list), axis = 0, arr = Smap))
        Sds.close()

S_       = np.nan_to_num(np.array(S_))

trend_S= np.zeros((np.shape(S_)[1:]))
p_S     = np.zeros((np.shape(S_)[1:]))

for l in range(np.shape(trend_S)[0]):
    for c in range(np.shape(trend_S)[1]):
        slope         = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S_[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_S[l,c]     = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S_[:,l,c],np.ones(12)/12,  mode = 'valid')).pvalue
        trend_S[l,c] = slope * (365)



trend_S[trend_S == 0] = np.nan
p_S[p_S == 1]         = np.nan
trend_S[trend_S<-0.1] = np.nan
trend_S[trend_S> 0.1] = np.nan




def WF(coords,file):
    WF = np.loadtxt(f'ORAS5/Data/sowaflup/{file}')
    WF = np.nansum(np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),
                            WF * surface,
                            np.nan))
    WF *= 60*60*24*30 # [kg/month]
    WF /= 1000 # m^3/month
    

    return WF

def SIV(coords,file):
    sic = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))
    sit = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
    siv = sit * sic * surface
    siv = np.nansum(np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),siv,np.nan))

    return siv

def S(coords,file,depth):
    Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
    Smap = Sds['vosaline'].where(Sds.deptht < depth, drop = True)[:]
    # height list = list og the height in meters of each cells
    depth_list = np.array(Sds.deptht)
    height_list = [2*depth_list[0]]
    for i in range(1,len(depth_list)):
        if depth_list[i] > depth:
            break
        height_list.append(2*(depth_list[i] - np.sum(height_list)))

    mean_depth_integrated_salinity = np.zeros(np.shape(Smap)[1:])
    for i in range(len(Smap[0])):
        for j in range(len(Smap[0][0])):
            mean_depth_integrated_salinity[i,j] = np.nansum(Smap[:,i,j] * height_list)/np.sum(height_list)
    
    
    S = np.nanmean(np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),
                            mean_depth_integrated_salinity,
                            np.nan))
    
    water_vol = np.nansum(np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),
                            surface * np.sum(height_list),
                            np.nan)) #m^3 of water in the area aboce the given depth

    return S, water_vol

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
    for lat_ in np.arange(coords[0,0]+0.5,coords[1,0]+0.5): # arange stops one before the last point
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


def S_change_due_to_SIM(S1,SIM1,SIM2, water_vol):
    water_mass = water_vol * 1000 #passing from m^3 to kg
    delta_sim = SIM2 - SIM1
    delta_sim_mass = delta_sim * 1000
    delta_S = ((S1*1e-3 * water_mass + 0.005 * delta_sim_mass)/(water_mass + delta_sim_mass) - S1*1e-3)*1e+3

    return delta_S


A1_corners = np.array([[76,-16.5],[79,-2.5]])
A2_corners = np.array([[68,-19.5],[71,-13.5]])


A1 = np.concatenate(([[A1_corners[0,0],i] for i in np.linspace(A1_corners[0,1],A1_corners[1,1],10)],[[A1_corners[1,0],A1_corners[1,1]]],[[A1_corners[1,0],i] for i in np.linspace(A1_corners[1,1],A1_corners[0,1],10)],[[A1_corners[0,0],A1_corners[0,1]]]))
A2 = np.concatenate(([[A2_corners[0,0],i] for i in np.linspace(A2_corners[0,1],A2_corners[1,1],10)],[[A2_corners[1,0],A2_corners[1,1]]],[[A2_corners[1,0],i] for i in np.linspace(A2_corners[1,1],A2_corners[0,1],10)],[[A2_corners[0,0],A2_corners[0,1]]]))

SIM_list  = [[],[]]
S_list    = [[],[]]
WF_list   = [[],[]]
SIV_list  = [[],[]]
water_vol = []

for coords,i in zip([A1_corners,A2_corners],range(2)):
    for file in os.listdir('ORAS5/Data/iicethic'):
        print(f'{i+1}/{len(SIM_list)} : {file}')
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
            SIM_list[i].append(SIM(coords,file)*1e-9)
            WF_list[i].append(WF(coords,file))
            SIV_list[i].append(SIV(coords,file))

            ####- mean salinity in the upper layer - ####
            Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
            Smap = Sds['vosaline'].where(Sds.deptht < depth, drop = True)[:]
            if file == '2011-01.txt':
                # height list = list og the height in meters of each cells
                depth_list = np.array(Sds.deptht)
                height_list = [2*depth_list[0]]
                for j in range(1,len(depth_list)):
                    if depth_list[j] > depth:
                        break

                    height_list.append(2*(depth_list[j] - np.sum(height_list)))
                water_vol.append(np.nansum(np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),
                            surface * np.sum(height_list),
                            np.nan))) #m^3 of water in the area above the given depth
                
            mean_depth_integrated_salinity = np.apply_along_axis(lambda x : np.nansum(x*height_list)/np.sum(height_list), axis = 0, arr = Smap)
            Sds.close()
            
            
            S_ = np.nanmean(np.where((lat > coords[0,0]) & (lat < coords[1,0]) & (lon < coords[1,1]) & (lon > coords[0,1]),
                                    mean_depth_integrated_salinity,
                                    np.nan))
            
            S_list[i].append(S_)

    SIM_list[i] = np.convolve(SIM_list[i], np.ones(12)/12, mode = 'valid')
    S_list[i]   = np.convolve(S_list[i]  , np.ones(12)/12, mode = 'valid')
    WF_list[i]  = np.convolve(WF_list[i] , np.ones(12)/12, mode = 'valid')
    SIV_list[i] = np.convolve(SIV_list[i] , np.ones(12)/12, mode = 'valid')







labelsize = 28
titlesize = 27
ticklabelsize = 25
legendsize = 20

fig = plt.figure(figsize = (22,10))
# add grid specifications
gs = fig.add_gridspec(2,2)
axs_map1 = fig.add_subplot(gs[:2,1], projection = ccrs.LambertConformal(central_longitude = -9))
axs_plot = []
for i in range(2):
    axs_plot.append(fig.add_subplot(gs[i,0]))
""" axs_plot.append(fig.add_subplot(gs[2,1]))
axs_plot.append(fig.add_subplot(gs[3,0]))
axs_plot.append(fig.add_subplot(gs[3,1])) """
plt.subplots_adjust(top = 0.88, bottom = 0.11,left = 0.07,right = 0.95,hspace = 0.2)#, hspace = 0.3, left = 0.01,wspace=0.1,right =0.95)

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
axs_map1.set_title('100m depth averaged\nsalinity trends',fontsize = 35)        

cs = axs_map1.pcolormesh(lon, lat, trend_S,vmin = -0.1, vmax = 0.1,cmap = "cmo.balance", transform=ccrs.PlateCarree())

coords_p_signi = []
for k in range(len(p_S)):
    for l in range(len(p_S[0])):
        if p_S[k,l] <= 0.05:
            coords_p_signi.append([lon[k,l],lat[k,l]])
coords_p_signi = np.array(coords_p_signi)
if len(coords_p_signi) != 0:
    cs_non_signi = axs_map1.scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())

cax = fig.add_axes([axs_map1.get_position().x0+0,axs_map1.get_position().y0+0,axs_map1.get_position().width,0.04])
cb = plt.colorbar(cs, cax = cax, ticks = [-0.1,0,0.1],orientation = 'horizontal')
cb.ax.tick_params(labelsize=labelsize)
cb.set_label(label = r'[$psu \cdot y^{-1}$]',size = labelsize)
### - Transect - ###

axs_map1.plot(np.array(A1)[:,1], np.array(A1)[:,0], linewidth = 2,color = 'black', linestyle = 'dashed',transform=ccrs.PlateCarree())
axs_map1.plot(np.array(A2)[:,1], np.array(A2)[:,0], linewidth = 2,color = 'black', linestyle = 'dashed',transform=ccrs.PlateCarree())

axs_map1.text(np.array(A1)[0,1]+3,np.array(A1)[0,0]+0.2,'B1',color = 'black',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
axs_map1.text(np.array(A2)[0,1]+1,np.array(A2)[0,0]+0.2,'B2',color = 'black',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')



#################################
########### - Plots - ###########
#################################

for i in range(2):
    

    axs_plot[i].plot(np.linspace(2012,2019,len(SIM_list[i])),SIM_list[i],linewidth = 2,label = f'A{i+1}', color = 'blue')
    

    #linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],MB_list[0])
    #axs_plot[i].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),linregress.intercept + linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'blue')

    r = stats.pearsonr(S_list[i],SIM_list[i])[0]
    p = stats.pearsonr(S_list[i],SIM_list[i])[1]
    if i == 0:
        axs_plot[i].set_title(f' B{i+1}, r: {round(r,2)}, p: {round(p,3)}',fontsize = titlesize)
    if i == 1:
        axs_plot[i].set_title(f' B{i+1}, r: {round(r,2)},'+r'$ p < 10^{-4}$',fontsize = titlesize)
    

    
    ax_sal = plt.twinx(axs_plot[i])
    ax_sal.plot(np.linspace(2012,2019,len(S_list[i])),S_list[i],linewidth = 2,label = f'A1', color = 'red')
    
    ax_sal.tick_params(labelsize = ticklabelsize, color = 'red')
    ax_sal.tick_params(axis='y', colors='red')
    axs_plot[i].grid()

    if i == 0:
        first = 39
        second = 52
        third = 69
        axs_plot[i].scatter(np.linspace(2012,2019,len(SIM_list[i]))[first],  SIM_list[i][first],  marker = 'v',s = 250, color = 'black')
        axs_plot[i].scatter(np.linspace(2012,2019,len(SIM_list[i]))[second], SIM_list[i][second], marker = '^',s = 250, color = 'black')
        axs_plot[i].scatter(np.linspace(2012,2019,len(SIM_list[i]))[third],  SIM_list[i][third],  marker = 'v',s = 250, color = 'black')

        axs_plot[i].text(2015.5,30.4,  r'$t^{B1}_1$',fontsize = 30, weight = 'bold',color = 'black')
        axs_plot[i].text(2016.5,28, r'$t^{B1}_2$',fontsize = 30, weight = 'bold',color = 'black')
        axs_plot[i].text(2018,30.4,  r'$t^{B1}_3$',fontsize = 30, weight = 'bold',color = 'black')

        axs_plot[i].plot([np.linspace(2012,2019,len(SIM_list[i]))[first],np.linspace(2012,2019,len(SIM_list[i]))[first]],[-10,50],linestyle = 'dashed', color = 'black')
        axs_plot[i].plot([np.linspace(2012,2019,len(SIM_list[i]))[second],np.linspace(2012,2019,len(SIM_list[i]))[second]],[-10,50],linestyle = 'dashed', color = 'black')
        axs_plot[i].plot([np.linspace(2012,2019,len(SIM_list[i]))[third],np.linspace(2012,2019,len(SIM_list[i]))[third]],[-10,50],linestyle = 'dashed', color = 'black')
        axs_plot[i].set_ylim(0,38)

        ax_sal.scatter(np.linspace(2012,2019,len(S_list[i]))[first],  S_list[i][first],  marker = '^',s = 250, color = 'black')
        ax_sal.scatter(np.linspace(2012,2019,len(S_list[i]))[second], S_list[i][second], marker = 'v',s = 250, color = 'black')
        ax_sal.scatter(np.linspace(2012,2019,len(S_list[i]))[third],  S_list[i][third],  marker = '^',s = 250, color = 'black')
        
        deltaS1   = S_list[i][second] - S_list[i][first]
        deltaS_induced_by_SIM1 = S_change_due_to_SIM(S_list[i][first],SIM_list[i][first]*1e+9,SIM_list[i][second]*1e+9,water_vol[i])

        print(f'### - B1 - ###\nSIM observed : {round(SIM_list[i][second] - SIM_list[i][first],3)}, S observed : {round(deltaS1,3)} calculated : {round(deltaS_induced_by_SIM1,3)} --- SIM explain S change by {round((deltaS_induced_by_SIM1/deltaS1)*100,2)}%')

        deltaS2   = S_list[i][third] - S_list[i][second]
        deltaS_induced_by_SIM2 = S_change_due_to_SIM(S_list[i][second],SIM_list[i][second]*1e+9,SIM_list[i][third]*1e+9,water_vol[i])

        print(f'SIM observed : {round(SIM_list[i][third] - SIM_list[i][second],3)}, S observed : {round(deltaS2,3)} calculated : {round(deltaS_induced_by_SIM2,3)} --- SIM explain S change by {round((deltaS_induced_by_SIM2/deltaS2)*100,2)}%')

    if i == 1:
        first  = 3
        second = 16
        third  = 30
        axs_plot[i].scatter(np.linspace(2012,2019,len(SIM_list[i]))[first],  SIM_list[i][first],  marker = '^',s = 250, color = 'black')
        axs_plot[i].scatter(np.linspace(2012,2019,len(SIM_list[i]))[second], SIM_list[i][second], marker = 'v',s = 250, color = 'black')
        axs_plot[i].scatter(np.linspace(2012,2019,len(SIM_list[i]))[third],  SIM_list[i][third],  marker = '^',s = 250, color = 'black')

        axs_plot[i].text(2012.5,14.8,  r'$t^{B2}_1$', fontsize = 30, weight = 'bold',color = 'black')
        axs_plot[i].text(2013.5,14.8,r'$t^{B2}_2$', fontsize = 30, weight = 'bold',color = 'black')
        axs_plot[i].text(2014.7,14.8,  r'$t^{B2}_3$', fontsize = 30, weight = 'bold',color = 'black')

        axs_plot[i].plot([np.linspace(2012,2019,len(SIM_list[i]))[first],np.linspace(2012,2019,len(SIM_list[i]))[first]],[-10,50],linestyle = 'dashed', color = 'black')
        axs_plot[i].plot([np.linspace(2012,2019,len(SIM_list[i]))[second],np.linspace(2012,2019,len(SIM_list[i]))[second]],[-10,50],linestyle = 'dashed', color = 'black')
        axs_plot[i].plot([np.linspace(2012,2019,len(SIM_list[i]))[third],np.linspace(2012,2019,len(SIM_list[i]))[third]],[-10,50],linestyle = 'dashed', color = 'black')
        axs_plot[i].set_ylim(0,18)
        ax_sal.scatter(np.linspace(2012,2019,len(S_list[i]))[first],  S_list[i][first],  marker = 'v',s = 250, color = 'black')
        ax_sal.scatter(np.linspace(2012,2019,len(S_list[i]))[second], S_list[i][second], marker = '^',s = 250, color = 'black')
        ax_sal.scatter(np.linspace(2012,2019,len(S_list[i]))[third],  S_list[i][third],  marker = 'v',s = 250, color = 'black')
        
        deltaS1   = S_list[i][second] - S_list[i][first]
        deltaS_induced_by_SIM1 = S_change_due_to_SIM(S_list[i][first],SIM_list[i][first]*1e+9,SIM_list[i][second]*1e+9,water_vol[i])

        print(f'### - B2 - ###\nSIM observed : {round(SIM_list[i][second] - SIM_list[i][first],3)}, S observed : {round(deltaS1,3)} calculated : {round(deltaS_induced_by_SIM1,3)} --- SIM explain S change by {round((deltaS_induced_by_SIM1/deltaS1)*100,2)}%')

        deltaS2   = S_list[i][third] - S_list[i][second]
        deltaS_induced_by_SIM2 = S_change_due_to_SIM(S_list[i][second],SIM_list[i][second]*1e+9,SIM_list[i][third]*1e+9,water_vol[i])

        print(f'SIM observed : {round(SIM_list[i][third] - SIM_list[i][second],3)}, S observed : {round(deltaS2,3)} calculated : {round(deltaS_induced_by_SIM2,3)} --- SIM explain S change by {round((deltaS_induced_by_SIM2/deltaS2)*100,2)}%')

    
    axs_plot[i].set_ylabel(r'$[km^3\cdot month^{-1}]$',size = labelsize, color = 'blue')
    ax_sal.set_ylabel(r'$[psu]$',size = labelsize, color = 'red')

    axs_plot[i].tick_params(labelsize = ticklabelsize, color = 'blue')
    axs_plot[i].tick_params(axis='y', colors='blue')
    if i == 0:
        axs_plot[i].set_xticklabels([])
    
    linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],SIM_list[i])
    print(f'slope:{linregress.slope *365}km^3/month/year, p_lin:{linregress.pvalue}, p_corr:{p}\n')
""" 
for i in [2,3]:
    

    axs_plot[i].plot(np.linspace(2012,2019,len(SIM_list[i-2])),SIM_list[i-2],linewidth = 2,label = f'A{i-1}', color = 'blue')
    

    #linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],MB_list[0])
    #axs_plot[i].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),linregress.intercept + linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'blue')

    r = stats.pearsonr(WF_list[i-2],SIM_list[i-2])[0]
    p = stats.pearsonr(WF_list[i-2],SIM_list[i-2])[1]
    axs_plot[i].set_title(f' A{i-1}, r: {round(r,2)},'+r'$ p < 10^{-18}$',fontsize = titlesize)
    

    
    ax_sal = plt.twinx(axs_plot[i])
    ax_sal.plot(np.linspace(2012,2019,len(WF_list[i-2])),WF_list[i-2],linewidth = 2,label = f'A1', color = 'red')
    
    ax_sal.tick_params(labelsize = ticklabelsize, color = 'red')
    ax_sal.tick_params(axis='y', colors='red')
    axs_plot[i].grid()


    axs_plot[i].tick_params(labelsize = ticklabelsize, color = 'blue')
    axs_plot[i].tick_params(axis='y', colors='blue')
    
    linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],WF_list[i-2])
    print(f'slope:{linregress.slope *365}km^3/month/year, p_lin:{linregress.pvalue}, p_corr:{p}\n')

for i in [4,5]:
    

    axs_plot[i].plot(np.linspace(2012,2019,len(SIV_list[i-4])),SIV_list[i-4],linewidth = 2,label = f'A{i-1}', color = 'blue')
    

    #linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],MB_list[0])
    #axs_plot[i].plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),linregress.intercept + linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'blue')

    r = stats.pearsonr(WF_list[i-4],SIV_list[i-4])[0]
    p = stats.pearsonr(WF_list[i-4],SIV_list[i-4])[1]
    axs_plot[i].set_title(f' A{i-1}, r: {round(r,2)},'+r'$ p < 10^{-18}$',fontsize = titlesize)
    

    
    ax_sal = plt.twinx(axs_plot[i])
    ax_sal.plot(np.linspace(2012,2019,len(WF_list[i-4])),WF_list[i-4],linewidth = 2,label = f'A1', color = 'red')
    
    ax_sal.tick_params(labelsize = ticklabelsize, color = 'red')
    ax_sal.tick_params(axis='y', colors='red')
    axs_plot[i].grid()


    axs_plot[i].tick_params(labelsize = ticklabelsize, color = 'blue')
    axs_plot[i].tick_params(axis='y', colors='blue')
    
    linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],WF_list[i-4])
    print(f'slope:{linregress.slope *365}km^3/month/year, p_lin:{linregress.pvalue}, p_corr:{p}\n')
 """

plt.savefig('MT_plot/correlation/A1A6_S_due_to_SIM.png')
plt.close()