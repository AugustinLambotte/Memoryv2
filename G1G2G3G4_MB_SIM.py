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

def SIM(file):
    iicevelv = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicevelv/{file}'))
    iicevelu = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicevelu/{file}'))
    sit      = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{file}'))
    sic      = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{file}'))
    # Net bilan of sea ice import in each cell
    transport = np.zeros(np.shape(iicevelv))
    for line in range(1,np.shape(iicevelv)[0]-1):
        for col in range(1,np.shape(iicevelv)[1]-1):

            #####################
            ### - TRANSPORT - ###
            #####################

            from_north_si_drift = - (iicevelv[line+1,col] + iicevelv[line,col])/2 # [m/s]
            from_west_si_drift  =   (iicevelu[line,col-1] + iicevelu[line,col])/2 # [m/s]
            from_east_si_drift  = - (iicevelu[line,col+1] + iicevelu[line,col])/2 # [m/s]
            from_south_si_drift =   (iicevelv[line-1,col] + iicevelv[line,col])/2 # [m/s]

            north_border_sit = (sit[line+1,col] + sit[line,col])/2
            south_border_sit = (sit[line-1,col] + sit[line,col])/2
            west_border_sit  = (sit[line,col-1] + sit[line,col])/2
            east_border_sit  = (sit[line,col+1] + sit[line,col])/2

            north_border_sic = (sic[line+1,col] + sic[line,col])/2
            south_border_sic = (sic[line-1,col] + sic[line,col])/2
            west_border_sic  = (sic[line,col-1] + sic[line,col])/2
            east_border_sic  = (sic[line,col+1] + sic[line,col])/2
             
            from_north_transport = from_north_si_drift * north_border_sit * north_border_sic * X_dist[line,col] #[m^3/s] 
            from_south_transport = from_south_si_drift * south_border_sit * south_border_sic * X_dist[line,col] #[m^3/s]
            from_east_transport  = from_east_si_drift  * east_border_sit  * east_border_sic  * Y_dist[line,col] #[m^3/s]
            from_west_transport  = from_west_si_drift  * west_border_sit  * west_border_sic  * Y_dist[line,col] #[m^3/s] 
                
            transport[line,col] = from_north_transport + from_south_transport + from_east_transport + from_west_transport #[m^3/s] of sea ice transport in the cell

    #################
    ### - BILAN - ###
    #################

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
    

    startMonthSIV = (previousSIT * previousSIC * surface + sit * sic * surface)/2
    endMonthSIV   = (nextSIT     * nextSIC     * surface + sit * sic * surface)/2

    bilan_siv = endMonthSIV - startMonthSIV #Positive when SIV increase over area

    MB = 30*24*60*60 * transport #[m^3/month]

            
    SIM = MB - bilan_siv

    return SIM

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


gate1_coord = np.array([[79,-19.5],[79,10.5]])
gate2_coord = np.array([[76,-19.5],[76,2.5]])
gate3_coord = np.array([[73,-20.5],[73,-3.5]])
gate4_coord = np.array([[70,-21.5],[70,-9.5]])

gate1 = [[gate1_coord[0,0],i] for i in np.linspace(gate1_coord[0,1],gate1_coord[1,1],10)]
gate2 = [[gate2_coord[0,0],i] for i in np.linspace(gate2_coord[0,1],gate2_coord[1,1],10)]
gate3 = [[gate3_coord[0,0],i] for i in np.linspace(gate3_coord[0,1],gate3_coord[1,1],10)]
gate4 = [[gate4_coord[0,0],i] for i in np.linspace(gate4_coord[0,1],gate4_coord[1,1],10)]

gate = [gate1,gate2,gate3,gate4]
MB_list  = [[],[],[],[]]
SIM_list = []
for coords,i in zip([gate1_coord,gate2_coord,gate3_coord,gate4_coord],range(4)):
    for file in os.listdir('ORAS5/Data/iicethic'):
        if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
            MB_list[i].append(MB(coords,file)*1e-9)
            if i == 0:
                SIM_list.append(SIM(file))
    #MB_list[i] = np.convolve(MB_list[i], np.ones(12)/12, mode = 'valid')

MB_list = np.array(MB_list)
SIM_list = np.array(SIM_list)

corr  = np.zeros((len(MB_list),len(SIM_list[0]),len(SIM_list[0][0])))
p_val = np.zeros((len(MB_list),len(SIM_list[0]),len(SIM_list[0][0])))

for i in range(len(MB_list)):
    for l in range(np.shape(SIM_list[0])[0]):
        for c in range(np.shape(SIM_list[0])[1]):
            corr[i,l,c]  = stats.pearsonr(np.convolve(MB_list[i],np.ones(12)/12, mode = 'valid'),np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid'))[0]
            p_val[i,l,c] = stats.pearsonr(np.convolve(MB_list[i],np.ones(12)/12, mode = 'valid'),np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid'))[1]


titles_list = ['',' ','','']

label_list = ['r','r','r','r']
v_list = [1,1,1,1]
titlesize = 30
labelsize = 22
fig,axes = plt.subplots( nrows = 2, ncols = 2,figsize = (16,16), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -9)})

plt.subplots_adjust(left=0.01,
                    bottom=0.05,  
                    top=0.93,
                    wspace = 0.17)

for map,p_value,i in zip(corr,p_val,range(4)):


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

    
    cs = axes.flatten()[i].pcolormesh(lon, lat,map,vmin = -v_list[i], vmax = v_list[i],cmap = "cmo.balance", transform=ccrs.PlateCarree())
    #contour = axes.flatten()[i].contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    
    cax = fig.add_axes([axes.flatten()[i].get_position().x1+0.005,axes.flatten()[i].get_position().y0,0.015,axes.flatten()[i].get_position().height])
    if i == 0:
        cb = plt.colorbar(cs, cax = cax,ticks = [-10,-5,0,5,10])
    else:
        cb = plt.colorbar(cs, cax = cax)

    cb.ax.tick_params(labelsize=labelsize)
    cb.set_label(label = label_list[i],size = labelsize) 

    coords_p_signi = []
    for k in range(len(p_value)):
        for l in range(len(p_value[0])):
            if p_value[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axes.flatten()[i].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())
    

    axes.flatten()[i].plot(np.array(gate[i])[:,1], np.array(gate[i])[:,0], linewidth = 3,color = 'purple',transform=ccrs.PlateCarree())
    
    axes.flatten()[i].text(np.array(gate[i])[-1,1]+3,np.array(gate[i])[0,0]+0.2,f'G{i+1}',color = 'purple',transform=ccrs.PlateCarree(),fontsize=25,weight='bold')
    axes.flatten()[i].set_title(titles_list[i],fontsize = titlesize)
    

plt.savefig(f'MT_plot/correlation/G1G2G3G4_MB_SIM.png')
plt.close() 