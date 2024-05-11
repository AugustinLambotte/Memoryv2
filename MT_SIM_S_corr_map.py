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


lon          = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat          = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")
surface      = np.loadtxt('ORAS5/Data/geo/Surface.txt')
X_dist       = np.loadtxt('ORAS5/Data/geo/X_dist.txt')
Y_dist       = np.loadtxt('ORAS5/Data/geo/Y_dist.txt')
earth_radius = 6370 * 1e3
EGC_path     = np.loadtxt('ORAS5/Data/geo/Current_path.txt')

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

def SIM_3months(file):
    """
        month of file have to be 1, 4, 7, 10
    """
    year  = int(file[:4])
    month = int(file[5:7])
    
    for month_ in range(month, month + 3):
        iicevelv = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicevelv/{year}-{month_:02d}.txt'))
        iicevelu = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicevelu/{year}-{month_:02d}.txt'))
        sit      = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{year}-{month_:02d}.txt'))
        sic      = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{year}-{month_:02d}.txt'))
        if month_ == month:
            transport = np.zeros(np.shape(iicevelv))
        # Net bilan of sea ice import in each cell
        
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
                    
                transport[line,col] += from_north_transport + from_south_transport + from_east_transport + from_west_transport #[m^3/s] of sea ice transport in the cell

    #################
    ### - BILAN - ###
    #################

    
    
    previousSIT = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{year}-{month:02d}.txt'))
    nextSIT     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/iicethic/{year}-{month+2:02d}.txt'))

    previousSIC = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{year}-{month:02d}.txt'))
    nextSIC     = np.nan_to_num(np.loadtxt(f'ORAS5/Data/ileadfra/{year}-{month+2:02d}.txt'))
    

    startMonthSIV = (previousSIT * previousSIC * surface + sit * sic * surface)/2
    endMonthSIV   = (nextSIT     * nextSIC     * surface + sit * sic * surface)/2

    bilan_siv = endMonthSIV - startMonthSIV #Positive when SIV increase over area

    MB = 3*30*24*60*60 * transport #[m^3/3months]

            
    SIM = MB - bilan_siv

    return SIM

def S(file,depth):
    Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:-3]}nc')
    Smap = np.array(Sds['vosaline'].sel(deptht = depth, method = 'nearest')[:])
    Smap[Smap<20] = np.nan
    Sds.close()
    return Smap
def S_3months(file,depth):
    """
        month of file have to be 1, 4, 7, 10
    """
    Smap = []
    month = int(file[5:7])
    for month_ in range(month, month + 3):
        Sds = xr.open_dataset(f'ORAS5/Data/vosaline/{file[:4]}-{month_:02d}.nc')
        Smap.append(np.array(Sds['vosaline'].sel(deptht = depth, method = 'nearest')[:]))
        Sds.close()

    Smap = np.array(Smap)
    Smap[Smap<20] = np.nan
    Smap = np.mean(Smap,axis = 0)
    return Smap
S_list   = []
SIM_list = []

day_btwn_2011_start_and_2018_end = []
for file in os.listdir('ORAS5\Data\iicethic'):
    if int(file[:4]) >= 2011:# and int(file[5:7]) in [1, 4, 7, 10]:
        print(file)
        S_list.append(S(file,0))
        SIM_list.append(SIM(file))
        day_btwn_2011_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)
    """ if int(file[:4]) == 2013:
        break """


S_list   = np.nan_to_num(np.array(S_list))
SIM_list = np.nan_to_num(np.array(SIM_list))
corr  = np.zeros(np.shape(S_list[0]))
p_val = np.zeros(np.shape(S_list[0]))

trend_sim = np.zeros(np.shape(S_list[0]))
p_sim     = np.zeros(np.shape(S_list[0]))
trend_sal = np.zeros(np.shape(S_list[0]))
p_sal     = np.zeros(np.shape(S_list[0]))
for l in range(np.shape(corr)[0]):
    for c in range(np.shape(corr)[1]):
        corr[l,c]  = stats.pearsonr(np.convolve(S_list[:,l,c],np.ones(12)/12, mode = 'valid'),np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid'))[0]
        p_val[l,c] = stats.pearsonr(np.convolve(S_list[:,l,c],np.ones(12)/12, mode = 'valid'),np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid'))[1]

        slope          = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sim[l,c]     = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sim[l,c] = slope * (365)

        slope          = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S_list[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_sal[l,c]     = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(S_list[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_sal[l,c] = slope * (365)

trend_sim[trend_sim == 0] = np.nan
trend_sal[trend_sal == 0] = np.nan
p_sim[p_sim == 1]         = np.nan
p_sal[p_sal == 1]         = np.nan
######################
#### - PLOTTING - ####
######################

labelsize = 21
titlesize = 30
ticklabelsize = 20
legendsize = 20
v = [1,0.4,0.1]
titles = ['','','','Correlation between sea ice melt \nand surface Salinity','Sea ice melt trend','Salinity trend']
cbar_labels = ['r',r'$[km^3*month^{-1}*y^{-1}]$',r'$[psu*y^{-1}]$']
fig,axs = plt.subplots(ncols = 2, figsize = (16,9), subplot_kw= {'projection' : ccrs.LambertConformal(central_longitude = -12)})
xlim = [-33, 12]
ylim = [65, 81]
plt.subplots_adjust(left=0.01,
                    bottom=0.05, 
                    top=0.93,)
lower_space = 3  
rect = mpath.Path([[xlim[0], ylim[0]],
                [xlim[1], ylim[0]],
                [xlim[1], ylim[1]],
                [xlim[0], ylim[1]],
                [xlim[0], ylim[0]],
                ]).interpolated(20)
for map,p_map,i in zip([corr,trend_sim*1e-9],[p_val,p_sim],range(2)):
#for map,p_map,i in zip([corr],[p_val],range(1)):
    proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axs.flatten()[i]) - axs.flatten()[i].transData
    rect_in_target = proj_to_data.transform_path(rect)
    axs.flatten()[i].set_boundary(rect_in_target)
    axs.flatten()[i].set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

    axs.flatten()[i].coastlines()
    axs.flatten()[i].gridlines()

    axs.flatten()[i].set_title(titles[i], fontsize = 21)
    #cs = axs.flatten()[i].contourf(lon, lat, map,levels = np.linspace(-v[i],v[i],100),extend = 'both',vmin = -v[i], vmax = v[i],cmap = "cmo.balance", transform=ccrs.PlateCarree())
    cs = axs.flatten()[i].pcolormesh(lon, lat, map,vmin = -v[i], vmax = v[i],cmap = "cmo.balance", transform=ccrs.PlateCarree())
    contour = axs.flatten()[i].contour(lon, lat, np.where(EGC_path == 1,1,0), levels = [0,1], colors = 'black', interpolation = 'none', transform=ccrs.PlateCarree())

    coords_p_signi = [] 
    for k in range(len(p_val)):
        for l in range(len(p_val[0])):
            if p_map[k,l] <= 0.05:
                coords_p_signi.append([lon[k,l],lat[k,l]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axs.flatten()[i].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())


    cax = fig.add_axes([axs.flatten()[i].get_position().x1,axs.flatten()[i].get_position().y0,0.015,axs.flatten()[i].get_position().height])
    cb = plt.colorbar(cs, cax = cax)
    cb.ax.tick_params(labelsize=15)
    cb.set_label(label = cbar_labels[i],size = labelsize)

plt.savefig('MT_plot/SIM_S_corr.png')
plt.close()