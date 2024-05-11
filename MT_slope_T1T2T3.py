import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import cartopy.crs as ccrs
import gsw
import seawater as sw
import os
import matplotlib.path as mpath


grid_transect1 = [[77.5,-8 ],[77.5,  0]]
grid_transect2 = [[76.5,-8 ],[76.5,-1 ]]
grid_transect3 = [[75.5,-13],[75.5,-6 ]]
grid_transect4 = [[69.5,-23],[69.5,-17]]
grid_transect5 = [[68.5,-25],[68.5,-18]]

grid_transect = [grid_transect1,grid_transect2,grid_transect3,grid_transect4,grid_transect5]
transectName  = ['Transect 1','Transect 2','Transect 3','Transect 4','Transect 5']
SSlope_list     = {}
TSlope_list     = {}
sigmaSlope_list = {}
 
depth_list = [0,30]

for d in depth_list:
    SSlope_list[f'{d}']     = {'Transect 1':[], 'Transect 2':[],'Transect 3':[],'Transect 4':[],'Transect 5':[]}
    TSlope_list[f'{d}']     = {'Transect 1':[], 'Transect 2':[],'Transect 3':[],'Transect 4':[],'Transect 5':[]}
    sigmaSlope_list[f'{d}'] = {'Transect 1':[], 'Transect 2':[],'Transect 3':[],'Transect 4':[],'Transect 5':[]}
m = 0

for file in os.listdir('ORAS5/Data/votemper'):
    m += 1
    
    print(file[:-3])
    dsS = xr.open_dataset(f"ORAS5/Data/vosaline/{file}")
    dsT = xr.open_dataset(f"ORAS5/Data/votemper/{file}")
    for transect_, transectName_ in zip(grid_transect,transectName):
        for i in range(len(depth_list)):
            
            S_west = dsS['vosaline'].sel(lon = transect_[0][1],deptht = depth_list[i], lat = transect_[0][0], method = 'nearest')
            S_east = dsS['vosaline'].sel(lon = transect_[1][1],deptht = depth_list[i], lat = transect_[1][0], method = 'nearest')

            T_west = dsT['votemper'].sel(lon = transect_[0][1],deptht = depth_list[i], lat = transect_[0][0], method = 'nearest')
            T_east = dsT['votemper'].sel(lon = transect_[1][1],deptht = depth_list[i], lat = transect_[1][0], method = 'nearest')

            sigma_west = sw.pden(S_west,T_west,sw.pres(depth_list[i],transect_[0][0]))-1000
            sigma_east = sw.pden(S_east,T_east,sw.pres(depth_list[i],transect_[1][0]))-1000

            dist = int(sw.extras.dist([transect_[0][0],transect_[1][0]],[transect_[0][1],transect_[1][1]])[0])
            
            SSlope_list[f'{depth_list[i]}'][f'{transectName_}'].append((S_east-S_west)/dist)
            TSlope_list[f'{depth_list[i]}'][f'{transectName_}'].append((T_east-T_west)/dist)
            sigmaSlope_list[f'{depth_list[i]}'][f'{transectName_}'].append((sigma_east-sigma_west)/dist)
    
    dsS.close()
    dsT.close()

color_list = ['red','blue','green','orange','purple', 'yellow']
for transectName_ in transectName:

    
    fig,axes = plt.subplots(nrows = 3, figsize = (10,11), sharex = True)
    labelsize = 18
    titlesize = 18
    ticklabelsize = 15
    legendsize = 10    

    fig.suptitle(f'{transectName_}',fontsize = 25)
    for j in range(len(depth_list)):
        depth_ = str(depth_list[j])
        ratioSalDist_ra = np.convolve(SSlope_list[depth_][transectName_], np.ones(12)/12, mode = 'valid')
        axes[0].plot([year for year in np.linspace(2011,2019,len(SSlope_list[depth_][transectName_]))], np.array(SSlope_list[depth_][transectName_])*1e3, alpha = 0.4,         color = color_list[j])
        axes[0].plot([year for year in np.linspace(2012,2019,len(ratioSalDist_ra))], np.array(ratioSalDist_ra)*1e3, label = f'depth = {depth_}m', linewidth = 3, color = color_list[j])
    axes[0].set_title(f'Salinity',fontsize = titlesize)
    axes[0].grid()
    axes[0].tick_params(labelsize = ticklabelsize)
    axes[0].set_ylabel(r'$10^{-3}psu/km$',size = labelsize)
    axes[0].legend(fontsize = legendsize)
    axes[0].set_xlim(2011,2019)
    #axes[0].set_ylim(-5,25)

    ratioTempDist = []
    for j in range(len(depth_list)):
        depth_ = str(depth_list[j])
        ratioTempDist_ra = np.convolve(TSlope_list[depth_][transectName_], np.ones(12)/12, mode = 'valid')
        axes[1].plot([year for year in np.linspace(2011,2019,len(TSlope_list[depth_][transectName_]))], np.array(TSlope_list[depth_][transectName_]), alpha = 0.4,      color = color_list[j])
        axes[1].plot([year for year in np.linspace(2012,2019,len(ratioTempDist_ra))], np.array(ratioTempDist_ra), label = f'depth = {depth_}m', linewidth = 3, color = color_list[j])
    axes[1].set_title(f'Potential temperature',fontsize = titlesize)
    axes[1].grid()
    axes[1].tick_params(labelsize = ticklabelsize)
    axes[1].set_ylabel(r'$°C/km$',size = labelsize)
    axes[1].legend(fontsize = legendsize)

    axes[1].set_xlim(2011,2019)
    #axes[1].set_ylim(-0.02,0.06)

    ratioSigmDist = []
    for j in range(len(depth_list)):
        depth_ = str(depth_list[j])
        ratioSigmDist_ra = np.convolve(sigmaSlope_list[depth_][transectName_], np.ones(12)/12, mode = 'valid')
        axes[2].plot([year for year in np.linspace(2011,2019,len(sigmaSlope_list[depth_][transectName_]))], np.array(sigmaSlope_list[depth_][transectName_]), alpha = 0.4,      color = color_list[j])
        axes[2].plot([year for year in np.linspace(2012,2019,len(ratioSigmDist_ra))], np.array(ratioSigmDist_ra), label = f'depth = {depth_}m', linewidth = 3, color = color_list[j])
    axes[2].set_title(r'$\sigma_{\theta}$',fontsize = titlesize)
    axes[2].grid()
    axes[2].tick_params(labelsize = ticklabelsize)
    axes[2].set_ylabel(r'$kg * m^{-3} * km^{-1}$',size = labelsize)
    axes[2].legend(fontsize = legendsize)
    axes[2].set_xlim(2011,2019)
    #axes[2].set_ylim(-0.0026,0.0205)
    plt.savefig(f'MT_plot/slope_{transectName_}.png')
    plt.close()
