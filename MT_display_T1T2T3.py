import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import cartopy.crs as ccrs
import gsw
import os
import matplotlib.path as mpath
from scipy import stats
from datetime import date
import seawater as sw



grid_transect1 = [[77.5,-8 ],[77.5,-7 ],[77.5,-6 ],[77.5,-5 ],[77.5,-4 ],[77.5,-3 ],[77.5,-2 ]]
grid_transect2 = [[76.5,-9 ],[76.5,-8 ],[76.5,-7 ],[76.5,-6 ],[76.5,-5 ],[76.5,-4 ],[76.5,-3 ]]
grid_transect3 = [[75.5,-14],[75.5,-13],[75.5,-12],[75.5,-11],[75.5,-10],[75.5,-9 ],[75.5,-8 ],[75.5,-7 ]]
grid_transect4 = [[69.5,-23],[69.5,-22],[69.5,-21],[69.5,-20],[69.5,-19],[69.5,-18]]
grid_transect5 = [[68.5,-26],[68.5,-25],[68.5,-24],[68.5,-23],[68.5,-22]]
transectName   = ['T1','T2','T3','T4','T5']
grid_transect = [grid_transect1,grid_transect2,grid_transect3,grid_transect4,grid_transect5]
Dens_grad_list = [[],[],[],[],[]]
day_btwn_2011_start_and_2018_end = []
for file in os.listdir('ORAS5/Data/votemper'):
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        dsS = xr.open_dataset(f"ORAS5/Data/vosaline/{file}")
        dsT = xr.open_dataset(f"ORAS5/Data/votemper/{file}")
        depth = np.array(dsS['deptht'])
        day_btwn_2011_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)
        for transect_, transectName_,i in zip(grid_transect,transectName,range(5)):
            SProfile = []
            TProfile = []
            sigmaProfile = []
            distance = []
            for point in range(len(transect_)):
                if point == 0:
                    distance.append(0) #km
                else:
                    distance.append(int(sw.extras.dist([transect_[0][0],transect_[point][0]],[transect_[0][1],transect_[point][1]])[0])) #km
                SProfile.append(np.array(dsS['vosaline'].sel(lon = transect_[point][1], lat = transect_[point][0], method = 'nearest')))
                TProfile.append(np.array(dsT['votemper'].sel(lon = transect_[point][1], lat = transect_[point][0], method = 'nearest')))
                sigmaProfile.append(sw.pden(SProfile[-1],TProfile[-1],sw.pres(depth,transect_[point][0]))-1000)
            sigmaProfile = np.array(sigmaProfile)
            Dens_grad = []
            for j in range(len(sigmaProfile[0])):
                Dens_grad.append(np.gradient(sigmaProfile[:,j],distance))
                
            Dens_grad_list[   i].append(np.array(Dens_grad))
            
        dsS.close()
        dsT.close()






lon = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/longitude.txt")
lat = np.loadtxt("C:/Users/Augustin/Desktop/Physique/Master_Climatologie/Memoire/Codev2/ORAS5/Data/geo/latitude.txt")
earth_radius = 6370 * 1e3
EGC_path = np.loadtxt('ORAS5/Data/geo/Current_path.txt')

ke              = []
day_btwn_2011_start_and_2018_end = []
for file in os.listdir('ORAS5/Data/ke'):
    if file == 'month_mean':
        continue
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        ke.append(np.loadtxt(f'ORAS5/Data/ke/{file}'))
        day_btwn_2011_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)

ke               = np.nan_to_num(np.array(ke))

ke[:,EGC_path != 1]               = np.nan

trend_ke          = np.zeros((np.shape(ke)[1:]))

p_ke          = np.zeros((np.shape(ke)[1:]))

for l in range(np.shape(trend_ke)[0]):
    for c in range(np.shape(trend_ke)[1]):
        slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(ke[:,l,c],np.ones(12)/12, mode = 'valid')).slope
        p_ke[l,c]              = stats.linregress(day_btwn_2011_start_and_2018_end[11:],np.convolve(ke[:,l,c],np.ones(12)/12, mode = 'valid')).pvalue
        trend_ke[l,c]          = slope * (365)


trend_ke[trend_ke == 0] = np.nan
p_ke[p_ke == 1]         = np.nan

month_name = ['January', 'February', 'March', 'April', 'May', 'June', 'July','August','September','October','November','December']


grid_transect1 = [[77.5,-8 ],[77.5,-7 ],[77.5,-6 ],[77.5,-5 ],[77.5,-4 ],[77.5,-3 ],[77.5,-2 ]]
grid_transect2 = [[76.5,-9 ],[76.5,-8 ],[76.5,-7 ],[76.5,-6 ],[76.5,-5 ],[76.5,-4 ],[76.5,-3 ]]
grid_transect3 = [[75.5,-14],[75.5,-13],[75.5,-12],[75.5,-11],[75.5,-10],[75.5,-9 ],[75.5,-8 ],[75.5,-7 ]]
grid_transect4 = [[69.5,-23],[69.5,-22],[69.5,-21],[69.5,-20],[69.5,-19],[69.5,-18]]
grid_transect5 = [[68.5,-26],[68.5,-25],[68.5,-24],[68.5,-23],[68.5,-22]]


labelsize = 30
titlesize = 18
ticklabelsize = 20
legendsize = 15
    
fig = plt.figure(figsize = (32,17))
# add grid specifications
gs = fig.add_gridspec(5,5)
axs_map = fig.add_subplot(gs[:,0:3], projection = ccrs.LambertConformal(central_longitude = -10))
axs_vert = []
for i in range(5):
    axs_vert.append(fig.add_subplot(gs[i,3:5]))
plt.subplots_adjust(top = 0.95, bottom = 0.1, hspace = 0.3, left = 0.01,wspace=0.3,right =0.9)

##############################
#### - plotting the map - ####
##############################
xlim = [-33, 12]
ylim = [65, 81]
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

cs = axs_map.pcolormesh(lon, lat, trend_ke*1e3,vmin = -1, vmax = 1,cmap = "cmo.balance", transform=ccrs.PlateCarree())
cax = fig.add_axes([axs_map.get_position().x1 - 0.44,axs_map.get_position().y0 - 0.02,axs_map.get_position().width*0.95,axs_map.get_position().height*0.05])
cb = plt.colorbar(cs, cax = cax, orientation = 'horizontal', ticks = [-1,-0.5,0,0.5,1])
cb.ax.tick_params(labelsize=35)
cb.set_label(label = r'$[mJ\cdot kg^{-1}\cdot y^{-1}]$',size = 35)
### - Transect - ###

axs_map.scatter(np.array(grid_transect1)[:,1], np.array(grid_transect1)[:,0],marker = 's', color = 'blue',transform=ccrs.PlateCarree())
axs_map.scatter(np.array(grid_transect2)[:,1], np.array(grid_transect2)[:,0],marker = 's', color = 'red',transform=ccrs.PlateCarree())
axs_map.scatter(np.array(grid_transect3)[:,1], np.array(grid_transect3)[:,0],marker = 's', color = 'purple',transform=ccrs.PlateCarree())
axs_map.scatter(np.array(grid_transect4)[:,1], np.array(grid_transect4)[:,0],marker = 's', color = 'yellow',transform=ccrs.PlateCarree())
axs_map.scatter(np.array(grid_transect5)[:,1], np.array(grid_transect5)[:,0],marker = 's', color = 'orange',transform=ccrs.PlateCarree())

axs_map.text(np.array(grid_transect1)[0,1]-3,np.array(grid_transect1)[0,0],'T1',transform=ccrs.PlateCarree(),fontsize=35,weight='bold')
axs_map.text(np.array(grid_transect2)[0,1]-3,np.array(grid_transect2)[0,0],'T2',transform=ccrs.PlateCarree(),fontsize=35,weight='bold')
axs_map.text(np.array(grid_transect3)[0,1]-3,np.array(grid_transect3)[0,0],'T3',transform=ccrs.PlateCarree(),fontsize=35,weight='bold')
axs_map.text(np.array(grid_transect4)[-1,1]+1,np.array(grid_transect4)[0,0],'T4',transform=ccrs.PlateCarree(),fontsize=35,weight='bold')
axs_map.text(np.array(grid_transect5)[-1,1]+1,np.array(grid_transect5)[0,0],'T5',transform=ccrs.PlateCarree(),fontsize=35,weight='bold')
################################
#### - Plotting dens grad - ####
################################
for trans,i in zip(Dens_grad_list,range(5)):
    trans = np.array(trans)
    p     = np.zeros(np.shape(trans[0]))
    trend = np.zeros(np.shape(trans[0]))
    for line in range(len(trans[0])):
        for col in range(len(trans[0][0])):
            averaged               = np.convolve(trans[:,line,col], np.ones(12)/12,mode = 'valid')
            slope                  = stats.linregress(day_btwn_2011_start_and_2018_end[11:],averaged).slope
            p[line,col]            = stats.linregress(day_btwn_2011_start_and_2018_end[11:],averaged).pvalue
            trend[line,col]        = slope * (365)

                
    
    axs_vert[i].set_title(f'{transectName[i]}',fontsize = 20,fontweight="bold")

   
    
    lonList = np.array(grid_transect[i])[:,1]
    cs = axs_vert[i].pcolormesh(lonList[:],depth,trend*1e3, cmap = 'cmo.balance', vmin = -1, vmax = 1)
    if i ==4:
        axs_vert[i].set_xlabel('Longitude [°]',size = labelsize)
    axs_vert[i].set_ylabel('Depth [m]',size = labelsize)
    axs_vert[i].tick_params(labelsize = ticklabelsize)
    axs_vert[i].set_ylim(200,0)
    axs_vert[i].grid()
    

    coords_p_signi = []
    for k in range(len(p)):
        for l in range(len(p[0])):
            if p[k,l] <= 0.05:
                coords_p_signi.append([lonList[l],depth[k]])
    coords_p_signi = np.array(coords_p_signi)
    if len(coords_p_signi) != 0:
        cs_non_signi = axs_vert[i].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black')


cax = fig.add_axes([axs_vert[-1].get_position().x1+0.02,axs_vert[-1].get_position().y0 + 0.05,0.02,axs_vert[-1].get_position().height * 5])
cb = plt.colorbar(cs, cax = cax,ticks = [-1,-0.5,0,0.5,1])
cb.ax.tick_params(labelsize=35)
cb.set_label(label = r'$[g\cdot m^{-3}\cdot km^{-1}\cdot y^{-1}]$',size = 35)    


plt.savefig(f'MT_plot/T1-T2-T3-T4-T5.png')
plt.close()
