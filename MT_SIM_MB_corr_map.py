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

            from_north_si_drift = np.nan_to_num(from_north_si_drift)
            from_west_si_drift  = np.nan_to_num(from_west_si_drift )
            from_east_si_drift  = np.nan_to_num(from_east_si_drift )
            from_south_si_drift = np.nan_to_num(from_south_si_drift)
             
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
    

    startMonthSIV = (previousSIT * previousSIC * surface + sit * sic * surface)/2
    endMonthSIV   = (nextSIT     * nextSIC     * surface + sit * sic * surface)/2

    bilan_siv = endMonthSIV - startMonthSIV #Positive when SIV increase over area

    MB = 30*24*60*60 * transport #[m^3/month]

            
    SIM = MB - bilan_siv

    return SIM

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

def SIM_A1(file):
    iicevelv = np.loadtxt(f'ORAS5/Data/iicevelv/{file}')
    iicevelu = np.loadtxt(f'ORAS5/Data/iicevelu/{file}')
    sit      = np.loadtxt(f'ORAS5/Data/iicethic/{file}')
    sic      = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')



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

    SIM = MB - bilanSIV

    return SIM

def SIM_A2(file):
    iicevelv = np.loadtxt(f'ORAS5/Data/iicevelv/{file}')
    iicevelu = np.loadtxt(f'ORAS5/Data/iicevelu/{file}')
    sit      = np.loadtxt(f'ORAS5/Data/iicethic/{file}')
    sic      = np.loadtxt(f'ORAS5/Data/ileadfra/{file}')



    fromNorth = 0
    for lon_ in range(-19,-12):
        sit_interp = (np.nan_to_num(float(sit[     (lat == 73.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == 72.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == 73.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == 72.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == 73.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == 72.5) & (lon == lon_)])))/2

        fromNorth -= sit_interp * sic_interp * v_interp * float(X_dist[(lat == 73.5) & (lon == lon_)])
   
    fromSouth = 0
    for lon_ in range(-19,-12):
        sit_interp = (np.nan_to_num(float(sit[     (lat == 68.5) & (lon == lon_)])) + np.nan_to_num(float(sit[     (lat == 67.5) & (lon == lon_)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == 68.5) & (lon == lon_)])) + np.nan_to_num(float(sic[     (lat == 67.5) & (lon == lon_)])))/2
        v_interp   = (np.nan_to_num(float(iicevelv[(lat == 68.5) & (lon == lon_)])) + np.nan_to_num(float(iicevelv[(lat == 67.5) & (lon == lon_)])))/2

        fromSouth += sit_interp * sic_interp * v_interp * float(X_dist[(lat == 67.5) & (lon == lon_)])

    fromWest = 0
    for lat_ in [68.5,69.5,70.5,71.5,72.5]:
        sit_interp = (np.nan_to_num(float(sit[     (lat == lat_) & (lon == -19)])) + np.nan_to_num(float(sit[     (lat == lat_) & (lon == -20)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == lat_) & (lon == -19)])) + np.nan_to_num(float(sic[     (lat == lat_) & (lon == -20)])))/2
        v_interp   = (np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == -19)])) + np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == -20)])))/2

        fromWest += sit_interp * sic_interp * v_interp * float(Y_dist[(lat == lat_) & (lon == lon_)])

    fromEast = 0
    for lat_ in [68.5,69.5,70.5,71.5,72.5]:
        sit_interp = (np.nan_to_num(float(sit[     (lat == lat_) & (lon == -12)])) + np.nan_to_num(float(sit[     (lat == lat_) & (lon == -13)])))/2
        sic_interp = (np.nan_to_num(float(sic[     (lat == lat_) & (lon == -12)])) + np.nan_to_num(float(sic[     (lat == lat_) & (lon == -13)])))/2
        v_interp   = (np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == -12)])) + np.nan_to_num(float(iicevelu[(lat == lat_) & (lon == -13)])))/2

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
    
    previousSIT = np.where((lat >= 68.5) & (lat <= 72.5) & (lon <= -12) & (lon >= -19),previousSIT,0)
    previousSIC = np.where((lat >= 68.5) & (lat <= 72.5) & (lon <= -12) & (lon >= -19),previousSIC,0)
    nextSIT     = np.where((lat >= 68.5) & (lat <= 72.5) & (lon <= -12) & (lon >= -19),nextSIT,    0)
    nextSIC     = np.where((lat >= 68.5) & (lat <= 72.5) & (lon <= -12) & (lon >= -19),nextSIC,    0)


    startMonthSIV = (np.nansum(previousSIT*previousSIC*surface) + np.nansum(np.where((lat >= 68.5) & (lat <= 72.5) & (lon <= -12) & (lon >= -19),sit*sic*surface,0)))/2
    endMonthSIV   = (np.nansum(nextSIT*nextSIC*surface)         + np.nansum(np.where((lat >= 68.5) & (lat <= 72.5) & (lon <= -12) & (lon >= -19),sit*sic*surface,0)))/2

    bilanSIV = endMonthSIV - startMonthSIV #Positive when SIV increase over area

    SIM = MB - bilanSIV

    return SIM

MBG1_list     = []
MBG2_list     = []
SIM_list      = []
sohtc300_list = []

day_btwn_2011_start_and_2018_end = []
for file in os.listdir('ORAS5\Data\iicethic'):
    if int(file[:4]) >= 2011:
        print(file)
        a,b = MB_G1G2(file)
        MBG1_list.append(a)
        MBG2_list.append(b)
        SIM_list.append(SIM(file))
        sohtc300_list.append(np.where((lat > 70) & (lat < 80) & (lon > -30) & (lat > 65.4 + (76.5-65.4)/(9+17) * (lon + 17)) & (lon < 7),np.loadtxt(f'ORAS5/Data/sohtc300/{file}'),np.nan))
        day_btwn_2011_start_and_2018_end.append((date(int(file[:4]),int(file[5:7]),1) - date(2011,1,1)).days)
        plot(sohtc300_list[-1])
    """ if int(file[:4]) == 2013:
        break """

        
MBG1_list   = np.nan_to_num(np.array(MBG1_list))
MBG2_list   = np.nan_to_num(np.array(MBG2_list))
SIM_list = np.nan_to_num(np.array(SIM_list))
sohtc300_list = np.nan_to_num(np.array(sohtc300_list))
np.where((lon > 70) & (lon < 80) & (lat > -20) & (lat > 65.4 + (76.5-65.4)/(9+17) * (lon + 17)) & (lat < 5),np.array(sohtc300_list),np.nan)
corr_MB  = np.zeros(np.shape(SIM_list[0]))
p_val_MB = np.zeros(np.shape(SIM_list[0]))

corr_sohtc300  = np.zeros(np.shape(SIM_list[0]))
p_val_sohtc300 = np.zeros(np.shape(SIM_list[0]))

mean_sohtc300_GS = np.nanmean()
for l in range(np.shape(corr_MB)[0]):
    for c in range(np.shape(corr_MB)[1]):
        corr_MB[l,c]  = stats.pearsonr(np.convolve(MBG1_list,np.ones(12)/12, mode = 'valid'),np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid'))[0]
        p_val_MB[l,c] = stats.pearsonr(np.convolve(MBG1_list,np.ones(12)/12, mode = 'valid'),np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid'))[1]

        corr_sohtc300[l,c]  = stats.pearsonr(np.convolve(sohtc300_list[:,l,c],np.ones(12)/12, mode = 'valid'),np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid'))[0]
        p_val_sohtc300[l,c] = stats.pearsonr(np.convolve(sohtc300_list[:,l,c],np.ones(12)/12, mode = 'valid'),np.convolve(SIM_list[:,l,c],np.ones(12)/12, mode = 'valid'))[1]

        
gate1 = np.concatenate(([[82,i] for i in range(-12,18)],[[81,17],[80,17]]))

######################
#### - PLOTTING - ####
######################

labelsize = 20
titlesize = 27
ticklabelsize = 20
legendsize = 20
    
fig = plt.figure(figsize = (21,10))
# add grid specifications
gs = fig.add_gridspec(2,2)
axs_map = fig.add_subplot(gs[:,0], projection = ccrs.LambertConformal(central_longitude = -3))
axs_plot = []
for i in range(2):
    axs_plot.append(fig.add_subplot(gs[i,1]))
plt.subplots_adjust(top = 0.95, bottom = 0.1, hspace = 0.3, left = 0.01,wspace=0.1,right =0.95)

##############################
#### - plotting the map - ####
##############################
xlim = [-25, 20]
ylim = [65, 82.5]
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

cs = axs_map.pcolormesh(lon, lat, corr_MB,vmin = -1, vmax = 1,cmap = "cmo.balance", transform=ccrs.PlateCarree())
#### - P_VALUE - ####
coords_p_signi = []
for k in range(len(p_val_MB)):
    for l in range(len(p_val_MB[0])):
        if p_val_MB[k,l] <= 0.05:
            coords_p_signi.append([lon[k,l],lat[k,l]])
coords_p_signi = np.array(coords_p_signi)
if len(coords_p_signi) != 0:
    cs_non_signi = axs_map.scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())

#### - DISPLAY GATE - ####
axs_map.plot(np.array(gate1)[:,1], np.array(gate1)[:,0], linewidth = 4,color = 'black',transform=ccrs.PlateCarree())

axs_map.text(-16,81.7,'G',transform=ccrs.PlateCarree(),fontsize=20,weight='bold')


cax = fig.add_axes([axs_map.get_position().x0,axs_map.get_position().y0 - 0.02,axs_map.get_position().width,0.023])
cb = plt.colorbar(cs, cax = cax,orientation = 'horizontal')
cb.ax.tick_params(labelsize=15)
cb.set_label(label = "r",size = 15)

#################################
########### - Plots - ###########
#################################
MBG1 = []
MBG2 = []
MBG3 = []

SIMA1 = []
SIMA2 = []

for file in os.listdir('ORAS5/Data/iicethic'):
    if int(file[:4]) >= 2011 and int(file[:4]) <= 2018:
        a,b = MB_G1G2(file)
        MBG1.append(a*1e-9)
        MBG2.append(b*1e-9)
        SIMA1.append(SIM_A1(file))
        SIMA2.append(SIM_A2(file))

MBG1 = np.convolve(MBG1,np.ones(12)/12, mode = 'valid')
MBG2 = np.convolve(MBG2,np.ones(12)/12, mode = 'valid')

SIMA1 = np.convolve(SIMA1,np.ones(12)/12, mode = 'valid')
SIMA2 = np.convolve(SIMA2,np.ones(12)/12, mode = 'valid')

MBG1linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],MBG1)
MBG2linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],MBG2)

SIMA1linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],SIMA1)
SIMA2linregress = stats.linregress(day_btwn_2011_start_and_2018_end[11:],SIMA2)


axs_plot[0].plot(np.linspace(2012,2019,len(SIMA1)),SIMA1*1e-9,linewidth = 2,label = 'A2', color = 'blue')

ax_MB = plt.twinx(axs_plot[0])
ax_MB.plot(np.linspace(2012,2019,len(MBG1)),MBG1,linewidth = 2,label = 'A2', color = 'red')
#ax_sal.plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),SA2linregress.intercept + SA2linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'red')

r = stats.pearsonr(MBG1,SIMA1)[0]
p = stats.pearsonr(MBG1,SIMA1)[1]
axs_plot[0].set_title(f'A1, r: {round(r,2)}, p : {round(p,3)}',fontsize = titlesize)
#axs_plot[0].legend(fontsize = legendsize)
axs_plot[0].grid()
axs_plot[0].set_ylabel(r'$[km^3*month^{-1}]$',size = labelsize, color = 'blue')
axs_plot[0].tick_params(labelsize = ticklabelsize, color = 'blue')
axs_plot[0].tick_params(axis='y', colors='blue')
ax_MB.set_ylabel(r'$[km^{3}*month^{-1}]$',size = labelsize, color = 'red')
ax_MB.tick_params(labelsize = ticklabelsize, color = 'red')
ax_MB.tick_params(axis='y', colors='red')
axs_plot[0].set_xticklabels([])

#############################

axs_plot[1].plot(np.linspace(2012,2019,len(SIMA2)),SIMA2*1e-9,linewidth = 2,label = 'A2', color = 'blue')

ax_MB = plt.twinx(axs_plot[1])
ax_MB.plot(np.linspace(2012,2019,len(MBG1)),MBG1,linewidth = 2,label = 'A2', color = 'red')
#ax_sal.plot(np.linspace(2012,2019,len(day_btwn_2011_start_and_2018_end[11:])),SA2linregress.intercept + SA2linregress.slope * np.array(day_btwn_2011_start_and_2018_end[11:]), color = 'red')

r = stats.pearsonr(MBG1,SIMA2)[0]
p = stats.pearsonr(MBG1,SIMA2)[1]
axs_plot[1].set_title(f'A2, r: {round(r,2)}, p : {round(p,3)}',fontsize = titlesize)
#axs_plot[1].legend(fontsize = legendsize)
axs_plot[1].grid()
axs_plot[1].set_ylabel(r'$[km^3*month^{-1}]$',size = labelsize, color = 'blue')
axs_plot[1].tick_params(labelsize = ticklabelsize, color = 'blue')
axs_plot[1].tick_params(axis='y', colors='blue')
ax_MB.set_ylabel(r'$[km^{3}*month^{-1}]$',size = labelsize, color = 'red')
ax_MB.tick_params(labelsize = ticklabelsize, color = 'red')
ax_MB.tick_params(axis='y', colors='red')
axs_plot[1].set_xticklabels([])


plt.savefig('MT_plot/SIM_MB_corr.png')
plt.close()



















######################
#### - PLOTTING - ####
######################

labelsize = 20
titlesize = 27
ticklabelsize = 20
legendsize = 20
    
fig,axs_map = plt.subplots(ncols = 2, figsize = (13,8), subplot_kw={'projection':ccrs.LambertConformal(central_longitude = -5)})
plt.subplots_adjust(top = 0.95, bottom = 0.1, hspace = 0.3, left = 0.01,wspace=0.1,right =0.95)

##############################
#### - plotting the map - ####
##############################
xlim = [-25, 20]
ylim = [65, 82.5]
lower_space = 3 
rect = mpath.Path([[xlim[0], ylim[0]],
                [xlim[1], ylim[0]],
                [xlim[1], ylim[1]],
                [xlim[0], ylim[1]],
                [xlim[0], ylim[0]],
                ]).interpolated(20)
proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axs_map[0]) - axs_map[0].transData
rect_in_target = proj_to_data.transform_path(rect)
axs_map[0].set_boundary(rect_in_target)
axs_map[0].set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

axs_map[0].coastlines()
axs_map[0].gridlines()        

cs = axs_map[0].pcolormesh(lon, lat, corr_MB,vmin = -1, vmax = 1,cmap = "cmo.balance", transform=ccrs.PlateCarree())
#### - P_VALUE - ####
coords_p_signi = []
for k in range(len(p_val_MB)):
    for l in range(len(p_val_MB[0])):
        if p_val_MB[k,l] <= 0.05:
            coords_p_signi.append([lon[k,l],lat[k,l]])
coords_p_signi = np.array(coords_p_signi)
if len(coords_p_signi) != 0:
    cs_non_signi = axs_map[0].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())

#### - DISPLAY GATE - ####
axs_map[0].plot(np.array(gate1)[:,1], np.array(gate1)[:,0], linewidth = 4,color = 'black',transform=ccrs.PlateCarree())

axs_map[0].text(-16,81.7,'G',transform=ccrs.PlateCarree(),fontsize=20,weight='bold')


#########################################
########## - CORR WITH THERM - ##########
#########################################

proj_to_data = ccrs.PlateCarree()._as_mpl_transform(axs_map[1]) - axs_map[1].transData
rect_in_target = proj_to_data.transform_path(rect)
axs_map[1].set_boundary(rect_in_target)
axs_map[1].set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

axs_map[1].coastlines()
axs_map[1].gridlines()        

cs = axs_map[1].pcolormesh(lon, lat, corr_sohtc300,vmin = -1, vmax = 1,cmap = "cmo.balance", transform=ccrs.PlateCarree())
#### - P_VALUE - ####
coords_p_signi = []
for k in range(len(p_val_sohtc300)):
    for l in range(len(p_val_sohtc300[1])):
        if p_val_sohtc300[k,l] <= 0.05:
            coords_p_signi.append([lon[k,l],lat[k,l]])
coords_p_signi = np.array(coords_p_signi)
if len(coords_p_signi) != 0:
    cs_non_signi = axs_map[1].scatter(coords_p_signi[:,0], coords_p_signi[:,1],marker = '.', color = 'black',transform=ccrs.PlateCarree())




cax = fig.add_axes([axs_map[-1].get_position().x1,axs_map[-1].get_position().y0,0.02,axs_map[-1].get_position().height])
cb = plt.colorbar(cs, cax = cax)
cb.ax.tick_params(labelsize=15)
cb.set_label(label = "r",size = 15)


plt.savefig('MT_plot/SIM_MB_TH_corr.png')
plt.close()

