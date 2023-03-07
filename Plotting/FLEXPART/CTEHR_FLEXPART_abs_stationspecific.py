# Daan Kivits, 2023

# A simple plotting script used to plot the FLEXPART simulation results for each station individually, by means of 
# plotting the background, and the observed and the simulated atmospheric mixing ratios. Note that we did not perform
# any interpolation on the data, and hence only data is available for daytime hours (12 to 16 UTC, local time differs
# per station) for lowland tower measurement sites, and nighttime hours high-altitude measurement sites (23 to 03 UTC).

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
from pandas import read_hdf, read_csv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd
import numpy as np
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU

stationcsv = read_csv('/projects/0/ctdas/dkivits/DATA/CTEHR_Observation_Sites.csv', sep=',')
stationlist, stationnamelist = ([] for i in range(2))
for stations in range(0, len(stationcsv)):
    stationlist.append(stationcsv["ID"][stations].lower())
    stationnamelist.append(stationcsv["station"][stations])

Obs = read_hdf("/projects/0/ctdas/dkivits/DATA/FLEXPART/code/Obsnew_newcode.hdf",key='observations')
#stationlist = ['htm','gat','tac','ope']
#stationnamelist = ['Hyltemossa (HTM, Sweden)','Gartow (GAT, Germany)','Tacolneston (TAC, UK)','Observatoire Pérenne de l\'Environnement (OPE, France)']

# Remove ZEP and IZO from observations, as they lie outside the CTE-HR domain
Obs.drop(Obs[(Obs.code == 'ZEP')| (Obs.code == 'IZO')].index)

for i in range(0, len(stationlist)):
    fig, ax = plt.subplots(figsize=(16,9))
    mask = Obs.code == str(stationlist[i])
    stationdata = Obs[mask]
    stationdata.loc[:, 'mix_background'] = (stationdata.background)
    stationheight = str(int(stationdata.loc[:,'height'].values[0]))
    stationname = (stationnamelist[i] + ' - ' + stationheight + 'm')
    
    # FOR (PARTIAL) TEMPORAL INTEPOLATION PURPOSES: REINDEX
    stationdata.index = pd.DatetimeIndex(stationdata.time)
    #stationdata = stationdata.reindex(pd.date_range(start = "2018-01-01", end = "2019-01-01", freq = '1H'), fill_value=np.nan)
    #stationdata.mix_background = stationdata.mix_background.interpolate(method='linear',limit=19)
    #stationdata.mix = stationdata.mix.interpolate(method='linear',limit=19)
    #stationdata.obs = stationdata.obs.interpolate(method='linear',limit=19)

    ctehr = plt.scatter(stationdata.index, stationdata['mix'], label = 'CTE-HR + FLEXPART + background')
    obs = plt.scatter(stationdata.index, stationdata['obs'], label = 'observed')
    backgr = plt.scatter(stationdata.index, stationdata['mix_background'], label = 'background TM5')
    
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    
    plt.ylabel('Atmospheric CO$_{2}$ concentration (ppm)')
    plt.xlabel('Time (month in 2018)')
    
    # Set the ticks and ticklabels for all axes
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right',rotation_mode='anchor')

    ax.set_xlim([datetime(2018,6,1), datetime(2018,9,1)])
    ax.set_ylim([380, 460])

    if i == 1:
        # Put a legend below current axis
        ax.legend(fancybox=True, shadow=True)
    
    ax.grid()
    ax.set_title(stationname)
    
    plt.savefig("/gpfs/work1/0/ctdas/dkivits/figures/CTEHR-FLEXPART/2018drought_CTEHR_FLEXPART_" + str(stationlist[i]),dpi=300,bbox_inches='tight')