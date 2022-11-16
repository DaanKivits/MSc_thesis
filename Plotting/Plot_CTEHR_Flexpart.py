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
 
#fig, ax = plt.subplots(nrows=2, ncols=2,sharex=True,sharey=True, figsize = (16, 9))
#axs = ax.ravel()
#plt.subplots_adjust(wspace=0.05, hspace=0.2)

# Remove ZEP and IZO from observations, as they lie outside the CTE-HR domain
Obs.drop(Obs[(Obs.code == 'ZEP')| (Obs.code == 'IZO')].index)

for i in range(0, len(stationlist)):
#for i in range(2):
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

    #ctehr = axs[i].plot(stationdata['time'], stationdata['mix'], label = 'mixed CTE-HR + background',linestyle='none', marker='o', markersize = 3)
    #obs = axs[i].plot(stationdata['time'], stationdata['obs'], label = 'observed',linestyle='none', marker='o', markersize = 3)
    #backgr = axs[i].plot(stationdata['time'], stationdata['mix_background'], label = 'background TM5',linestyle='none', marker='o', markersize=3)
    ctehr = plt.scatter(stationdata.index, stationdata['mix'], label = 'CTE-HR + FLEXPART + background')
    obs = plt.scatter(stationdata.index, stationdata['obs'], label = 'observed')
    backgr = plt.scatter(stationdata.index, stationdata['mix_background'], label = 'background TM5')
    
    ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=MO))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    
    plt.ylabel('Atmospheric CO$_{2}$ concentration (ppm)')
    plt.xlabel('Time (month in 2018)')
    
    #plt.setp(ax[-1, :], xlabel='Time (months)')
    #plt.setp(ax[:, 0], ylabel='Atmospheric CO$_{2}$ concentration (ppm)')

    # Set the ticks and ticklabels for all axes
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right',rotation_mode='anchor')

    ax.set_xlim([datetime(2018,6,1), datetime(2018,9,1)])
    ax.set_ylim([380, 460])

    #if i == 2 or 3:
    # Shrink current axis's height by 10% on the bottom
    #box = axs[i].get_position()
    #axs[i].set_position([box.x0, box.y0 + box.height * 0.10,
    #        box.width, box.height * 0.9])

    #if i == 1:
    #    # Put a legend below current axis
    #    axs[i].legend(fancybox=True, shadow=True)
    
    ax.legend(fancybox=True, shadow=True)
    ax.grid()
    ax.set_title(stationname)
    
    plt.savefig("/gpfs/work1/0/ctdas/dkivits/figures/CTEHR-FLEXPART/2018drought_CTEHR_FLEXPART_" + str(stationlist[i]),dpi=300,bbox_inches='tight')
    #plt.show()
