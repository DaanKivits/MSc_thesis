# Daan Kivits, 2023

# A simple plotting script used to visually compare the FLEXPART simulation results for three measurement sites, by means of 
# plotting the anomalies between the observed and the simulated atmospheric mixing ratios. The three measurement sites can be
# chosen by the user by changing the station-specific codes in the variable 'stationnames'. Note that we did not perform
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
import seaborn as sns
from scipy.stats import norm
import sys
sys.path.insert(0, '/projects/0/ctdas/dkivits/scripts/Plotting/functions')
import plotting_functions

stationcsv = read_csv('/projects/0/ctdas/dkivits/DATA/CTEHR_Observation_Sites.csv', sep=',', encoding = 'latin1')
#stationlist, stationnamelist = ([] for i in range(2))
#for stations in range(0, len(stationcsv)):
#    stationlist.append(stationcsv["ID"][stations].lower())
#    stationnamelist.append(stationcsv["station"][stations])

Obs = read_hdf("/projects/0/ctdas/dkivits/DATA/FLEXPART/code/ObsPack2018drought_2022update_withbackground.hdf",key='observations')
stationnames = ['smr','ope','ohp']
stationnamelist = ['Hyytiälä (SMR, Finland)','Observatoire pérenne de l\'environnement (OPE, France)','Observatoire de Haute Provence (OHP, France)']

plt.rcParams['axes.grid'] = True
fig, ax = plt.subplots(nrows=3, ncols=2,gridspec_kw={'width_ratios': [4.5, 1]}, figsize=set_size(483.69687,subplots=(3,2)), sharex='col')
plt.subplots_adjust(bottom=0.08,top=0.95,left=0.1,right=0.95,wspace=0.05, hspace=0.12)

ax[2,0].set_xlabel('Time (month in 2018)')
ax[1,0].set_ylabel('Atmospheric CO$_{2}$\n concentration (ppm)')

# Remove ZEP and IZO from observations, as they lie outside the CTE-HR domain
Obs.drop(Obs[(Obs.code == 'ZEP')| (Obs.code == 'IZO')].index)

for i in range(3):
    mask = Obs.code == stationnames[i]
    stationdata = Obs[mask]
    stationdata.loc[:, 'mix_background'] = (stationdata.background)
    stationheight = str(int(stationdata.loc[:,'height'].values[0]))
    stationname = (stationnamelist[i] + ' - ' + stationheight + 'm')
       
    stationdata['dif_2018'] = stationdata['mix_2018']-stationdata['obs']
    stationdata['dif_2017'] = stationdata['mix_2017']-stationdata['obs']
     
    # FOR (PARTIAL) TEMPORAL INTEPOLATION PURPOSES: REINDEX
    stationdata.index = pd.DatetimeIndex(stationdata.time)
    stationdata = stationdata.reindex(pd.date_range(start = "2018-01-01", end = "2019-01-01", freq = '1H'), fill_value=np.nan)
    stationdata.mix_background = stationdata.mix_background.interpolate(method='linear',limit=19)
    stationdata.mix_2017 = stationdata.mix_2017.interpolate(method='linear',limit=100)
    stationdata.mix_2018 = stationdata.mix_2018.interpolate(method='linear',limit=100)
    stationdata.obs = stationdata.obs.interpolate(method='linear',limit=100)
    
    stationdata['dif_2018_int'] = stationdata['mix_2018']-stationdata['obs']
    stationdata['dif_2017_int'] = stationdata['mix_2017']-stationdata['obs']

    stationdata['dif_runmean_2017'] = stationdata['dif_2017_int'].rolling(window=168).mean()
    stationdata['dif_runmean_2018'] = stationdata['dif_2018_int'].rolling(window=168).mean()

    difplot_2017 = ax[i,0].scatter(stationdata.index, stationdata['dif_2017'], label = 'CTE-HR 2017 + FLEXPART model residual', color = 'blue', alpha = 0.8, s = 5)
    difplot_2018 = ax[i,0].scatter(stationdata.index, stationdata['dif_2018'], label = 'CTE-HR 2018 + FLEXPART model residual', color = 'red', alpha = 0.8, s = 5)
    difline_2017 = ax[i,0].plot(stationdata.index, stationdata['dif_runmean_2017'], color = 'blue', alpha = 1, linewidth = 2)
    difline_2018 = ax[i,0].plot(stationdata.index, stationdata['dif_runmean_2018'], color = 'red', alpha = 1, linewidth = 2)

    ax[i,0].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax[i,0].xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    
    plt.setp(ax[-1], xlabel='Time (months)')
    plt.setp(ax[:], ylabel='Atmospheric CO$_{2}$ concentration (ppm)')
    
    # Set the ticks and ticklabels for all axes
    plt.setp(ax[i,0].xaxis.get_majorticklabels(), rotation=45, ha='right',rotation_mode='anchor')

    ax[i,0].set_xlim([datetime(2018,4,1), datetime(2018,10,1)])
    
    ax[i,0].set_ylim([-20, 20])
    ax[i,1].set_ylim([-20, 20])
    
    start, end = [-20, 25]
    ax[i,0].yaxis.set_ticks(np.arange(start, end, 5))
    ax[i,1].yaxis.set_ticks(np.arange(start, end, 10))
    
    bin_size = int(1)
    min_edge = int(-15)
    max_edge = int(15)
    N = (max_edge-min_edge)/bin_size
    Nplus1 = int(N) + 1
    bins = np.linspace(min_edge, max_edge, Nplus1)
    
    difhist_2017 = sns.histplot(data=stationdata['dif_2017'], y=stationdata['dif_2017'], bins = bins, ax=ax[i,1], kde = True, color = 'blue', alpha=0.8)
    difhist_2018 = sns.histplot(data=stationdata['dif_2018'], y=stationdata['dif_2018'], bins = bins, ax=ax[i,1], kde = True, color = 'red', alpha=0.8)
    
    if i == 0:
        # Put a legend below current axis
        ax[i,0].legend(fancybox=True, prop={'size': 8}, borderaxespad = 0.75, framealpha=1)
        ax[i,1].legend(['dif_2017','dif_2018'], fancybox=True, prop={'size': 8}, borderaxespad = 0.75, framealpha=1)

plt.savefig("/gpfs/work1/0/ctdas/dkivits/figures/CTEHR-FLEXPART/2018drought_CTEHR_FLEXPART_DIF" + str(stationnames[0]) + '_' + str(stationnames[1]) + '_' + str(stationnames[2]),dpi=300,bbox_inches='tight')
plt.show()
