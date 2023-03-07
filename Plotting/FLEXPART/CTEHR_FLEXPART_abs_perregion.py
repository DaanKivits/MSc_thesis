# Daan Kivits, 2023

# A simple script that is used to calculate the monthly RMSE values between the observations and 2018 CTE-HR 
# flux sets, and plot the absolute timeseries of the transported 2018 CTE-HR fluxes compared to the observations in a timeseries
# plot that is subdivided into the Northern, Temperate, and Mediterranean climate regions (see thesis work for the included
# stations in this analysis). We included a TM5 background run as a reference, which is a 7-day mean of the TM5 simulation results.
# A 7-day running mean is also calculated and shown for both the simulated as well as the observed mixing ratios.

# This variation of the script returns a figure with a 16:9 aspect ratio, but this can be changed in line 33.

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
from scipy import stats
import netCDF4 as nc
import sys
sys.path.insert(0, '/projects/0/ctdas/dkivits/scripts/Plotting/functions')
import plotting_functions

RMSElist_2017, RMSElist_2018, RMSElist_lu = ([] for i in range(3))

Obs2018 = read_hdf("/projects/0/ctdas/dkivits/DATA/FLEXPART/code/ObsPack2018drought_2022update_withbackground.hdf",key='observations')

plt.rcParams['axes.grid'] = True
fig, ax = plt.subplots(nrows=3, ncols=2,gridspec_kw={'width_ratios': [6, 1]}, figsize=(16,9), sharex='col')
plt.subplots_adjust(bottom=0.08,top=0.95,left=0.1,right=0.95,wspace=0.08, hspace=0.17)

ax[2,0].set_xlabel('Time (month in 2018)', size=10)
ax[1,0].set_ylabel('Atmospheric CO$_{2}$ concentration (ppm)', size = 10)

regionlist = ['North','Temperate','Mediterranean']

for region in range(0,len(Obs2018.name.unique())):
    Obsdata_2018 = Obs2018[(Obs2018.name == str(regionlist[region])) & (Obs2018.type == 'tall-tower')]

    Obsdata_2018['dif_2018'] = Obsdata_2018['mix_2018']-Obsdata_2018['obs']
    dif_2018 = Obsdata_2018[['time', 'dif_2018']]
    Obsdata_2018['dif_2017'] = Obsdata_2018['mix_2017']-Obsdata_2018['obs']
    dif_2017 = Obsdata_2018[['time', 'dif_2017']]
    Obsdf = Obsdata_2018[['time','obs','mix_2018','background']]

    df_permonth_2018 = Obsdata_2018.copy()
    df_permonth_2018.index = pd.DatetimeIndex(df_permonth_2018.time)
    df_permonth_2018["SE_mix_2017"] = ((df_permonth_2018["mix_2017"] - df_permonth_2018["obs"]) ** 2) 
    df_permonth_2018["SE_mix_2018"] = ((df_permonth_2018["mix_2018"] - df_permonth_2018["obs"]) ** 2)
    df_permonth_2018 = df_permonth_2018.groupby(by=[df_permonth_2018.index.month]).mean()
    
    RMSE_2017, RMSE_2018, RMSE_lu = ([] for i in range(3))
    
    for month in range (1,13):
        RMSE_2017.append((df_permonth_2018["SE_mix_2017"] ** 0.5)[month])
        RMSE_2018.append((df_permonth_2018["SE_mix_2018"] ** 0.5)[month])
 
    RMSElist_2017.append(RMSE_2017)
    RMSElist_2018.append(RMSE_2018)

    # FOR (PARTIAL) TEMPORAL INTEPOLATION PURPOSES: REINDEX
    data_2018 = Obsdata_2018.groupby('time').mean()
    data_2018 = data_2018.reindex(pd.date_range(start = "2018-01-01", end = "2019-01-01", freq = '1H'), fill_value=np.nan)
    data_2018.mix_2017 = data_2018.mix_2017.interpolate(method='linear',limit=100)
    data_2018.mix_2018 = data_2018.mix_2018.interpolate(method='linear',limit=100)
    data_2018.background = data_2018.background.interpolate(method='linear',limit=100)
    data_2018.obs = data_2018.obs.interpolate(method='linear',limit=100)
    data_2018['dif_2018_int'] = data_2018['mix_2018']-data_2018['obs']
    data_2018['dif_2017_int'] = data_2018['mix_2017']-data_2018['obs']
    data_2018['dif_runmean_2017'] = data_2018['dif_2017_int'].rolling(window=168).mean()
    data_2018['dif_runmean_2018'] = data_2018['dif_2018_int'].rolling(window=168).mean()
    data_2018['runmean_obs'] = data_2018['obs'].rolling(window=168).mean()
    data_2018['runmean_mix_2018'] = data_2018['mix_2018'].rolling(window=168).mean()
    data_2018['runmean_background'] = data_2018['background'].rolling(window=168).mean()

    obs = ax[region,0].scatter(Obsdf['time'], Obsdf['obs'], label = 'observations', color = 'green', alpha = 0.2, s = 5, zorder = 1)
    ctehr2018 = ax[region,0].scatter(Obsdf['time'], Obsdf['mix_2018'], label = '2018 CTE-HR', color = 'orange', alpha = 0.2, s = 5, zorder = 1)
    #backgr = ax[region,0].plot(Obsdf['time'], Obsdf['background'], label = 'background TM5')
    obsline = ax[region,0].plot(data_2018.index, data_2018['runmean_obs'], label = 'observations runmean', color = 'green', alpha = 1, linewidth = 3, zorder = 5)
    ctehr2018line = ax[region,0].plot(data_2018.index, data_2018['runmean_mix_2018'], label = '2018 CTE-HR runmean', color = 'orange', alpha = 1, linewidth = 3, zorder = 5) 
    backgrline = ax[region,0].plot(data_2018.index, data_2018['runmean_background'], label = 'TM5 background', color = 'grey', alpha = 1, linewidth = 3, zorder = 5)
    
    count = dif_2018.describe().loc['count'][0].astype(int)
    mean2018 = dif_2018.describe()['dif_2018'].loc['mean'].round(2)
    mean2017 = dif_2017.describe()['dif_2017'].loc['mean'].round(2)
    std2017 = dif_2017.describe()['dif_2017'].loc['std'].round(2)
    std2018 = dif_2018.describe()['dif_2018'].loc['std'].round(2)

    # Make sure grid is below graph
    ax[region,0].set_axisbelow(True)
    ax[region,1].set_axisbelow(True)

    ax[region,0].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax[region,0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    
    # Set the ticks and ticklabels for all axes
    plt.setp(ax[region,0].xaxis.get_majorticklabels(), rotation=45, ha='right',rotation_mode='anchor')

    ax[region,0].set_xlim([datetime(2018,4,1), datetime(2018,10,1)])
    ax[region,0].set_ylim([380, 440])
    ax[region,1].set_ylim([-20, 20])
    
    ax[region,0].yaxis.set_ticks(np.arange(380, 440, 10))
    ax[region,1].yaxis.set_ticks(np.arange(-20, 25, 10))
    
    ax[region,0].set_title(regionlist[region])

    bin_size = int(1)
    min_edge = int(-15)
    max_edge = int(15)
    N = (max_edge-min_edge)/bin_size
    Nplus1 = int(N) + 1
    bins = np.linspace(min_edge, max_edge, Nplus1)
    
    difhist_2018 = sns.histplot(data=Obsdata_2018['dif_2018'], y=Obsdata_2018['dif_2018'], bins = bins, ax=ax[region,1], kde = True, color = 'orange', alpha=0.8)
    ax[region,1].set(ylabel='')
    
    if region == 0:
        # Put a legend below current axis
        leg = ax[region,0].legend(fancybox=True, prop={'size': 9}, borderaxespad = 0.75, framealpha=1, ncol = 3)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        ax[region,1].legend(['2018 CTE-HR'], fancybox=True, prop={'size': 9}, borderaxespad = 0.75, framealpha=1)
    
    # Include basic statistics!
    ax[region,1].text(0.95, 0.11, f"{mean2018} \u00B1 {std2018}", color = 'orange', horizontalalignment='right', verticalalignment='center', transform=ax[region,1].transAxes)
    ax[region,1].text(0.95, 0.04, f"N = {count}", color = 'k', horizontalalignment='right', verticalalignment='center', transform=ax[region,1].transAxes)
    
    # Include RMSE per month
    for months in range(0,6):
        ax[region,0].text(0.05 + ((1/6) * months), 0.07, "RMSE = " + str(round(RMSE_2018[months+3],2)), color = 'orange', verticalalignment = 'center', transform=ax[region,0].transAxes)    
    
    fignumlist = ['a','b','c']
    ax[region,1].text(1.05, 0.5, fignumlist[region], size = 12, fontweight = 'bold', verticalalignment = 'center', transform=ax[region,1].transAxes)

plt.savefig("/gpfs/work1/0/ctdas/dkivits/figures/CTEHR-FLEXPART/2018drought_CTEHR_FLEXPART_perregion.pdf",dpi=300,bbox_inches='tight')
plt.savefig("/gpfs/work1/0/ctdas/dkivits/figures/CTEHR-FLEXPART/2018drought_CTEHR_FLEXPART_perregion.png",dpi=300,bbox_inches='tight')
plt.show()

df = pd.DataFrame(RMSElist_2017)
df.to_csv("/projects/0/ctdas/dkivits/figures/CTEHR-FLEXPART/RMSElist_2017.csv", index=False)
df = pd.DataFrame(RMSElist_2018)
df.to_csv("/projects/0/ctdas/dkivits/figures/CTEHR-FLEXPART/RMSElist_2018.csv", index=False)