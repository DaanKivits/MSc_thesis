# Daan Kivits, 2023

# A simple script that is used to calculate the monthly RMSE values between the observations, FLEXPART-transported 2018 CTE-HR, 
# and 2018 - 2017 PFT-based anthropogenic (fossil fuel) injected flux sets. This is only done for the Temperate region.

# The script also compares these PFT-specific experiments to the observations in an anomalies timeseries plot that is 
# subdivided into the Northern, Temperate, and Mediterranean climate regions (see thesis work for the included stations 
# in this analysis). A 7-day running mean is also calculated and shown for both the simulated as well as the observed 
# mixing ratios.

# Lastly, the script calculates the contribution of the biospheric fluxes to the mixed flux signal for each of the simulated
# months.

# This script returns the CSV file that contains the monthly change in RMSE that is needed in the RMSE barplot script.

# This script returns a figure with a 12:9 aspect ratio.

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

luvar_nc = nc.Dataset('/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/nep.202208.nc','r', format='NETCDF3_CLASSIC')
fluxname = list(luvar_nc.variables.keys())[3]
array = luvar_nc.variables[fluxname]
lu = get_lu(array)

Obs2018 = read_hdf("/projects/0/ctdas/dkivits/DATA/FLEXPART/code/ObsPack2018drought_2022update_withbackground.hdf",key='observations') 

RMSElist_2017, RMSElist_2018, RMSElist_lu, RMSElist_2017temp, RMSElist_2018temp, RMSElist_lutemp, contributionlist_biosphere_2018temp, contributionlist_biosphere_lutemp = ([] for i in range(8))

for types in np.unique(lu):
    print("WORKING ON ... " + str(types))

    Obs = read_hdf("/projects/0/ctdas/dkivits/DATA/FLEXPART/code/landuse_experiments_fossil/ObsPack2018drought_2022update_withbackground_fossil_per_lu_" + str(types),key='observations')

    plt.rcParams['axes.grid'] = True
    fig, ax = plt.subplots(nrows=3, ncols=2,gridspec_kw={'width_ratios': [6, 1]}, figsize=(16,9), sharex='col')
    plt.subplots_adjust(bottom=0.08,top=0.95,left=0.1,right=0.95,wspace=0.08, hspace=0.17)

    ax[2,0].set_xlabel('Time (month in 2018)', size=10)
    ax[1,0].set_ylabel('Atmospheric CO$_{2}$ concentration (ppm)', size = 10)

    contribution = len(lu[lu == 1]) / (np.shape(lu)[0] * np.shape(lu)[1])

    regionlist = ['North','Temperate','Mediterranean']
    for region in range(0,len(Obs.name.unique())):
        Obsdata = Obs[(Obs.name == str(regionlist[region])) & (Obs.type == 'tall-tower')]
        Obsdata_2018 = Obs2018[(Obs2018.name == str(regionlist[region])) & (Obs2018.type == 'tall-tower')]

        Obsdata_2018['dif_2018'] = Obsdata_2018['mix_2018']-Obsdata_2018['obs']
        dif_2018 = Obsdata_2018[['time', 'dif_2018']]
        Obsdata_2018['dif_2017'] = Obsdata_2018['mix_2017']-Obsdata_2018['obs']
        dif_2017 = Obsdata_2018[['time', 'dif_2017']]
        Obsdata['mix'] = Obsdata['mix'] + Obsdata['background']
        Obsdata['dif_lu'] = Obsdata['mix'] - Obsdata['obs']
        dif_lu = Obsdata[['time','dif_lu']]

        df_permonth = Obsdata.copy()
        df_permonth_2018 = Obsdata_2018.copy()
        df_permonth_2018.index = pd.DatetimeIndex(df_permonth_2018.time)
        df_permonth.index = pd.DatetimeIndex(df_permonth.time)
        df_permonth["SE_mix_lu"] = ((df_permonth["mix"] - df_permonth["obs"]) ** 2)
        df_permonth["contribution_biosphere"] = abs(df_permonth["mix_biosphere"])/(abs(df_permonth["mix_biosphere"]) + abs(df_permonth["mix_fossil"]) + abs(df_permonth["mix_ocean"])+ abs(df_permonth["mix_fire"]))
        df_permonth_2018["SE_mix_2017"] = ((df_permonth_2018["mix_2017"] - df_permonth_2018["obs"]) ** 2) 
        df_permonth_2018["SE_mix_2018"] = ((df_permonth_2018["mix_2018"] - df_permonth_2018["obs"]) ** 2)
        df_permonth_2018["contribution_biosphere"] = abs(df_permonth_2018["mix_biosphere_2018"])/(abs(df_permonth_2018["mix_biosphere_2018"]) + abs(df_permonth_2018["mix_fossil_2018"]) + abs(df_permonth_2018["mix_ocean_2018"])+ abs(df_permonth_2018["mix_fire_2018"]))
        df_permonth = df_permonth.groupby(by=[df_permonth.index.month]).mean()
        df_permonth_2018 = df_permonth_2018.groupby(by=[df_permonth_2018.index.month]).mean()
        
        if region == 1:
            RMSE_2017temp, RMSE_2018temp, RMSE_lutemp, contribution_biosphere_lutemp, contribution_biosphere_2018temp = ([] for i in range(5))
            for month in range (1,13):
                RMSE_2017temp.append((df_permonth_2018["SE_mix_2017"] ** 0.5)[month])
                RMSE_2018temp.append((df_permonth_2018["SE_mix_2018"] ** 0.5)[month])
                RMSE_lutemp.append((df_permonth["SE_mix_lu"] ** 0.5)[month])
                contribution_biosphere_2018temp.append((df_permonth_2018["contribution_biosphere"])[month])
                contribution_biosphere_lutemp.append((df_permonth["contribution_biosphere"])[month])
        
            contributionlist_biosphere_2018temp.append(contribution_biosphere_2018temp)
            contributionlist_biosphere_lutemp.append(contribution_biosphere_lutemp)
            RMSElist_lutemp.append(RMSE_lutemp)
            RMSElist_2017temp.append(RMSE_2017temp)
            RMSElist_2018temp.append(RMSE_2018temp)
         
        # FOR (PARTIAL) TEMPORAL INTEPOLATION PURPOSES: REINDEX
        data = Obsdata.groupby('time').mean()
        data_2018 = Obsdata_2018.groupby('time').mean()
        data = data.reindex(pd.date_range(start = "2018-01-01", end = "2019-01-01", freq = '1H'), fill_value=np.nan)
        data_2018 = data_2018.reindex(pd.date_range(start = "2018-01-01", end = "2019-01-01", freq = '1H'), fill_value=np.nan)
        data.mix = data.mix.interpolate(method='linear',limit=100)
        data_2018.mix_2017 = data_2018.mix_2017.interpolate(method='linear',limit=100)
        data_2018.mix_2018 = data_2018.mix_2018.interpolate(method='linear',limit=100)
        data_2018.obs = data_2018.obs.interpolate(method='linear',limit=100)
        data.obs = data.obs.interpolate(method='linear',limit=100)
        data['dif_lu_int'] = data['mix']-data['obs']
        data_2018['dif_2018_int'] = data_2018['mix_2018']-data_2018['obs']
        data_2018['dif_2017_int'] = data_2018['mix_2017']-data_2018['obs']
        data['dif_runmean_lu'] = data['dif_lu_int'].rolling(window=168).mean()
        data_2018['dif_runmean_2017'] = data_2018['dif_2017_int'].rolling(window=168).mean()
        data_2018['dif_runmean_2018'] = data_2018['dif_2018_int'].rolling(window=168).mean()

        difplot_lu = ax[region,0].scatter(dif_lu['time'], dif_lu['dif_lu'], label = 'model residuals 2017-2018 mixed CTE-HR + FLEXPART', color = 'blue', alpha = 0.3, s = 5, zorder = 1)
        difplot_2018 = ax[region,0].scatter(dif_2018['time'], dif_2018['dif_2018'], label = 'model residuals 2018 CTE-HR + FLEXPART', color = 'orange', alpha = 0.3, s = 5, zorder = 1)
        difline_lu = ax[region,0].plot(data.index, data['dif_runmean_lu'], color = 'blue', alpha = 1, linewidth = 3, zorder = 5)
        difline_2018 = ax[region,0].plot(data_2018.index, data_2018['dif_runmean_2018'], color = 'orange', alpha = 1, linewidth = 3, zorder = 5)
        zeroline = ax[region,0].plot([datetime(2018,1,1),datetime(2019,1,1)], [0,0], color = 'k', linewidth = 1, zorder = 3)

        count = dif_2018.describe().loc['count'][0]
        mean2018 = dif_2018.describe()['dif_2018'].loc['mean'].round(2)
        std2018 = dif_2018.describe()['dif_2018'].loc['std'].round(2)
        meanlu = dif_lu.describe()['dif_lu'].loc['mean'].round(2)
        stdlu = dif_lu.describe()['dif_lu'].loc['std'].round(2)

        # Make sure grid is below graph
        ax[region,0].set_axisbelow(True)
        ax[region,1].set_axisbelow(True)

        ax[region,0].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
        ax[region,0].xaxis.set_major_formatter(mdates.DateFormatter('%m'))
        
        # Set the ticks and ticklabels for all axes
        plt.setp(ax[region,0].xaxis.get_majorticklabels(), rotation=45, ha='right',rotation_mode='anchor')

        ax[region,0].set_xlim([datetime(2018,4,1), datetime(2018,10,1)])
        ax[region,0].set_title(regionlist[region])
        ax[region,0].set_ylim([-25, 25])
        ax[region,1].set_ylim([-20, 20])
        
        start, end = [-25, 30]
        ax[region,0].yaxis.set_ticks(np.arange(start, end, 5))
        ax[region,1].yaxis.set_ticks(np.arange(-20, 25, 10))
        
        bin_size = int(1)
        min_edge = int(-15)
        max_edge = int(15)
        N = (max_edge-min_edge)/bin_size
        Nplus1 = int(N) + 1
        bins = np.linspace(min_edge, max_edge, Nplus1)
        
        difhist_lu = sns.histplot(data=Obsdata['dif_lu'], y=Obsdata['dif_lu'], bins = bins, ax=ax[region,1], kde = True, color = 'blue', alpha=0.8)
        difhist_2018 = sns.histplot(data=Obsdata_2018['dif_2018'], y=Obsdata_2018['dif_2018'], bins = bins, ax=ax[region,1], kde = True, color = 'orange', alpha=0.8)
        ax[region,1].set(ylabel='')
        
        if region == 0:
            # Put a legend below current axis
            ax[region,0].legend(fancybox=True, prop={'size': 9}, borderaxespad = 0.75, framealpha=1)
            ax[region,1].legend(['mixed fluxes','2018 fluxes'], fancybox=True, prop={'size': 9}, borderaxespad = 0.75, framealpha=1)
        
        # Include basic statistics!
        ax[region,1].text(0.95, 0.18, f"{meanlu} \u00B1 {stdlu}", color = 'blue', horizontalalignment='right', verticalalignment='center', transform=ax[region,1].transAxes)
        ax[region,1].text(0.95, 0.11, f"{mean2018} \u00B1 {std2018}", color = 'orange', horizontalalignment='right', verticalalignment='center', transform=ax[region,1].transAxes)
        ax[region,1].text(0.95, 0.04, f"N = {count}", color = 'k', horizontalalignment='right', verticalalignment='center', transform=ax[region,1].transAxes)
        
        # Include RMSE per month
        for months in range(0,6):
            ax[region,0].text(0.05 + ((1/6) * months), 0.12, "RMSE = " + str(round(RMSE_lu[months+3],2)), color = 'blue', verticalalignment = 'center', transform=ax[region,0].transAxes)
            ax[region,0].text(0.05 + ((1/6) * months), 0.07, "RMSE = " + str(round(RMSE_2018[months+3],2)), color = 'orange', verticalalignment = 'center', transform=ax[region,0].transAxes)
    
    fignumlist = ['a','b','c']
    ax[region,1].text(1.05, 0.5, fignumlist[region], size = 12, fontweight = 'bold', verticalalignment = 'center', transform=ax[region,1].transAxes)

    plt.savefig("/gpfs/work1/0/ctdas/dkivits/figures/CTEHR-FLEXPART/2018drought_CTEHR_FLEXPART_DIF_perregion_perlu_" + str(types) + "_fossil.pdf",dpi=300,bbox_inches='tight')
    plt.savefig("/gpfs/work1/0/ctdas/dkivits/figures/CTEHR-FLEXPART/2018drought_CTEHR_FLEXPART_DIF_perregion_perlu_" + str(types) + "_fossil.png",dpi=300,bbox_inches='tight')

df = pd.DataFrame(RMSElist_2017temp)
df.to_csv("/projects/0/ctdas/dkivits/figures/CTEHR-FLEXPART/RMSElist_2017.csv", index=False)
df = pd.DataFrame(RMSElist_2018temp)
df.to_csv("/projects/0/ctdas/dkivits/figures/CTEHR-FLEXPART/RMSElist_2018.csv", index=False)
df = pd.DataFrame(RMSElist_lutemp)
df.to_csv("/projects/0/ctdas/dkivits/figures/CTEHR-FLEXPART/RMSElist_lu.csv", index=False)
df = pd.DataFrame(contributionlist_biosphere_2018temp)
df.to_csv("/projects/0/ctdas/dkivits/figures/CTEHR-FLEXPART/contributionlist_biosphere_2018.csv", index=False)
df = pd.DataFrame(contributionlist_biosphere_lutemp)
df.to_csv("/projects/0/ctdas/dkivits/figures/CTEHR-FLEXPART/contributionlist_biosphere_lu.csv", index=False)

