from pandas import read_hdf, read_csv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd
import numpy as np
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU

def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float
            Document textwidth or columnwidth in pts
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1]) * 0.5

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim

stationcsv = read_csv('/projects/0/ctdas/dkivits/DATA/CTEHR_Observation_Sites.csv', sep=',', encoding = 'latin1')
#stationlist, stationnamelist = ([] for i in range(2))
#for stations in range(0, len(stationcsv)):
#    stationlist.append(stationcsv["ID"][stations].lower())
#    stationnamelist.append(stationcsv["station"][stations])

Obs = read_hdf("/projects/0/ctdas/dkivits/DATA/FLEXPART/code/Obsnew_newcode.hdf",key='observations')
stationnames = ['smr','ope','ohp']
stationnamelist = ['Hyytiälä (SMR, Finland)','Observatoire pérenne de l\'environnement (OPE, France)','Observatoire de Haute Provence (OHP, France)']

fig, ax = plt.subplots(nrows=3, ncols=1,sharex=True,sharey=True,figsize=set_size(483.69687,subplots=(3,1)))
axs = ax.ravel()
plt.subplots_adjust(bottom=0.16,top=0.95,left=0.1,right=0.95)

#plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
axs[2].set_xlabel('Time (month in 2018)')
axs[1].set_ylabel('Atmospheric CO$_{2}$\n concentration (ppm)')

# Remove ZEP and IZO from observations, as they lie outside the CTE-HR domain
Obs.drop(Obs[(Obs.code == 'ZEP')| (Obs.code == 'IZO')].index)

for i in range(3):
    mask = Obs.code == stationnames[i]
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

    ctehr = axs[i].plot(stationdata.index, stationdata['mix'], label = 'mixed CTE-HR + background',linestyle='none', marker='o', markersize = 1.5)
    obs = axs[i].plot(stationdata.index, stationdata['obs'], label = 'observed',linestyle='none', marker='o', markersize = 1.5)
    backgr = axs[i].plot(stationdata.index, stationdata['mix_background'], label = 'background TM5',linestyle='none', marker='o', markersize= 1.5)
    #ctehr = plt.scatter(stationdata.index, stationdata['mix'], label = 'CTE-HR + FLEXPART + background')
    #obs = plt.scatter(stationdata.index, stationdata['obs'], label = 'observed')
    #backgr = plt.scatter(stationdata.index, stationdata['mix_background'], label = 'background TM5')
    
    axs[i].xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=MO))
    axs[i].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    
    #plt.setp(ax[-1], xlabel='Time (months)')
    #plt.setp(ax[:], ylabel='Atmospheric CO$_{2}$ concentration (ppm)')
    
    # Set the ticks and ticklabels for all axes
    plt.setp(axs[i].xaxis.get_majorticklabels(), rotation=45, ha='right',rotation_mode='anchor')

    axs[i].set_xlim([datetime(2018,6,1), datetime(2018,9,1)])
    axs[i].set_ylim([380, 460])

    #if i == 2 or 3:
    # Shrink current axis's height by 10% on the bottom
    #box = axs[i].get_position()
    #axs[i].set_position([box.x0, box.y0 + box.height * 0.10,
    #        box.width, box.height * 0.9])

    if i == 0:
        # Put a legend below current axis
        axs[i].legend(fancybox=True, prop={'size': 6}, borderaxespad = 0.75, framealpha=1)
    
    #ax.legend(fancybox=True, shadow=True)
    axs[i].grid()
    axs[i].set_title(stationname, size=8)

plt.savefig("/gpfs/work1/0/ctdas/dkivits/figures/CTEHR-FLEXPART/2018drought_CTEHR_FLEXPART_" + str(stationnames[0]) + '_' + str(stationnames[1]) + '_' + str(stationnames[2]),dpi=300,bbox_inches='tight')
plt.show()

fig, ax = plt.subplots(figsize=set_size(483.69687))
plt.subplots_adjust(bottom=0.12,top=0.95,left=0.1,right=0.95)
dif = stationdata['mix']-stationdata['obs']

binwidth = 1
difscatter = ax.hist(dif, bins = (round(min(dif)), round(max(dif)) + binwidth, binwidth), 
        edgecolor = 'gray', color = 'gray')
plt.show()

