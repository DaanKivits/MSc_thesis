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

stationcsv = read_csv('/projects/0/ctdas/dkivits/DATA/CTEHR_Observation_Sites.csv', sep=',', encoding = 'latin1')
Obs = read_hdf("/projects/0/ctdas/dkivits/DATA/FLEXPART/code/ObsPack2018drought_2022update_withbackground.hdf",key='observations')

timedf = pd.DataFrame({'time': (Obs.groupby(Obs["time"], as_index=False).mean().time)})
#lenlist = range(len(Obs.site.unique()),0)
lenlist = np.arange(42, -1, -1).tolist()
colorlist = []
for station in Obs.site.unique():
    stationdata = Obs[Obs.site == station]
    name = stationdata.name.unique()[0]
    if name == 'North':
        color = 'aqua'
    if name == 'Temperate':
        color = 'springgreen'
    if name == 'Mediterranean':
        color = 'y'
    colorlist.append(color)
    stationdf = pd.DataFrame({'time': stationdata.time, station:stationdata.obs})
    timedf = pd.merge(timedf, stationdf, how="left")
    timedf.loc[~np.isnan(timedf[station]), station] = lenlist[Obs.site.unique().tolist().index(station)]

#timedf['time'] = df['Date'].dt.strftime('%Y-%m-%d')
#dates = timedf['time'].date
timedf.set_index('time', inplace = True)

ax = timedf.plot(legend=False, figsize = (12,9), color = colorlist)
ax.set_yticks(range(len(Obs.site.unique())))

namelist = []
for names in Obs.site.unique():
    namelist.append(names.upper())

namelist.sort(reverse=True)
ax.set_yticklabels(namelist)

ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))

ax.set_xlim([datetime(2018,1,1), datetime(2019,1,1)])
ax.set_ylim([0.5, 42.5])

#ax.set_xticks(Obs.site.unique())
#ax.set_xticks(range(len(timedf['time'])))
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right',rotation_mode='anchor')
plt.tight_layout()

plt.savefig('/projects/0/ctdas/dkivits/figures/Data_availability_table.pdf', dpi = 300, bbox_inches = 'tight')
plt.savefig('/projects/0/ctdas/dkivits/figures/Data_availability_table.png', dpi = 300, bbox_inches = 'tight')
plt.show()

