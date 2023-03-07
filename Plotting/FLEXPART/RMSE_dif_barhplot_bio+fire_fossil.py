# Daan Kivits, 2023

# A simple plotting script used to plot the PFT-based 2017 flux injection experiment results, by means of visualisation
# of a monthly change in RMSE compared to a 'more correct' 2018 flux map. This way, we can see what effect the different
# PFTs have on the performance of the FLEXPART-based transportation of the CTE-HR fluxes. 

# In this script, show a side-by-side comparison of the change in RMSE that results from exchanging the biospheric 
# (bio+fire; left side of plot) and the anthropogenic (fossil fuel; right side of plot) fluxes with 2017 CTE-HR
# fluxes for each of the pixels that belong to a certain PFT.

# The CSV files that contain the simulation results are obtained by running the lu_specific plotting scripts and 
# manual re-structuring of the resulting CSV files.

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
from pandas import read_csv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
path = '/gpfs/work1/0/ctdas/dkivits/DATA/'

RMSE_fos = read_csv(path + 'RMSE_change_fossil.csv')
RMSE = read_csv(path + 'RMSE_change.csv')
con = read_csv(path + 'mean_contribution_change.csv')

RMSEfosdf = pd.DataFrame(RMSE_fos)
RMSEdf = pd.DataFrame(RMSE)
condf = pd.DataFrame(con)

RMSEfosdf.drop(columns=['Shrublands (Non-tundra)','Evergreen Broadleaf Forest ', 'C4 Crops'],inplace=True)
RMSEdf.drop(columns=['Shrublands (Non-tundra)','Evergreen Broadleaf Forest ', 'C4 Crops'], inplace=True)

mpl.rcParams.update({'axes.grid': True})
plt.rcParams['axes.axisbelow'] = True

fig, ax = plt.subplots(figsize = (10,7), nrows=1, ncols=2,gridspec_kw={'width_ratios': [1, 1]}, sharey = True, sharex = True) 
plt.subplots_adjust(bottom=0.1,top=0.8,left=0.1,right=0.98,wspace=0.05, hspace=0.12)

ax[0] = RMSEdf.iloc[3:10].plot.barh(legend = True, ax = ax[0], width = 0.85)
ax[1] = RMSEfosdf.iloc[3:10].plot.barh(legend=False, ax = ax[1], width = 0.85)

ax[0].set_xlabel('Change in RMSE (ppm)')
ax[1].set_xlabel('Change in RMSE (ppm)')
#ax[1].set_xlabel('Change in mean contribution to CTE-HR flux signal (%)')
ax[0].set_ylabel('Time (month in 2018)')

monthlist = []
for month in RMSEdf.iloc[3:10,0]:
    monthlist.append(str(month))

ax[0].set_xticks(np.arange(-0.3,0.4,0.1))
#ax[1].set_xticks(np.arange(-3,4,1))
ax[1].set_xticks(np.arange(-0.3,0.4,0.1))


ax[0].set_xlim([-0.35,0.35])
ax[1].set_xlim([-0.35,0.35])
#ax[1].set_xlim([-3,3])
ax[0].set_yticks(np.arange(0,7,1))
ax[1].set_yticks(np.arange(0,7,1))
ax[0].set_yticklabels(monthlist)
ax[1].set_yticklabels(monthlist)

ax[0].grid(axis = 'y')
ax[1].grid(axis = 'y')

ax[0].text(0.95, 0.95, 'a', size = 12, fontweight = 'bold', verticalalignment = 'center', transform=ax[0].transAxes)
ax[0].text(0.28, 1.14, '2017 biosphere fluxes', size = 11, fontweight = 'bold', verticalalignment = 'center', transform=ax[0].transAxes)
ax[1].text(0.95, 0.95, 'b', size = 12, fontweight = 'bold', verticalalignment = 'center', transform=ax[1].transAxes)
ax[1].text(0.28, 1.14, '2017 fossil fluxes', size = 11, fontweight = 'bold', verticalalignment = 'center', transform=ax[1].transAxes)

ax[0].legend(bbox_to_anchor=(0., 1.02, 2.04, .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0)
plt.savefig('/gpfs/work1/0/ctdas/dkivits/figures/RMSE_bio+fire_fossil_barhplot.png',dpi=300, bbox_inches= 'tight')
plt.savefig('/gpfs/work1/0/ctdas/dkivits/figures/RMSE_bio+fire_fossil_barhplot.pdf',dpi=300, bbox_inches= 'tight')
plt.show()
