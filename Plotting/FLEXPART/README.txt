# Daan Kivits, 2023

The 'FLEXPART' directory contains Python scripts that were used in the various FLEXPART analyses. All the scripts are therefore
used to plot results in the CONCENTRATION space, not the FLUX space.

Here is a short explanation of the functionality of the scripts included in this directory:

CTEHR_FLEXPART_abs_perregion.py: A script that plots the simulated mixing ratios of the 2018 CTE-HR fluxes transported by FLEXPART
    with a TM5 background simulation as a reference, divided into the Northern, Temperate, and Mediterannean climate regions.
    The script plots both a 7-day running mean as well as all the available observations at each timestamp. This could be multiple 
    observations, as multiple stations are included in each climate region.
CTEHR_FLEXPART_abs_stationspecific.py: Similar to 'CTEHR_FLEXPART_abs_perregion.py', but this script plots the simulated mixing
    ratios for each of the measurement sites individually. Contrary to 'CTEHR_FLEXPART_abs_perregion.py', the script does not plot
    any running means, and thus only point measurements are plotted. Remember, for lowland measurement sites data is only available
    from 12 to 16 UTC, and for mountaineous / high-altitude measurement sites from 23 to 03 UTC.
CTEHR_FLEXPART_dif_perregion.py: Similar to 'CTEHR_FLEXPART_abs_perregion.py', but this script plots an anomaly timeseries in which 
    the simulation results of the transported 2018 CTE-HR fluxes are subtracted from the observed mixing ratios at a measurement site
    level. Moreover, this script also plots the results from a 2017 CTE-HR transport simulation, to compare to the 2018 CTE-HR
    simulation results.
CTEHR_FLEXPART_dif_perregion_perLU_bio.py: This script takes the flux sets in which the 2018 CTE-HR biospheric flux fields have been filled with
    2017 CTE-HR biospheric fluxes for each of the PFTs / land use types, and creates the same plots as the 'CTEHR_FLEXPART_dif_perregion.py' script.
CTEHR_FLEXPART_dif_perregion_perLU_fossil.py: This script does the same as 'CTEHR_FLEXPART_dif_perregion_perLU_bio.py', but now for flux sets in
    which the flux fields are filled with 2017 CTE-HR anthropogenic (fossil fuel) fluxes.
CTEHR_FLEXPART_dif_threestations.py: This script does the same as 'CTEHR_FLEXPART_dif_perregion.py', but instead of visualizing the three climate 
    regions it shows three seperate measurement sites that can be chosen by the user. 
RMSE_dif_barhplot_bio+fire_fossil.py: This script creates a side-by-side comparison of the change in RMSE that results from exchanging the biospheric 
    (bio+fire; left side of plot) and the anthropogenic (fossil fuel; right side of plot) fluxes with 2017 CTE-HR fluxes for each of the pixels that 
    belong to a certain PFT.
        NOTE: in order for this script to run we need the results from the 'CTEHR_FLEXPART_dif_perregion_perLU_bio' and 'CTEHR_FLEXPART_dif_perregion_
        perLU_fossil' scripts. 