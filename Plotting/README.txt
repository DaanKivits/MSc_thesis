# Daan Kivits, 2023

The 'Plotting' directory contains multiple sub-directories with scripts that have been used to visualize the CTE-HR fluxes (flux space) and the 
FLEXPART and WRF-CHEM (concentration space) atmospheric transport simulation results. 

Here is a short explanation of the contents of the subdirectories:

- checkWRFruns: contains Python scripts that were used in to check the emission flux files provided by Friedemann Reum
(personal communications, 2022) and the CTE-HR fluxes. 
- CTE-HR_Fluxes: contains Python scripts that were used to spatially check the CTE-HR fluxes over Europe, to detect any
noticeable trends.
- FLEXPART: contains Python scripts that were used in the various FLEXPART analyses.
- Observations: contains a script that we used to create and visualize a data availability table of the atmospheric 
measurement sites that we used.
- functions: contains a Python function file with functions that were used to plot the results of our 
FLEXPART analyses, and a Python script that is used by the 'compareWRFoutput_visually.sh' Bash script.