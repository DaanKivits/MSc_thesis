# Daan Kivits, 2023

The 'WRF_CHEM' directory contains the scripts we used to test the fluxes that were provided by Friedemann Reum (2022, personal communications),
and transform the CTE-HR fluxes to the WRF-CHEM spatial grid. 

Here is a short explanation of the functionality of the scripts included in this directory:

- FillWRFChemiWithCTEHR.py: This script takes CTE-HR flux data and fills (already existing) WRFchemi emission files with this data for 
    every hourly timestep at a certain vertical injection level (that ranges from 0-8 in these example WRFchemi emission files).
- FillWRFChemiWithHighEmissions.py: This script fills (already existing) WRFchemi emission files with a constant high-emission 
    flux field of 2 mol m-2 s-1 (or 7200 mol km-2 hr-1). The rest of the fields are set to 0!
- FillWRFChemiWithZeros.py: This script fills (already existing) WRFchemi emission files with a constant flux field of 0.
- MakeWRFChemiFiles.py: This script creates WRFchemi flux files based on the grid definitions that results from running the WPS geogrid algorithm,
    and fills these files with CTE-HR high-resolution fluxes for the biosphere (NEP and wildfires), ocean, and anthropogenic sectors.
    The script therefore does not need any existing WRFChemi emission files as input.
- SjoerdBarten_examplescript.py: This script was provided by Sjoerd Barten (personal communications, 2022) and serves as a starting point
    for using the Python Basemap package to transform the CTE-HR flux fields, and visualizing the results.
- submit_CreateWRFChemiFiles.py: This is a script that can be used to submit to the HPC cluster and let the WRF-CHEM emission flux file
    creation run remotely.