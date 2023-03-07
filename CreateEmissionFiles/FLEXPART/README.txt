# Daan Kivits, 2023

The FLEXPART directory contains the scripts we used to transform the CTE-HR fluxes to the FLEXPART spatial grid. 

Here is a short explanation of the functionality of the scripts included in this directory:

- MakeFLEXPARTChemiFiles.py: This script creates FLEXPART-compatible flux files based on the spatial grid characteristics of
    the FLEXPART footprints, and fills these files with CTE-HR high-resolution fluxes for the biosphere (NEP and wildfires), ocean, 
    and anthropogenic sectors.
- MakeFLEXPARTChemiFiles_1x1.py: This script creates FLEXPART-compatible flux files with a 1 x 1 degree spatial resolution, and 
    fills these files with CTE-HR high-resolution fluxes for the biosphere (NEP and wildfires), ocean, and anthropogenic sectors.
- MakeFLEXPARTChemiFiles_perlu_bio.py: Same as 'MakeFLEXPARTChemiFiles.py', but this specific script creates flux files in which
     the 2018 CTE-HR biosphere flux fields (NEP and fire fields) are exchanged by the 2017 CTE-HR biosphere fluxes for each 
     of the specific landuse types.
- MakeFLEXPARTChemiFiles_perlu_fossil.py: Same as 'MakeFLEXPARTChemiFiles.py', but this specific script creates flux files in which
     the 2018 CTE-HR anthropogenic flux fields (fossil fuels) are exchanged by the 2017 CTE-HR anthropogenic fluxes for each 
     of the specific landuse types.
- submit_CreateFLEXPARTChemiFiles.py: This is a script that can be used to submit to the HPC cluster and let the FLEXPART emission flux file
    creation run remotely.

The directory also contains the 'notebooks' and 'function' sub-directories. The former contains a Jupyter notebook that I used to test the areas 
that were masked when selecting only CTE-HR gridcells that were classified under a certain PFT / landuse, and the latter contains a Python 
function file with functions that were used to extract the landuse data from the CORINE-based SiB4 PFT map and check the area covered by each of 
the landuse types.