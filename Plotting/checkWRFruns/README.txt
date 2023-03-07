# Daan Kivits, 2023

The 'checkWRFruns' folder contains Python-based scripts that were used in to check both the emission flux files provided by Friedemann Reum
(personal communications, 2022) and the CTE-HR fluxes. 

Here is a short explanation of the functionality of the scripts included in this directory:

- check_CTEHRfluxset.py: A simple plotting script used to check the CTE-HR fluxes.
- check_interpolation_methods.py: A simple plotting script used to check the effect of different interpolation methods (bilinear
    and nearest neighbour interpolation) on the re-gridding interpolation of CTE-HR fluxes.
- check_originalfluxset.py: A simple plotting script used to check the flux set provided by Friedemann Reum 
    (personal communication, 2022) and the same flux set that has been filled with high emissions 
    (2 mol m-2 s-1 (or 7200 mol km-2 hr-1)).
- compareWRFoutput_visually.sh: The script creates a temporal average mixing ratio over the testrun period of WRF-CHEM results from
    simulations that use four different flux sets, and visualizes these results side-by-side to enable a visual comparison 
    between these different simulations setups.