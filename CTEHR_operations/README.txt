# Daan Kivits, 2023

The 'CTEHR_operations' directory contains scripts that were used to create the temporal average fluxmaps using a multitude of CDO commands and bash scripting for-loops.
Some points of attention:
- 2022 was not included in most of these analyses, since at the time of writing the scripts the CTE-HR growing season flux data was still incomplete.
- We used CDO commands here for their ease-of-use. Documentation can be found here: https://code.mpimet.mpg.de/projects/cdo/embedded/cdo_refcard.pdf
    Because the use of CDO commands proved to be very cumbersome when a hierarchical (HDF-like) and grouped file structure is neccessary, 
    we chose to move to simpler Python array and raster transformation operations later down the road.

Here is a short explanation of the functionality of the scripts included in this directory:

- monthdif.sh: (1) copies CTE-HR NEP and wildfire fluxes for the growing season of 2018 and for a set reference period (which could be multiple years);
(2) combines the NEP and wildfire fluxes into a single signal (the biosphere flux signal); (3) merges the monthly flux files into a multi-year monthly flux file
(for both the 2018 and non-2018 flux sets); (4) performs a yearly mean on these multiyear flux files to create a multi-year monthly average 
(for both the 2018 and non-2018 flux sets); (5) subtracts the 2018 average monthly flux fields with the climatology to create a monthly anomaly biospheric flux field.
- seasondif.sh: Similar to 'monthdif.sh', but creates a seasonal anomaly of 2018 versus the climatology.    
- monthdif_lu.sh: Performs the operation described by 'monthdif.sh' for every landuse type included in the CORINE-based SiB4 Plant Functional Type (PFT) map, 
which classifies the landuse over Europe into multiple categories.
- seasondif_lu.sh: Performs the operation described by 'seasondif.sh' for every landuse type included in the CORINE-based SiB4 Plant Functional Type (PFT) map, 
which classifies the landuse over Europe into multiple categories.
- year_of_fluxes.sh: A seperate script that creates a yearly mean biospheric flux field for the year of 2018.

