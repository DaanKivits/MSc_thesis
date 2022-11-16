#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

- Computes xwavenum and ywavenum for WRF spectral nudging based on namelist.wps file
- Based on Gomez et al 2017 (10.1002/qj.3032) and others

Created on Mon Nov 11 09:53:11 2019

@author: friedemann
"""

import f90nml
import copy
import numpy as np


# Options
#namelist_filepath = "/Users/friedemann/projects/CHE/wrf_meteo_runs/preliminary_namelists/namelist.che_shanghai.wps"
#namelist_filepath = "/Users/friedemann/projects/CHE/wrf_meteo_runs/preliminary_namelists/namelist.che_beijing.wps"
#namelist_filepath = "/Users/friedemann/projects/CHE/wrf_meteo_runs/preliminary_namelists/namelist.che_berlin.wps"
#namelist_filepath = "/Users/friedemann/projects/SCARBO/WP17/simulations/namelist.domv4.wps"
#nml_kind = "wps"

namelist_filepath = "/Users/friedemann/wrf_run_ctdas/testrun12/namelist.input"
nml_kind = "input"


# Ivo
#namelist_filepath = "/Users/friedemann/projects/Ivo/nudging_namelist/namelist.input.txt"
#nml_kind = "input"


nudge_length_scale = 1000. # km (1000. = recommended by Gomez et al. 2017)


# Definitions
def parent_domains(namelist, nml_kind, domain):
    """
    Figure out the domain parents of this domain.
    An actually useful application of recursion!
    """

    if domain == 1:
        return []
    else:
        if nml_kind == "wps":
            parent = namelist["geogrid"]["parent_id"][domain-1] # "-1" because domain numbers start at 1
        elif nml_kind == "input":
            parent = namelist["domains"]["parent_id"][domain-1] # "-1" because domain numbers start at 1
        parents_of_parent = parent_domains(namelist, nml_kind, parent)
        parents_of_parent.append(parent)
        return parents_of_parent

# Read namelist
namelist = f90nml.read(namelist_filepath)

# Some definitions and initializations

if nml_kind == "wps":
    max_dom = namelist["share"]["max_dom"]
    dx = namelist["geogrid"]["dx"]/1000. # km
    dy = namelist["geogrid"]["dy"]/1000. # km
elif nml_kind == "input":
    max_dom = namelist["domains"]["max_dom"]
    dx = namelist["domains"]["dx"]
    dy = namelist["domains"]["dy"]
    if isinstance(dx, list):
        dx = dx[0]
    if isinstance(dy, list):
        dy = dy[0]
    dx /= 1000 # km
    dy /= 1000 # km


domain_extent = np.ndarray((max_dom, 2))

# Compute domain extent
for dom in range(1, max_dom + 1):
    # Compute resolution
      # Get all parents
    domain_list = parent_domains(namelist, nml_kind, dom)
    domain_list.append(dom)
    # Compute ratio of resolution with respect to mother of all domains
    if nml_kind == "wps":
        ratio_dom = np.array(namelist["geogrid"]["parent_grid_ratio"], ndmin=1)[np.array(domain_list)-1].prod()
    elif nml_kind == "input":
        ratio_dom = np.array(namelist["domains"]["parent_grid_ratio"], ndmin=1)[np.array(domain_list)-1].prod()
    dx_dom = dx/ratio_dom
    dy_dom = dy/ratio_dom
    # Compute and save domain extent
    if nml_kind == "wps":
        domain_extent[dom-1, 0] = dx_dom*(np.array(namelist["geogrid"]["e_we"], ndmin=1)[dom-1]-1)
        domain_extent[dom-1, 1] = dy_dom*(np.array(namelist["geogrid"]["e_sn"], ndmin=1)[dom-1]-1)
    elif nml_kind == "input":
        domain_extent[dom-1, 0] = dx_dom*(np.array(namelist["domains"]["e_we"], ndmin=1)[dom-1]-1)
        domain_extent[dom-1, 1] = dy_dom*(np.array(namelist["domains"]["e_sn"], ndmin=1)[dom-1]-1)


# Compute and print wave numbers
# Gomez 2017, Eq. (6)
xwavenum = np.round(domain_extent[:, 0]/nudge_length_scale + 1)
ywavenum = np.round(domain_extent[:, 1]/nudge_length_scale + 1)

print("xwavenum = ", xwavenum)
print("ywavenum = ", ywavenum)
print("Caution: Might want to set wavenums from 1 to 0?")


