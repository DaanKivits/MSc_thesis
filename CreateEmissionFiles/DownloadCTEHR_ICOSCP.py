import netCDF4 as nc
from netCDF4 import stringtochar
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time
import glob
import xarray as xr
import os
from datetime import datetime, timedelta, time, date
import shutil
from dateutil.relativedelta import relativedelta
from icoscp.cpb.dobj import Dobj
from icoscp.collection import collection

CTEHR_linklist = ['https://meta.icos-cp.eu/collections/D1sJyzhK2lYg7X-Fthrn0iuk',
		'https://meta.icos-cp.eu/collections/FqzyT_IakfEPFGKFrhmVODLq',
		'https://meta.icos-cp.eu/collections/IzTvhhAXhrGTP5ui60FLiQey',
		'https://meta.icos-cp.eu/collections/HY8YHNa0uEWQCRC8iBpkGs04',
		'https://meta.icos-cp.eu/collections/xdLizVi0s1b79v8Y48qaz6G6',
		'https://meta.icos-cp.eu/collections/-QvhFaiazNYgREEPT1OQoIu4']

monthlydatalist = []
for month in CTEHR_linklist:
	try:
		monthlydata = collection.get(month)
	#object = [Dobj(datalink)]
	#monthlydatalist += object