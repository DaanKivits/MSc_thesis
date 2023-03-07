#!/bin/bash
fluxdirectory='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'
droughtdirectory='/projects/0/ctdas/dkivits/DATA/Fluxes_CDO/2018/2018'

if [ ! -d $droughtdirectory ]; then
  mkdir $droughtdirectory
fi

mv "$droughtdirectory"/* "$droughtdirectory/backup"

for month in {1..12}
do
	if (($month < 10))
	then
                cp $fluxdirectory/nep.20180$month.nc $droughtdirectory/nep.20180$month.nc
                cp $fluxdirectory/fire.20180$month.nc $droughtdirectory/fire.20180$month.nc
	else
		cp $fluxdirectory/nep.2018$month.nc $droughtdirectory/nep.2018$month.nc
                cp $fluxdirectory/fire.2018$month.nc $droughtdirectory/fire.2018$month.nc
	fi
	
done

cdo mergetime $droughtdirectory/nep.2018*.nc $droughtdirectory/nep.merged.2018.nc
cdo mergetime $droughtdirectory/fire.2018*.nc $droughtdirectory/fire.merged.2018.nc
cdo add $droughtdirectory/nep.merged.2018.nc $droughtdirectory/fire.merged.2018.nc $droughtdirectory/combined.merged.2018.nc
cdo yearmean $droughtdirectory/combined.merged.2018.nc $droughtdirectory/yearmean.2018.nc

