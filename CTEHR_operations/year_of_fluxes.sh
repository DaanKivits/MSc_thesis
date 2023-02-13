#!/bin/bash
fluxfolder='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'

droughtfolder='/projects/0/ctdas/dkivits/DATA/fluxes/2018'

if [ ! -d $droughtfolder ]; then
  mkdir $droughtfolder
fi

mv "$droughtfolder"/* "$droughtfolder/backup"

for month in {1..12}
do
	if (($month < 10))
	then
                cp $fluxfolder/nep.20180$month.nc $droughtfolder/nep.20180$month.nc
                cp $fluxfolder/fire.20180$month.nc $droughtfolder/fire.20180$month.nc
	else
		cp $fluxfolder/nep.2018$month.nc $droughtfolder/nep.2018$month.nc
                cp $fluxfolder/fire.2018$month.nc $droughtfolder/fire.2018$month.nc
	fi
	
done

cdo mergetime $droughtfolder/nep.2018*.nc $droughtfolder/nep.merged.2018.nc
cdo mergetime $droughtfolder/fire.2018*.nc $droughtfolder/fire.merged.2018.nc
cdo add $droughtfolder/nep.merged.2018.nc $droughtfolder/fire.merged.2018.nc $droughtfolder/combined.merged.2018.nc
cdo yearmean $droughtfolder/combined.merged.2018.nc $droughtfolder/yearmean.2018.nc

