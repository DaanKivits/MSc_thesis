#!/bin/bash
fluxdirectory='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'
outdirectory_first='/projects/0/ctdas/dkivits/DATA/Fluxes_CDO/2018'
droughtdirectory=$outdirectory_first/"2018"

if [ ! -d $droughdirectory ]; then
  mkdir $droughtdirectory
fi

for landuse_type in {0,1,2,5,8,11,14,17,18}
do

	outdir="$outdirectory_first/2018-2017_per_lu_$landuse_type"

	if [ ! -d $outdir ]; then
	  mkdir $outidr
	fi

	for month in {5..10}
	do	
		if (($month!=10))
		then
			cp $fluxdirectory/nep.20180${month}.nc $droughtdirectory/nep.20180${month}.nc
			cp $fluxdirectory/fire.20180${month}.nc $droughtdirectory/fire.20180${month}.nc

			cdo -O add $droughtdirectory/nep.20180$month.nc $droughtdirectory/fire.20180$month.nc $droughtdirectory/combined.20180$month.nc
			cdo -O add $outdir/nep.20180$month.nc $outdir/fire.20180$month.nc $outdir/combined.20180$month.mixed.nc
			cdo -O mergetime -cat $outdir/combined.*0$month.mixed.nc $outdir/merged.0$month.mixed.nc
			cdo -O yearmean $droughtdirectory/combined.20180$month.nc $droughtdirectory/combined.avg.20180$month.nc
			cdo -O timmean -cat $outdir/merged.0$month.mixed.nc $outdir/merged.avg.0$month.mixed.nc
			cdo -O sub $droughtdirectory/combined.avg.20180$month.nc $outdir/merged.avg.0$month.mixed.nc $outdirectory_first/nepfire.dif.per_lu_${landuse_type}.2018_2017.0$month.nc
	
		else
			cp $fluxdirectory/nep.2018${month}.nc $droughtdirectory/nep.2018${month}.nc
			cp $fluxdirectory/fire.2018${month}.nc $droughtdirectory/fire.2018${month}.nc

			cdo -O add $droughtdirectory/nep.2018$month.nc $droughtdirectory/fire.2018$month.nc $droughtdirectory/combined.2018$month.nc
			cdo -O add $outdir/nep.2018$month.nc $outdir/fire.2018$month.nc $outdir/combined.2018$month.mixed.nc
			cdo -O mergetime -cat $outdir/combined.*$month.mixed.nc $outdir/merged.$month.mixed.nc
			cdo -O yearmean $droughtdirectory/combined.2018$month.nc $droughtdirectory/combined.avg.2018$month.nc
			cdo -O timmean -cat $outdir/merged.$month.mixed.nc $outdir/merged.avg.$month.mixed.nc
			cdo -O sub $droughtdirectory/combined.avg.2018$month.nc $outdir/merged.avg.$month.mixed.nc $outdirectory_first/nepfire.dif.per_lu_${landuse_type}.2018_2017.$month.nc

		fi
	done
done
