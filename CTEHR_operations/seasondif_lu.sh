#!/bin/bash
fluxdirectory='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'
outdirectory_first='/projects/0/ctdas/dkivits/DATA/Fluxes_CDO/2018'

for landuse_type in {0,1,2,5,8,11,14,17,18}
do

	droughtdirectory="$outdirectory_first/2018"
	outdir="$outdirectory_first/2018-2017_per_lu_$landuse_type"

	if [ ! -d $outdir ]; then
	  mkdir $outdir
	fi

	for month in {7..9}
	do
		cp $fluxdirectory/nep.20180${month}.nc $droughtdirectory/nep.20180${month}.nc
		cp $fluxdirectory/fire.20180${month}.nc $droughtdirectory/fire.20180${month}.nc
	done
	
	cdo -O mergetime -cat $droughtdirectory/nep.201807.nc $droughtdirectory/nep.201808.nc $droughtdirectory/nep.201809.nc $droughtdirectory/nep.merged.2018.nc
	cdo -O mergetime -cat $droughtdirectory/fire.201807.nc $droughtdirectory/fire.201808.nc $droughtdirectory/fire.201809.nc $droughtdirectory/fire.merged.2018.nc
		
	cdo -O mergetime -cat $outdir/nep.201807.nc $outdir/nep.201808.nc $outdir/nep.201809.nc $outdir/nep.merged.2018.mixed.nc
	cdo -O mergetime -cat $outdir/fire.201807.nc $outdir/fire.201808.nc $outdir/fire.201809.nc $outdir/fire.merged.2018.mixed.nc
	
	cdo -O add $outdir/nep.merged.2018.mixed.nc $outdir/fire.merged.2018.mixed.nc $outdir/combined.merged.2018.mixed.nc
	cdo -O add $droughtdirectory/nep.merged.2018.nc $droughtdirectory/fire.merged.2018.nc $droughtdirectory/combined.merged.2018.nc
	
	cdo -O timmean $outdir/combined.merged.2018.mixed.nc $outdir/average.JAS.2018.mixed.nc
	cdo -O timmean $droughtdirectory/combined.merged.2018.nc $droughtdirectory/average.JAS.2018.nc

	cdo -O sub $droughtdirectory/average.JAS.2018.nc $outdir/average.JAS.2018.mixed.nc $outdirectory_first/nepfire.dif.JAS.per_lu_${landuse_type}.2018_2017.nc
done
