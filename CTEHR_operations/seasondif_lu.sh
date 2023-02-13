#!/bin/bash
fluxfolder='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'
outfolder_first='/projects/0/ctdas/dkivits/DATA/fluxes'

#echo 'give landuse type (CORINE PFT code): '
#read landuse_type
for landuse_type in {0,1,2,5,8,11,14,17,18}
do

	droughtfolder="$outfolder_first/2018"
	outdir="$outfolder_first/2018-2017_per_lu_$landuse_type"

	if [ ! -d $outdir ]; then
	  mkdir $outdir
	fi

	for month in {7..9}
	do
		cp $fluxfolder/nep.20180${month}.nc $droughtfolder/nep.20180${month}.nc
		cp $fluxfolder/fire.20180${month}.nc $droughtfolder/fire.20180${month}.nc
	done
	
	cdo -O mergetime -cat $droughtfolder/nep.201807.nc $droughtfolder/nep.201808.nc $droughtfolder/nep.201809.nc $droughtfolder/nep.merged.2018.nc
	cdo -O mergetime -cat $droughtfolder/fire.201807.nc $droughtfolder/fire.201808.nc $droughtfolder/fire.201809.nc $droughtfolder/fire.merged.2018.nc
		
	cdo -O mergetime -cat $outdir/nep.201807.nc $outdir/nep.201808.nc $outdir/nep.201809.nc $outdir/nep.merged.2018.mixed.nc
	cdo -O mergetime -cat $outdir/fire.201807.nc $outdir/fire.201808.nc $outdir/fire.201809.nc $outdir/fire.merged.2018.mixed.nc
	
	cdo -O add $outdir/nep.merged.2018.mixed.nc $outdir/fire.merged.2018.mixed.nc $outdir/combined.merged.2018.mixed.nc
	cdo -O add $droughtfolder/nep.merged.2018.nc $droughtfolder/fire.merged.2018.nc $droughtfolder/combined.merged.2018.nc
	
	cdo -O timmean $outdir/combined.merged.2018.mixed.nc $outdir/average.JAS.2018.mixed.nc
	cdo -O timmean $droughtfolder/combined.merged.2018.nc $droughtfolder/average.JAS.2018.nc

	cdo -O sub $droughtfolder/average.JAS.2018.nc $outdir/average.JAS.2018.mixed.nc $outfolder_first/nepfire.dif.JAS.per_lu_${landuse_type}.2018_2017.nc
done
