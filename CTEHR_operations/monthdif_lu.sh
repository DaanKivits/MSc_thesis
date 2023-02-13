#!/bin/bash
fluxfolder='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'
outfolder_first='/projects/0/ctdas/dkivits/DATA/fluxes'
droughtfolder=$outfolder_first/"2018"

if [ ! -d $droughfolder ]; then
  mkdir $droughtfolder
fi

#echo 'give landuse type (CORINE PFT code): '
#read landuse_type

for landuse_type in {0,1,2,5,8,11,14,17,18}
do

	outdir="$outfolder_first/2018-2017_per_lu_$landuse_type"

	if [ ! -d $outdir ]; then
	  mkdir $outidr
	fi

	for month in {5..10}
	do	
		if (($month!=10))
		then
			cp $fluxfolder/nep.20180${month}.nc $droughtfolder/nep.20180${month}.nc
	                cp $fluxfolder/fire.20180${month}.nc $droughtfolder/fire.20180${month}.nc

			cdo -O add $droughtfolder/nep.20180$month.nc $droughtfolder/fire.20180$month.nc $droughtfolder/combined.20180$month.nc
			cdo -O add $outdir/nep.20180$month.nc $outdir/fire.20180$month.nc $outdir/combined.20180$month.mixed.nc
			cdo -O mergetime -cat $outdir/combined.*0$month.mixed.nc $outdir/merged.0$month.mixed.nc
			cdo -O yearmean $droughtfolder/combined.20180$month.nc $droughtfolder/combined.avg.20180$month.nc
			cdo -O timmean -cat $outdir/merged.0$month.mixed.nc $outdir/merged.avg.0$month.mixed.nc
			cdo -O sub $droughtfolder/combined.avg.20180$month.nc $outdir/merged.avg.0$month.mixed.nc $outfolder_first/nepfire.dif.per_lu_${landuse_type}.2018_2017.0$month.nc
	
		else
			cp $fluxfolder/nep.2018${month}.nc $droughtfolder/nep.2018${month}.nc
	                cp $fluxfolder/fire.2018${month}.nc $droughtfolder/fire.2018${month}.nc

			cdo -O add $droughtfolder/nep.2018$month.nc $droughtfolder/fire.2018$month.nc $droughtfolder/combined.2018$month.nc
			cdo -O add $outdir/nep.2018$month.nc $outdir/fire.2018$month.nc $outdir/combined.2018$month.mixed.nc
			cdo -O mergetime -cat $outdir/combined.*$month.mixed.nc $outdir/merged.$month.mixed.nc
			cdo -O yearmean $droughtfolder/combined.2018$month.nc $droughtfolder/combined.avg.2018$month.nc
			cdo -O timmean -cat $outdir/merged.$month.mixed.nc $outdir/merged.avg.$month.mixed.nc
			cdo -O sub $droughtfolder/combined.avg.2018$month.nc $outdir/merged.avg.$month.mixed.nc $outfolder_first/nepfire.dif.per_lu_${landuse_type}.2018_2017.$month.nc

		fi
	done
done
