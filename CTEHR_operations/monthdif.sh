#!/bin/bash

fluxdirectory='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'

echo 'give output directory name: '
read outdir

if [ ! -d "$outdir/backup" ]; then
  mkdir -p "$outdir/backup"
fi

droughtdirectory='/projects/0/ctdas/dkivits/DATA/Fluxes_CDO/2018'
fluxdir='/projects/0/ctdas/dkivits/DATA/Fluxes_CDO'

if [ ! -d "$droughtdirectory/backup" ]; then
  mkdir -p "$droughtdirectory"
fi

mv "$outdir"/* "$outdir/backup"
mv "$droughtdirectory"/* "$droughtdirectory/backup"

#years={2017..2021}
years=2017

for month in {5..10}
do	
	if (($month!=10))
        then
			cp $fluxdirectory/nep.20180$month.nc $droughtdirectory/nep.20180$month.nc
			cp $fluxdirectory/fire.20180$month.nc $droughtdirectory/fire.20180$month.nc
			cdo add $droughtdirectory/nep.20180$month.nc $droughtdirectory/fire.20180$month.nc $droughtdirectory/combined.20180$month.nc
		else
			cp $fluxdirectory/nep.2018$month.nc $droughtdirectory/nep.2018$month.nc
			cp $fluxdirectory/fire.2018$month.nc $droughtdirectory/fire.2018$month.nc
			cdo add $droughtdirectory/nep.2018$month.nc $droughtdirectory/fire.2018$month.nc $droughtdirectory/combined.2018$month.nc
	fi
done

#for year in {2017..2021}
for year in $years
do
	if (($year!=2018)) && (($year!=2022))
	then
		for month in {5..10}
		do
			if (($month!=10))
			then
				cp $fluxdirectory/nep.${year}0$month.nc $outdir/nep.${year}0$month.nc
				cp $fluxdirectory/fire.${year}0$month.nc $outdir/fire.${year}0$month.nc
				cdo add $outdir/nep.${year}0$month.nc $outdir/fire.${year}0$month.nc $outdir/combined.${year}0$month.nc
			else
				cp $fluxdirectory/nep.${year}$month.nc $outdir/nep.${year}$month.nc
				cp $fluxdirectory/fire.${year}$month.nc $outdir/fire.${year}$month.nc
				cdo add $outdir/nep.${year}$month.nc $outdir/fire.${year}$month.nc $outdir/combined.${year}$month.nc
			fi
		done
	fi

	if (($year==2022))
	then
		for month in {5..8}
		do
			cp $fluxdirectory/nep.${year}0$month.nc $outdir/nep.${year}0$month.nc
			cp $fluxdirectory/fire.${year}0$month.nc $outdir/fire.${year}0$month.nc
			cdo add $outdir/nep.${year}0$month.nc $outdir/fire.${year}0$month.nc $outdir/combined.${year}0$month.nc
		done
	fi
done

for month in {5..10}
do	
	if (($month!=10))
        then
			cdo mergetime -cat $outdir/combined.*0$month.nc $outdir/merged.0${month}.nc
			cdo yearmean $droughtdirectory/combined.20180$month.nc $droughtdirectory/combined.avg.20180$month.nc
			cdo timmean -cat $outdir/merged.0$month.nc $outdir/merged.avg.$years.0$month.nc
			cdo sub $droughtdirectory/combined.avg.20180$month.nc $outdir/merged.avg.$years.0$month.nc $fluxdir/nepfire.dif.2018_$years.0$month.nc
		else
			cdo mergetime -cat $outdir/combined.*$month.nc $outdir/merged.${month}.nc
			cdo yearmean $droughtdirectory/combined.2018$month.nc $droughtdirectory/combined.avg.2018$month.nc
			cdo timmean -cat $outdir/merged.$month.nc $outdir/merged.avg.$years.$month.nc
			cdo sub $droughtdirectory/combined.avg.2018$month.nc $outdir/merged.avg.$years.$month.nc $fluxdir/nepfire.dif.2018_$years.$month.nc
	fi
done


