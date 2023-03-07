#!/bin/bash
fluxdirectory='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'

echo 'give output directory name: '
read outdir

if [ ! -d $outdir ]; then
  mkdir $outdir
fi

droughtdirectory='/projects/0/ctdas/dkivits/DATA/Fluxes_CDO/2018'
fluxdir='/projects/0/ctdas/dkivits/DATA/Fluxes_CDO'

if [ ! -d $droughtdirectory ]; then
  mkdir $droughtdirectory
fi

mv "$outdir"/* "$outdir/backup"
mv "$droughtdirectory"/* "$droughtdirectory/backup"

years={2017..2021}
#years=2020

for month in {7..9}
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

for year in {2017..2021}
#for year in $years
do
        if (($year!=2018)) && (($year!=2022))
	then
		for month in {7..9}
                do
			if (($month<10))
			then
				cp $fluxdirectory/nep.${year}0$month.nc $outdir/nep.${year}0$month.nc
				cp $fluxdirectory/fire.${year}0$month.nc $outdir/fire.${year}0$month.nc
			else
				cp $fluxdirectory/nep.${year}$month.nc $outdir/nep.${year}$month.nc
				cp $fluxdirectory/fire.${year}$month.nc $outdir/fire.${year}$month.nc
			fi
		done
	fi
	
	if (($year==2022))
	then
		for month in {7..8}
		do
			if (($month<10))
                        then
							cp $fluxdirectory/nep.${year}0$month.nc $outdir/nep.${year}0$month.nc
							cp $fluxdirectory/fire.${year}0$month.nc $outdir/fire.${year}0$month.nc
                        else
							cp $fluxdirectory/nep.${year}$month.nc $outdir/nep.${year}$month.nc
							cp $fluxdirectory/fire.${year}$month.nc $outdir/fire.${year}$month.nc
                        fi
		done
	fi
	
	if (($year!=2018))
	then
		cdo mergetime -cat $outdir/nep.${year}*.nc $outdir/nep.merged.${year}.nc
		cdo mergetime -cat $outdir/fire.${year}*.nc $outdir/fire.merged.${year}.nc
		cdo add $outdir/nep.merged.${year}.nc $outdir/fire.merged.${year}.nc $outdir/combined.merged.{$year}.nc
	fi
done

cdo mergetime $droughtdirectory/nep.2018*.nc $droughtdirectory/nep.merged.2018.nc
cdo mergetime $droughtdirectory/fire.2018*.nc $droughtdirectory/fire.merged.2018.nc
cdo add $droughtdirectory/nep.merged.2018.nc $droughtdirectory/fire.merged.2018.nc $droughtdirectory/combined.merged.2018.nc
cdo mergetime -cat $outdir/combined.merged.*.nc $outdir/combined.merged.nc

cdo yearmean $outdir/combined.merged.nc $outdir/average.JAS.multiyear.nc
cdo timmean $outdir/average.JAS.multiyear.nc $outdir/average.JAS.$years.nc
cdo yearmean $droughtdirectory/combined.merged.2018.nc $droughtdirectory/average.JAS.2018.nc

cdo sub $droughtdirectory/average.JAS.2018.nc $outdir/average.JAS.$years.nc $fluxdir/nepfire.dif.JAS.2018_$years.nc

