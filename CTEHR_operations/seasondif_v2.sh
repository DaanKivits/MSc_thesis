#!/bin/bash
fluxfolder='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'

echo 'give output directory name: '
read outdir

if [ ! -d $outdir ]; then
  mkdir $outdir
fi

droughtfolder='/projects/0/ctdas/dkivits/DATA/fluxes/2018'
fluxdir='/projects/0/ctdas/dkivits/DATA/fluxes/'

if [ ! -d $droughtfolder ]; then
  mkdir $droughtfolder
fi

mv "$outdir"/* "$outdir/backup"
mv "$droughtfolder"/* "$droughtfolder/backup"

years={2017..2021}
#years=2020

for month in {7..9}
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

for year in {2017..2021}
#for year in $years
do
        if (($year!=2018)) && (($year!=2022))
	then
		for month in {7..9}
                do
			if (($month<10))
			then
				cp $fluxfolder/nep.${year}0$month.nc $outdir/nep.${year}0$month.nc
                		cp $fluxfolder/fire.${year}0$month.nc $outdir/fire.${year}0$month.nc
			else
				cp $fluxfolder/nep.${year}$month.nc $outdir/nep.${year}$month.nc
                                cp $fluxfolder/fire.${year}$month.nc $outdir/fire.${year}$month.nc
			fi
		done
	fi
	
	if (($year==2022))
	then
		for month in {7..8}
		do
			if (($month<10))
                        then
                                cp $fluxfolder/nep.${year}0$month.nc $outdir/nep.${year}0$month.nc
                                cp $fluxfolder/fire.${year}0$month.nc $outdir/fire.${year}0$month.nc
                        else
                                cp $fluxfolder/nep.${year}$month.nc $outdir/nep.${year}$month.nc
                                cp $fluxfolder/fire.${year}$month.nc $outdir/fire.${year}$month.nc
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

cdo mergetime $droughtfolder/nep.2018*.nc $droughtfolder/nep.merged.2018.nc
cdo mergetime $droughtfolder/fire.2018*.nc $droughtfolder/fire.merged.2018.nc
cdo add $droughtfolder/nep.merged.2018.nc $droughtfolder/fire.merged.2018.nc $droughtfolder/combined.merged.2018.nc
cdo mergetime -cat $outdir/combined.merged.*.nc $outdir/combined.merged.nc

cdo yearmean $outdir/combined.merged.nc $outdir/average.JAS.multiyear.nc
cdo timmean $outdir/average.JAS.multiyear.nc $outdir/average.JAS.$years.nc
cdo yearmean $droughtfolder/combined.merged.2018.nc $droughtfolder/average.JAS.2018.nc

cdo sub $droughtfolder/average.JAS.2018.nc $outdir/average.JAS.$years.nc $fluxdir/nepfire.dif.JAS.2018_$years.nc

