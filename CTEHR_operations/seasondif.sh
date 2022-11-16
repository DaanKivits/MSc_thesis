#!/bin/bash
fluxfolder='/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT'

echo 'give output directory name: '
read outdir

if [ ! -d $outdir ]; then
  mkdir $outdir
fi

if [ ! -d backup ]; then
  mkdir backup
fi

mv "$outdir"/* "backup"

droughtfolder='/projects/0/ctdas/dkivits/DATA/fluxes/2018'

#years={2017..2022}
years=2017

for month in {7..9}
do
                cp $fluxfolder/nep.20180$month.nc $droughtfolder/nep.20180$month.nc
                cp $fluxfolder/fire.20180$month.nc $droughtfolder/fire.20180$month.nc
		cdo add $droughtfolder/nep.20180$month.nc $droughtfolder/fire.20180$month.nc $droughtfolder/combined.20180$month.nc
		cdo timmean $droughtfolder/combined.20180$month.nc $droughtfolder/combined.avg.20180$month.nc
done

#for year in {2017..2022}
for year in $years
do
        if (($year!=2018)) && (($year!=2022))
	then
		for month in {7..9}
                do
			cp $fluxfolder/nep.${year}0$month.nc $outdir/nep.${year}0$month.nc
                	cp $fluxfolder/fire.${year}0$month.nc $outdir/fire.${year}0$month.nc
                	cdo add $outdir/nep.${year}0$month.nc $outdir/fire.${year}0$month.nc $outdir/combined.${year}0$month.nc
                	cdo timmean $outdir/combined.${year}0$month.nc $outdir/combined.avg.${year}0$month.nc
		done
	fi
	
	if (($year==2022))
	then
		for month in {7..8}
		do
			cp $fluxfolder/nep.${year}0$month.nc $outdir/nep.${year}0$month.nc
                        cp $fluxfolder/fire.${year}0$month.nc $outdir/fire.${year}0$month.nc
                        cdo add $outdir/nep.${year}0$month.nc $outdir/fire.${year}0$month.nc $outdir/combined.${year}0$month.nc
                        cdo timmean $outdir/combined.${year}0$month.nc $outdir/combined.avg.${year}0$month.nc
		done
	fi
	
	if (($year!=2018))
	then
		cdo timmean -cat $outdir/combined.avg.${year}*.nc $outdir/combined.avg.peryear.${year}.nc
	fi
done


cdo timmean -cat $outdir/combined.avg.peryear.*.nc $outdir/average.JAS.$years.nc 
cdo timmean -cat $droughtfolder/combined.avg.2018*.nc $droughtfolder/average.JAS.2018.nc
cdo sub $droughtfolder/average.JAS.2018.nc $outdir/average.JAS.$years.nc nepfire.dif.JAS.2018_$years.nc

