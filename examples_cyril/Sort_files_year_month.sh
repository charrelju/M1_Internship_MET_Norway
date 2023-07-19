path_data="/lustre/storeB/project/copernicus/cosi/WP3/Data/ECMWF"
#
for year in {2012..2015}
do
   for mo in {1..12}
   do
      if [ $mo -lt 10 ]
      then
        month=0$mo
      else
        month=$mo
      fi
   echo $year, $month
   mkdir -p $path_data"/"$year"/"$month
   mv $path_data"/"$year"/ECMWF_operational_forecasts_T2m_10mwind_"$year$month* $path_data"/"$year"/"$month
   done
done
#
for year in {2012..2015}
do
	cd $path_data"/"$year"/"
	rmdir *
done
