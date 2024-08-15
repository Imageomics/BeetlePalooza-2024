gr3=/home/shared/group3
for img in A* ; do 
   domain=$(grep $img $gr3/beetleMeasurementsDomain.csv | tr -d '\r' | cut -d, -f 38 | uniq)
   if [ "$domain" == "NA" ] ; then
      echo "No domain found for image $img"
   else
      classes=$(grep $domain $gr3/taxonLists/foundtaxasubset.csv | tr -d '\r' | cut -d, -f 3 | sed -e 's/"//g' | python $gr3/names2list.py)
      #echo $classes
      echo "Running bioclip on the images in $img"
      bioclip predict --device cuda --cls "$classes" $img/*.png > $gr3/preds/${img}-preds.csv
   fi 
done
