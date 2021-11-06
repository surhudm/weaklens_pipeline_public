#### This script shows how to link the files downloaded from Xiangchong
#### Please change the linked files according to where you downloaded them

mkdir DataStore
cd DataStore

mkdir Blinded_S19A_v2.0/ 
mkdir S19A_v2.0/ 

# Blinded_S19A_v2.0 : Directory structure
#
# Link all the _m.fits files
# http://gfarm.ipmu.jp/~xiangchong.li/S19ACatalogs/v2/catalog_obs_reGaus_blinded_m/
#
cd Blinded_S19A_v2.0/
for file in /gpfs02/work/xiangchong.li/S19ACatalogs/v2/catalog_obs_reGaus_blinded_m/surhud-3blind/*_m.fits
do
    ln -svf $file
done
cd ..

#
# Link all the tractwise shape catalog files
#
# http://gfarm.ipmu.jp/~xiangchong.li/S19ACatalogs/v2/catalog_tracts/
cd S19A_v2.0
for dir in /gpfs02/work/xiangchong.li/S19ACatalogs/v2/catalog_tracts/*_tracts
do
    ln -svf $dir
done

# tract_list.dat
rm tract_list.dat
echo field tract > tract_list.dat
for file in *_tracts/*_no_m.fits
do
    echo $file | sed 's/_tracts\// /g' | sed 's/_no_m.fits//g' >> tract_list.dat
done
#python ../../generate_mapping.py S19A_v2.0
cd ..

# Back up to the main directory
cd ..
