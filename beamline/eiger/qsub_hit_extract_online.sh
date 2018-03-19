#$ -S /bin/bash
#$ -j y
#$ -V

export LANG=C

echo "Start: `date`"
hostname
echo "In $PWD"

#yamtbx.python /isilon/blconfig/bl32xu/local_bss/yam/hit_extract_online.py \
yamtbx.python /oys/xtal/yamtbx/bl32xu/eiger/hit_extract_online.py \
  --min-spots=3 --ctime-master=$ctime_master --tmpdir=$tmpdir \
  $master_h5 $dbfile 

echo "Finished: `date`"
echo "SECONDS= $SECONDS"
