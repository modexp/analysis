#!/bin/sh 
export PATH=/project/datagrid/anaconda/bin:$PATH 
source activate joran 
cd /dcache/xenon/jorana/log_files/ 
rm -f done_mc_2_0_0.0000000 
cd /data/xenon/mod_admin/Modulation/loglikelihood/ 
exec 1> /dcache/xenon/jorana/log_files/mc_2_0_0.0000000 
exec 2>&1 
python /data/xenon/mod_admin/Modulation/loglikelihood/make_and_fit.py -ch 2 -pk 0 -tau 59.100000 -amp 0.0000000 -phi 1.000000 -num 10 -fixa -1.000000
cd /dcache/xenon/jorana/log_files/ 
mv mc_2_0_0.0000000 done_mc_2_0_0.0000000 
