#!/bin/sh 
export PATH=/project/datagrid/anaconda/bin:$PATH 
source activate joran 
cd /dcache/xenon/jorana/log_files/ 
rm -f done_mc_4_2_0.0000000 
cd /data/xenon/mod_admin/Modulation/loglikelihood/ 
exec 1> /dcache/xenon/jorana/log_files/mc_4_2_0.0000000 
exec 2>&1 
python /data/xenon/mod_admin/Modulation/loglikelihood/make_and_fit.py -ch 4 -pk 2 -tau 5.271100 -amp 0.0000000 -phi 1.000000 -num 10 -fixa -1.000000
cd /dcache/xenon/jorana/log_files/ 
mv mc_4_2_0.0000000 done_mc_4_2_0.0000000 
