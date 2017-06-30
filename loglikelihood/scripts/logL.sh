#!/bin/sh 
export PATH=/project/datagrid/anaconda/bin:$PATH 
source activate joran 
python /data/xenon/mod_admin/Modulation/loglikelihood/logL.py 
