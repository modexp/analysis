#!/usr/bin/python



# This is the script that submits the jobs that make the log-likelihood analysis.
# Joran Angevaare 29-06-2017


import os, glob
import time

print("SUBMITSCRIPT:: Run log-likelihood analysis do the logL.py")

          #  Ch0, Ch1,  Ch2 , Ch3 , Ch4   , Ch5   , Ch6  , Ch7
halflife  = [1e9,  1e9, 59.1, 59.1, 5.2711, 5.2711, 30.08, 30.08]
dhalflife = [1e9,  1e9, 0.3,  0.3,  0.0004, 0.0004,  0.09,  0.09]
logL_dir  = "/data/xenon/mod_admin/Modulation/loglikelihood/"
logfiles  = "/dcache/xenon/jorana/log_files/"
N_DATASETS= 10
username  = 'jorana'

def check_qstat():
    '''Reruns how may files are running, in que and completed '''
    TMPFILE = 'qlist.tmp'
    cmd = 'qstat -u ' + username + ' | grep mc_  > ' + TMPFILE
    # print(cmd)
    os.system(cmd)
    f = open(TMPFILE, 'r')
    n, m, o = 0, 0, 0
    for line in f:
        job    = line.split()[0]
        status = line.split()[9]
        # cmd = 'qdel '+job
        if status == 'R': n = n + 1
        if status == 'Q': m = m + 1
        if status == 'C': o = o + 1
    cmd = 'rm ' + TMPFILE
    os.system(cmd)
    return [n, m, o]

def submit_fakes(ch,pk,tau,amp,phi,num, fix_a):
    '''Submit a job to the stoomboot que'''
    print("SUBMITSCRIPT:: submit fitting make_and_fit.py")
    if a > 0 : return 0
    queue = 'long'
    run   = 'mc_%i_%i_%.7f'%(ch, pk, amp)
    scriptfile = logL_dir + 'scripts/'+ run + '.sh'
    job = 'python ' + logL_dir + 'make_and_fit.py -ch %i -pk %i -tau %f -amp %.7f -phi %f -num %i -fixa %f'%(ch, pk, tau, amp, phi, num, fix_a)

    fout = open(scriptfile,'w')
    fout.write('#!/bin/sh \n')
    fout.write('export PATH=/project/datagrid/anaconda/bin:$PATH \n')
    fout.write('source activate joran \n')
    if os.path.exists(logfiles + run) or os.path.exists(logfiles + 'done_' + run):  
        print('MAIN::remove ' + logfiles + run)
        fout.write('cd ' +  logfiles  +' \n')
        if os.path.exists(logfiles + run): fout.write('rm -f ' +  run +' \n')
        if os.path.exists(logfiles + 'done_' + run): fout.write('rm -f ' + 'done_' + run +' \n')
        fout.write('cd ' +  logL_dir +' \n')
    fout.write('exec 1> ' + logfiles + run + ' \n' )
    fout.write('exec 2>&1 \n')
    fout.write(job + '\n')
    fout.write('cd '+ logfiles +' \n')
    fout.write('mv ' + run + ' done_' + run + ' \n') 
    fout.close()

    # submit the job to stoomboot
    cmd = 'qsub -W group_list=xenon -l walltime=96:00:00 -e logfiles -o logfiles -j eo -q '+ queue + ' ' + scriptfile 
    print("SUBMITSCRIPT:: Submit %s"%scriptfile)
    os.system(cmd)

def string(integer):
    integer = "000000000000000000000000000%i"%integer
    return integer[-3:]


print('SUBMITSCRIPT::  submit jobs')

# These are the modulatiopn amplitudes that I want to simulate since numpy is not standard to stoomboot I hardcorded this:
# import numpy as np
# mod_a = np.concatenate((np.linspace(0, 1e-5, 21), np.linspace(1e-5, 1e-4, 31), np.linspace(1e-4, 1e-3, 16)))
mod_a = [0.00000000e+00,   5.00000000e-07,   1.00000000e-06,         1.50000000e-06,   2.00000000e-06,   2.50000000e-06,         3.00000000e-06,   3.50000000e-06,   4.00000000e-06,
         4.50000000e-06,   5.00000000e-06,   5.50000000e-06,         6.00000000e-06,   6.50000000e-06,   7.00000000e-06,         7.50000000e-06,   8.00000000e-06,   8.50000000e-06,
         9.00000000e-06,   9.50000000e-06,   1.00000000e-05,         1.00000000e-05,   1.30000000e-05,   1.60000000e-05,         1.90000000e-05,   2.20000000e-05,   2.50000000e-05,
         2.80000000e-05,   3.10000000e-05,   3.40000000e-05,         3.70000000e-05,   4.00000000e-05,   4.30000000e-05,         4.60000000e-05,   4.90000000e-05,   5.20000000e-05,
         5.50000000e-05,   5.80000000e-05,   6.10000000e-05,         6.40000000e-05,   6.70000000e-05,   7.00000000e-05,         7.30000000e-05,   7.60000000e-05,   7.90000000e-05,
         8.20000000e-05,   8.50000000e-05,   8.80000000e-05,         9.10000000e-05,   9.40000000e-05,   9.70000000e-05,         1.00000000e-04,   1.00000000e-04,   1.60000000e-04,
         2.20000000e-04,   2.80000000e-04,   3.40000000e-04,         4.00000000e-04,   4.60000000e-04,   5.20000000e-04,         5.80000000e-04,   6.40000000e-04,   7.00000000e-04,
         7.60000000e-04,   8.20000000e-04,   8.80000000e-04,         9.40000000e-04,   1.00000000e-03]

# Submit a job for every modulation amplitude
for j, a in enumerate(mod_a):
    if j%10 == 0:
        # There is a maximum of ~ 200 jobs that can be in the que. Therefore wait if there are many jobs in the que, otherwise submit anouther 70 jobs into the que
        for i in range(10 * 60):
            res = check_qstat()
            if res[1] < 200 - 70: continue 
            else: print("SUBMITSCRIPT:: There are still %s files in the que (R = %s, Q = %s, C = %s) , waited %i minutes"%(string(sum(res)), string(res[0]), string(res[1]), string(res[2]), i)); time.sleep(60)
    for ch in [2,4,6]:
        t12 = halflife[ch]
        for pk in [0,1,2]:
            if ch >4 and pk > 0: continue
            submit_fakes(ch, pk , t12,  a, 1, N_DATASETS, -1)
    

queue = 'long'
scriptfile = 'scripts/logL' + '.sh'
fout = open(scriptfile,'w')
fout.write('#!/bin/sh \n')
fout.write('export PATH=/project/datagrid/anaconda/bin:$PATH \n')
fout.write('source activate joran \n')
fout.write('python ' + logL_dir + 'logL.py \n')
fout.close()

# submit the job to stoomboot
cmd = 'qsub -W group_list=xenon -l walltime=10:00:00 -e logfiles -o logfiles -j eo -q '+ queue + ' ' + scriptfile 
os.system(cmd)

print('SUBMITSCRIPT:: Exit from the stoomboot submit script. bye-bye.')
print('SUBMITSCRIPT:: check the stats of your job with qstat')

