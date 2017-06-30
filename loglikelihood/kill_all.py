#!/usr/bin/python

import os, glob
#
# kill all the stoomboot jobs
#
TMPFILE = '/tmp/qlist.tmp'
cmd = 'qstat | grep jorana > ' + TMPFILE
os.system(cmd)

f = open(TMPFILE, 'r')
for line in f:
   job    = line.split()[0]
   # if int(job.split('.')[0]) < 18401667: continue 
   status = line.split()[4]
   cmd = 'qdel '+job
   print(status, job)
   if status == 'R' or status == 'Q':
     print('kill_all:: Killing job - '+job)
     os.system(cmd)

cmd = 'rm ' + TMPFILE
os.system(cmd)
