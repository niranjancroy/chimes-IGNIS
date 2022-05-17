# >>> execfile('bhfeedback_multigrids_submit.py')

import numpy as np
import gadget as g
import os as os
import sys as sys
from glob import glob
import daa_lib as daa




IC_list = [ #'h113_HR_sn155',
            #'h113_HR_sn214',
            'h113_HR_sn277'
          ]

sim_list = [ #'spawn_s9e1v30f05_w100dt1000m100_TshVsh'#,
             #'spawn_s8e1v30f05_w100dt1000m100_TshVsh',
             #'spawn_s7e1v30f05_w100dt1000m100_TshVsh',
             #'nof_s8e1_n128',
             'spawn_s9e1v30f05_w100dt1000m100_TshVsh_MERGE'
           ]


#isim = 15
#sim_list = [ sim_list[isim] ]


nsnap_per_task = 2
#nsnap_per_task = 1

nsim = len(sim_list)


for ICdir in IC_list:
 for sname in sim_list:

  run_name = ICdir + '/' + sname
  log_name = ICdir + '_' + sname

  bash_job_str = \
  '#!/bin/bash' + '\n' + \
  '#SBATCH --job-name="' +  log_name +'"' + '\n' + \
  '#SBATCH --output="gordon.log"' + '\n' + \
  '#SBATCH --partition=general' + '\n' + \
  '#SBATCH --nodes=1' + '\n' + \
  '#SBATCH --ntasks-per-node=16' + '\n' + \
  '#SBATCH --mem=60000' + '\n' + \
  '#SBATCH --exclusive' + '\n' + \
  '#SBATCH --export=ALL' + '\n' + \
  '#SBATCH -t 168:00:00' + '\n' + \
  '###SBATCH --mail-type=ALL' + '\n' + \
  '###SBATCH --mail-user=danglesalcazar@flatironinstitute.org' + '\n\n' + \
  'cd /home/dangles/repos/daa_python' + '\n\n' + \
  'module load python' + '\n' + \
  'module load hdf5' + '\n' + \
  'module load scipy' + '\n\n' 

#  snapdir = glob( '/simons/scratch/dangles/FIRE/bhfeedback/' + run_name + '/snapdir*')
#  nsnap = len(snapdir)
#  snaplist = np.zeros(nsnap,dtype='int')
#  for i in range(nsnap):
#     snaplist[i] = int( snapdir[i][ snapdir[i].find('_', -5, -1)+1: ] )
#  snaplist.sort()

  bhdata = daa.load_dict_from_hdf5( '/simons/scratch/dangles/FIRE/bhfeedback/' + run_name + '/centralBH.hdf5')
  snaplist = bhdata['snapnum']
  nsnap = snaplist.size
  print snaplist

  for n, i in enumerate(range(0,nsnap,nsnap_per_task)):
     i_ini = snaplist[i]
     i_end = i+nsnap_per_task
     if i_end > nsnap-1:
       i_end = nsnap-1
       i_end = snaplist[i_end] + 1
     else:
       i_end = snaplist[i_end]
     py_str = 'python bhfeedback_multigrids.py  ' + run_name + '  ' + str(i_ini) + '  ' + str(i_end) + '  > log/' + log_name + '.' + str(n) + '.log\n\n'
     print py_str
     sh_name =  'tmp.sh'
     f = open(sh_name,'w')
     f.write(bash_job_str + py_str)
     f.close()

     os.system('sbatch ' + sh_name)
     os.remove(sh_name)


print '\n\nDone!\n'


