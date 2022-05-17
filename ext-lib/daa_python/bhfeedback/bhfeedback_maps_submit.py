# >>> execfile('bhfeedback_maps_submit.py')
# python bhfeedback_maps_submit.py

import numpy as np
import gadget as g
import os as os
import sys as sys



IC_list = [ #'h113_HR_sn155',
            #'h113_HR_sn214',
            #'h113_HR_sn277',
            'h113_HR_sn152'
          ]

sim_list = [ #'spawn_s9e1v30f05_w100dt1000m100_TshVsh',
             #'spawn_s8e1v30f05_w100dt1000m100_TshVsh',
             #'spawn_s7e1v30f05_w100dt1000m100_TshVsh',
             #'nof_s8e1_n128'
             #'spawn_s9e1v30f05_w100dt1000m100_TshVsh_MERGE',
             #'spawn_s9e1v30f05_w10dt100m1000_TshVsh',
             #'spawn_s9e1v30f05_w100dt1000m100',
             #'spawn_s9e1v30f033_w100dt1000m100_TshVsh',
             #'spawn_s9e1v30f009_w100dt1000m100_TshVsh',
             #'kick_s9e1v15000p2_n128',
             #'kick_s9e1v1500p20_n128',
             #'bal_s9e1v30f05_n128',
             #'kick_COL_s9e1v15000p2_n128',
             #'kick_COL_s9e1v1500p20_n128'
#             'nof_s8e1_n128',
#             'spawn_res251_s9e01v30f05_w10dt100m1000_TshVsh',
#             'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh',
#             'spawn_res251_s9e1v30f009_w10dt100m1000_TshVsh',
#             'TEST_mainbranch_e1f05',
#             'TEST_mainbranch_e1f009'
             'spawn_res251_s9e01v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f033_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f02_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f009_w10dt100m1000_TshVsh'
           ]


#isim = 15
#sim_list = [ sim_list[isim] ]

nsim = len(sim_list)

for ICdir in IC_list:
 for sname in sim_list:

  run_name = ICdir + '/' + sname
  log_name = ICdir + '_' + sname

  bash_job_str = \
  '#!/bin/bash' + '\n' + \
  '#SBATCH --job-name="' + log_name +'"' + '\n' + \
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

  py_str = 'python bhfeedback_maps_6p.py  ' + run_name + ' 0  2000  > log/map6p_' + log_name + '.log\n\n'
  print py_str
  sh_name = 'tmp.sh'
  f = open(sh_name,'w')
  f.write(bash_job_str + py_str)
  f.close()
  os.system('sbatch ' + sh_name)
  os.remove(sh_name)


print '\n\nDone!\n'


