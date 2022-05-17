# run as >>> execfile('bhfeedback_makemovie.py')

import numpy as np
import os as os
from shutil import copyfile
import daa_lib as daa



#######################################
 ### DEFINE SIMULATION DIRECTORIES ###
#######################################


IC_list = [ #'h113_HR_sn155',
            #'h113_HR_sn214',
            #'h113_HR_sn277',
            'h113_HR_sn152'
          ]

sim_list = [ #'spawn_s9e1v30f05_w100dt1000m100_TshVsh',
             #'spawn_s8e1v30f05_w100dt1000m100_TshVsh',
             ##'spawn_s7e1v30f05_w100dt1000m100_TshVsh',
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
#             'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh'#,
#            'spawn_res251_s9e1v30f009_w10dt100m1000_TshVsh'
#             'TEST_mainbranch_e1f05',
#             'TEST_mainbranch_e1f009' 
             'spawn_res251_s9e01v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f033_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f02_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f009_w10dt100m1000_TshVsh'
           ]


bhfeedback_dir = '/home/dangles/repos/daa_python/bhfeedback/'
#bhfeedback_dir = '/simons/scratch/dangles/FIRE/bhfeedback/'


nsim = len(sim_list)


for ICdir in IC_list:
 for sname in sim_list:

  run_name = ICdir + '/' + sname
  log_name = ICdir + '_' + sname

  simdir = bhfeedback_dir + run_name + '/'

  
  mapdir = simdir + 'map6p/'
  movie_name = log_name + '_map6p.mov'
  filelist, filestr = daa.file_list( mapdir + 'map6p', ext='.png' )
  ffmpeg_str = 'ffmpeg -f image2 -framerate 24 -start_number ' + filestr[0] + ' -i ' + mapdir + 'map6p_%04d.png -r 24 -c:v h264 -crf 1 -pix_fmt yuv420p -y ' + simdir + movie_name
  #ffmpeg_str = 'ffmpeg -f image2 -framerate 24 -start_number 0000 -i ' + mapdir + 'map6p_%04d.png -r 24 -c:v h264 -crf 1 -pix_fmt yuv420p -y ' + simdir + movie_name
  #ffmpeg_str = 'ffmpeg -f image2 -framerate 12 -start_number 0000 -i ' + mapdir + 'map6p_%04d.png -r 24 -c:v h264 -crf 1 -pix_fmt yuv420p -y ' + simdir + movie_name

  """
  mapdir = simdir + 'tzgrid/'
  movie_name = log_name + '_windzoom.mov'
  filelist, filestr = daa.file_list( mapdir + 'windzoomplot', ext='.png' )
  ffmpeg_str = 'ffmpeg -f image2 -framerate 24  -start_number ' + filestr[0] + ' -i ' + mapdir + 'windzoomplot_%04d.png -r 24 -c:v h264 -crf 1 -pix_fmt yuv420p -vf scale=3000:-2  -y ' + simdir + movie_name 
  """
   
  """
  mapdir = simdir + 'tzgrid/'
  movie_name = log_name + '_tzgrid.mov'
  #ffmpeg_str = 'ffmpeg -f image2 -framerate 24 -start_number 0000 -i ' + mapdir + 'tzgrid_%04d.png -r 24 -c:v h264 -crf 1 -pix_fmt yuv420p -y ' + simdir + movie_name
  ffmpeg_str = 'ffmpeg -f image2 -framerate 48 -start_number 0000 -i ' + mapdir + 'tzgrid_%04d.png -r 48 -c:v h264 -crf 1 -pix_fmt yuv420p -y ' + simdir + movie_name
  """ 

  print '\n'+ffmpeg_str+'\n'
  os.system(ffmpeg_str)
  try:
    copyfile( simdir + movie_name, '/home/dangles/repos/daa_python/bhfeedback/movies/' + movie_name )
  except IOError:
    print '  --> file not found: ', simdir + movie_name, '/home/dangles/repos/daa_python/bhfeedback/movies/' + movie_name


print '\n\nDone!!'





