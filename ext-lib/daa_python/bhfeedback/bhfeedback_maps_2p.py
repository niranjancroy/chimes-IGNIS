#>>> execfile('bhfeedback_maps_2p.py')

# python bhfeedback_maps_2p.py  h113_HR_sn153/nofback  0  200  > log/map2p_nofback.log 2>&1 &


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
import os as os
import sys as sys
import numpy as np
from scipy.stats import binned_statistic_2d
import h5py as h5py
import pandas as pd
from glob import glob
import gadget as g
import utilities as util
import visualization as vis
import daa_lib as daa
from daa_constants import *

import imp
imp.reload(daa)
imp.reload(util)


narg = len(sys.argv)
if narg == 5:
  basedir = sys.argv[1]
  sname = sys.argv[2]
  snap_ini = int(sys.argv[3])
  snap_end = int(sys.argv[4])
else:
  print 'syntax: python bhfeedback_maps_2p.py basedir sname snap_ini snap_end'
  sys.exit()

"""
snap_ini = 0
snap_end = 2
sname = 'h113_HR_sn153/nofback'
basedir = '/simons/scratch/dangles/FIRE/bhfeedback/old/'
#movie_name = 'h113_HR_sn153_2p.mov'
movie_name = 'h113_HR_sn153_2pRot.mov'
"""

SkipIfDone = 0


simdir = basedir + sname

#outdir = simdir + '/map3p'
outdir = './bhfeedback/' + sname + '/map2p'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir


acirc = np.linspace(0,2*np.pi,100)
xcirc = np.cos(acirc)
ycirc = np.sin(acirc)
R0 = 0.01


for snapnum in range(snap_ini,snap_end,1):

 outfile = outdir + '/map2p_%03d.png' % snapnum
 if SkipIfDone and os.path.exists( outfile ):
   print "\nfile exists:  %s\n" % outfile
   continue


 #gridfile = simdir + '/grids/grid_%03d.hdf5' % snapnum
 gridfile = simdir + '/grids/mgrid_%04d.hdf5' % snapnum
 if not os.path.exists( gridfile ):
   print "\nfile not found:  %s\n" % gridfile
   continue
 grd = daa.load_dict_from_hdf5( gridfile )
 
 grd = grd['edgeon']

 dx1 = grd['T0']['dx']
 xrg1 = [ -dx1, dx1 ]
 yrg1 = xrg1
 dx2 = grd['T1']['dx']
 xrg2 = [ -dx2, dx2 ]
 yrg2 = xrg2
 #dx3 = grd['T2']['dx']
 #xrg3 = [ -dx3, dx3 ]
 #yrg3 = xrg3


 print '\n', gridfile, ' --> ', outfile
 fig = plt.figure( figsize=(10.,5.), dpi=300 )
 fig.subplots_adjust(left=0., right=1., bottom=0., top=1., hspace=0., wspace=0.)
 ax1 = fig.add_subplot(121)
 ax2 = fig.add_subplot(122)
 #ax3 = fig.add_subplot(133)

 vmin = [ 6.25, 7 ]
 vmax = [ 9.5, 9.5 ]
 for j, ax in enumerate([ax1, ax2]):
   jp = str(j)
   dx = grd['T'+jp]['dx']
   xrg = [ -dx, dx ]
   extent = [xrg[0],xrg[1],xrg[0],xrg[1]]
   ax.set_xlim(xrg); ax.set_ylim(xrg); 
   ax.set_aspect('equal'); ax.axis('off')
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_visible(False)
   im_temp = ax.imshow( np.log10( grd['mgas'+jp]['g'][:,:] ), vmin=vmin[j], vmax=vmax[j], cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent )
   if j==0:
     #papertxt = r'$\rm{Angl\'es-Alc\'azar \, et \, al.\, 2017, \, MNRAS, \, 472, \, L109}$'
     papertxt = r'$\rm{Angl\'es}$-${\rmAlc\'azar \, et \, al.\, 2017, \, MNRAS, \, 472, \, L109}$'
     ax.text( 0.01, 0.005, papertxt, {'color':'orange', 'fontsize':8}, transform=ax.transAxes, horizontalalignment='left', verticalalignment='bottom' )
     x_sq = np.array( [ xrg2[0], xrg2[1], xrg2[1], xrg2[0], xrg2[0] ] )
     y_sq = np.array( [ yrg2[0], yrg2[0], yrg2[1], yrg2[1], yrg2[0] ] )
     ax.plot( x_sq, y_sq, 'k:', lw=1 )
     x_line = np.array( [ xrg2[1], xrg1[1] ] )
     y_line = np.array( [ yrg2[1], yrg1[1] ] )
     ax.plot( x_line, y_line, 'k:', lw=1 )
     x_line = np.array( [ xrg2[1], xrg1[1] ] )
     y_line = np.array( [ yrg2[0], yrg1[0] ] )
     ax.plot( x_line, y_line, 'k:', lw=1 )
     
 fig.savefig( outfile, dpi=300 )
 plt.close()


#ffmpeg_str = 'ffmpeg -f image2 -framerate 24 -start_number 000 -i ' + outdir + '/map2p_%03d.png -r 24 -c:v h264 -crf 1 -pix_fmt yuv420p -y ' + outdir + '/' + movie_name
#print '\n'+ffmpeg_str+'\n'
#os.system(ffmpeg_str)

print '\nDone!\n'


