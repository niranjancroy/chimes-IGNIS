# run from python command line
#>>> import sys
#>>> sys.argv = ['bhfeedback_maps_6p.py', 'h113_HR_sn153/kick_s8e1v4p1', 0, 10 ]
#>>> execfile('bhfeedback_maps_6p.py')

#python bhfeedback_maps_6p.py  h113_HR_sn214/spawn_s9e1v30f05_w1m100  0  500  > log/map6p.0.log 2>&1 &


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
if narg == 4:
  sname = sys.argv[1]
  snap_ini = int(sys.argv[2])
  snap_end = int(sys.argv[3])
else:
  print 'syntax: python bhfeedback_maps_3p.py sname snap_ini snap_end'
  sys.exit()


#sname = 'h113_HR_sn214/spawn_s9e1v30f05_w1m1000_TshVsh_MERGETEST'
#snap_ini = 1
#snap_end = 2



SkipIfDone = 0


simdir = '/simons/scratch/dangles/FIRE/bhfeedback/' + sname

outdir = './bhfeedback/' + sname + '/map6p'
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

 outfile = outdir + '/map6p_%04d.png' % snapnum
 if SkipIfDone and os.path.exists( outfile ):
   print "\nfile exists:  %s\n" % outfile
   continue


 gridfile = simdir + '/grids/mgrid_%04d.hdf5' % snapnum
 if not os.path.exists( gridfile ):
   print "\nfile not found:  %s\n" % gridfile
   continue
 grd = daa.load_dict_from_hdf5( gridfile )

 dx1 = grd['faceon']['mgas0']['dx']
 xrg1 = [ -dx1, dx1 ]
 yrg1 = xrg1
 dx2 = grd['faceon']['mgas1']['dx']
 xrg2 = [ -dx2, dx2 ]
 yrg2 = xrg2
 dx3 = grd['faceon']['mgas2']['dx']
 xrg3 = [ -dx3, dx3 ]
 yrg3 = xrg3


 print '\n', gridfile, ' --> ', outfile
 fig = plt.figure( figsize=(15,10) )
 fig.subplots_adjust(left=0, right=1, bottom=0, top=1, hspace=0., wspace=0.)
 ax1 = fig.add_subplot(231)
 ax2 = fig.add_subplot(232)
 ax3 = fig.add_subplot(233)
 ax4 = fig.add_subplot(234)
 ax5 = fig.add_subplot(235)
 ax6 = fig.add_subplot(236)

 vmin = [ 6, 6, 6, 6, 6, 6 ]
 vmax = [ 10, 10, 10, 10, 10, 10 ]
 for j, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6]):
   if j < 3:
     ostr = 'faceon'
     jp = str(j)
   else:
     ostr = 'edgeon'
     jp = str(j-3)
   dx = grd[ostr]['mgas'+jp]['dx']
   xrg = [ -dx, dx ]
   extent = [xrg[0],xrg[1],xrg[0],xrg[1]]
   ax.set_xlim(xrg); ax.set_ylim(xrg); 
   ax.set_aspect('equal'); ax.axis('off')
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_visible(False)
   im_temp = ax.imshow( np.log10( grd[ostr]['mgas'+jp]['g'][:,:] ), vmin=vmin[j], vmax=vmax[j], cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent )
   if j in [0,3]:
     ax.text( 0.02, 0.98, r'$t = '+ '%03d$' % snapnum, {'color':'black', 'fontsize':16}, bbox={'color':'white', 'pad':0.4},
            transform=ax.transAxes, horizontalalignment='left', verticalalignment='top' )
     x_sq = np.array( [ xrg2[0], xrg2[1], xrg2[1], xrg2[0], xrg2[0] ] )
     y_sq = np.array( [ yrg2[0], yrg2[0], yrg2[1], yrg2[1], yrg2[0] ] )
     ax.plot( x_sq, y_sq, 'k-', lw=1 )
     x_line = np.array( [ xrg2[1], xrg1[1] ] )
     y_line = np.array( [ yrg2[1], yrg1[1] ] )
     ax.plot( x_line, y_line, 'k-', lw=0.5 )
     x_line = np.array( [ xrg2[1], xrg1[1] ] )
     y_line = np.array( [ yrg2[0], yrg1[0] ] )
     ax.plot( x_line, y_line, 'k-', lw=0.5 )
#     bar_size = 200
#     xleft = -300.
#     ybar = -350
#     ax.plot( [xleft,xleft+bar_size], [ybar,ybar], '-', color='white', linewidth=2 )
#     ax.text( (xleft+0.5*bar_size+400.)/800., (ybar+400.)/800., r'200 kpc', {'color':'white', 'fontsize':14},
#               transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom' )
   elif j in [1,4]:
     x_sq = np.array( [ xrg3[0], xrg3[1], xrg3[1], xrg3[0], xrg3[0] ] )
     y_sq = np.array( [ yrg3[0], yrg3[0], yrg3[1], yrg3[1], yrg3[0] ] )
     ax.plot( x_sq, y_sq, 'k-', lw=1 )
     x_line = np.array( [ xrg3[1], xrg2[1] ] )
     y_line = np.array( [ yrg3[1], yrg2[1] ] )
     ax.plot( x_line, y_line, 'k-', lw=0.5 )
     x_line = np.array( [ xrg3[1], xrg2[1] ] )
     y_line = np.array( [ yrg3[0], yrg2[0] ] )
     ax.plot( x_line, y_line, 'k-', lw=0.5 )
#     bar_size = 10
#     xleft = -15.
#     ybar = -17.5
#     ax.plot( [xleft,xleft+bar_size], [ybar,ybar], '-', color='black', linewidth=2 )
#     ax.text( (xleft+0.5*bar_size+20.)/40., (ybar+20.)/40., r'10 kpc', {'color':'black', 'fontsize':14},
#               transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom' )
   elif j in [2,5]:
     ax.plot( R0*xcirc, R0*ycirc, '-', color='black', linewidth=1.)
#     bar_size = 0.5
#     xleft = -0.75
#     ybar = -0.875
#     ax.plot( [xleft,xleft+bar_size], [ybar,ybar], '-', color='black', linewidth=2 )
#     ax.text( (xleft+0.5*bar_size+1.)/2., (ybar+1)/2, r'500 pc', {'color':'black', 'fontsize':14},
#            transform=ax.transAxes, horizontalalignment='center', verticalalignment='bottom' )
     
 fig.savefig( outfile )
 plt.close()





print '\nDone!\n'
