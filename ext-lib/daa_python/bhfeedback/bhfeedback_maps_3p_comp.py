#>>> execfile('bhfeedback_maps_3p_comp.py')


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



basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
ICdir = 'h113_HR_sn214/'
simdir = basedir + ICdir
sims = [ 'bal_s7e1v30f05_n256',
         'bal_s8e1v30f05_n256', 
         'bal_s9e1v30f05_n256'
       ]

snap_ini = 1
snap_end = 291



#outdir = simdir + '/map3p'
outdir = './bhfeedback/' + ICdir + 'com3p_bal'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir


SkipIfDone = 0

acirc = np.linspace(0,2*np.pi,100)
xcirc = np.cos(acirc)
ycirc = np.sin(acirc)
R0 = 0.01

vmin = 6.3
vmax = 9.5


for snapnum in range(snap_ini,snap_end,1):
#for snapnum in [snap_ini,snap_end]:

 outfile = outdir + '/com3p_%03d.png' % snapnum
 if SkipIfDone and os.path.exists( outfile ):
   print "\nfile exists:  %s\n" % outfile
   continue

 print '\n --> ', outfile
 fig = plt.figure( figsize=(15.,5.) )
 fig.subplots_adjust(left=0., right=1., bottom=0., top=1., hspace=0., wspace=0.)
 ax1 = fig.add_subplot(131)
 ax2 = fig.add_subplot(132)
 ax3 = fig.add_subplot(133)

 for j, ax in enumerate([ax1, ax2, ax3]):

   gridfile = simdir + sims[j] + '/grids/grid_%03d.hdf5' % snapnum
   if not os.path.exists( gridfile ):
     print "\nfile not found:  %s\n" % gridfile
     continue
   grd = daa.load_dict_from_hdf5( gridfile )

   bh_file = simdir + sims[j] + '/bh_rbins.hdf5'
   bh = daa.load_dict_from_hdf5( bh_file )
   ind_snap = bh['snapnum'] == snapnum
   time = bh['time'][ind_snap] - bh['time'][0]
   print snapnum, bh['snapnum'][ind_snap], '  time =', time

   dx = grd['mgas1']['dx']
   xrg = [ -dx, dx ]
   extent = [xrg[0],xrg[1],xrg[0],xrg[1]]
   ax.set_xlim(xrg); ax.set_ylim(xrg); 
   ax.set_aspect('equal'); ax.axis('off')
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_visible(False)
   
   themap = np.log10( grd['mgas1']['g'][:,:] )
   im_temp = ax.imshow( themap, vmin=vmin, vmax=vmax, cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent )
   if j == 0:
     ax.text( 0.02, 0.98, r'$t = '+ '%05.2f$ Myr' % time, {'color':'black', 'fontsize':16}, bbox={'color':'white', 'pad':10.0},
           transform=ax.transAxes, horizontalalignment='left', verticalalignment='top' )
     
 fig.savefig( outfile )
 plt.close()





print '\nDone!\n'
