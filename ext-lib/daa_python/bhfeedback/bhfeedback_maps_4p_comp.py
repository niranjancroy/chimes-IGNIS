#>>> execfile('bhfeedback_maps_4p_comp.py')


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
sims = [ 'nof_s8e1_n256', 
         'bal_s9e1v30f05_n256',
         'bal_s8e1v30f05_n256',
         'bal_s7e1v30f05_n256',
         'kick_s8e1v15000p2_n256',
         'kick_s8e1v1500p20_n256',
         'kick_s7e1v15000p2_n256',
         'kick_s7e1v1500p20_n256' ]

snap_ini = 1
snap_end = 101



#outdir = simdir + '/map3p'
outdir = './bhfeedback/' + ICdir + sims[0] + '/com4p'
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

vmin = 6
vmax = 9
dx = 100.


for snapnum in range(snap_ini,snap_end,1):
#for snapnum in range(snap_ini,snap_end,20):

 outfile = outdir + '/com4p_%03d.png' % snapnum
 if SkipIfDone and os.path.exists( outfile ):
   print "\nfile exists:  %s\n" % outfile
   continue

 print '\n --> ', outfile
 fig = plt.figure( figsize=(10.,10.) )
 fig.subplots_adjust(left=0., right=1., bottom=0., top=1., hspace=0., wspace=0.)
 ax1 = fig.add_subplot(221)
 ax2 = fig.add_subplot(222)
 ax3 = fig.add_subplot(223)
 ax4 = fig.add_subplot(224)


 for j, ax in enumerate([ax1, ax2, ax3, ax4]):

   gridfile = simdir + sims[j] + '/grids/grid_%03d.hdf5' % snapnum
   if not os.path.exists( gridfile ):
     print "\nfile not found:  %s\n" % gridfile
     continue
   grd = daa.load_dict_from_hdf5( gridfile )

   dx0 = grd['mgas0']['dx']
   npix = grd['mgas0']['g'][:,0].size   
   dpix = 2*dx0/npix
   i = int( (dx0 - dx) / dpix )
   xrg = [ -dx, dx ]
   extent = [xrg[0],xrg[1],xrg[0],xrg[1]]
   ax.set_xlim(xrg); ax.set_ylim(xrg); 
   ax.set_aspect('equal'); ax.axis('off')
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_visible(False)
   
   if i > 0:
     themap = np.log10( grd['mgas0']['g'][i:-i,i:-i] )
   else:
     themap = np.log10( grd['mgas0']['g'][:,:] )
   im_temp = ax.imshow( themap, vmin=vmin, vmax=vmax, cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent )
   if j == 0:
     ax.text( 0.02, 0.98, r'$t = '+ '%03d$ Myr' % snapnum, {'color':'black', 'fontsize':16}, bbox={'color':'white', 'pad':0.4},
           transform=ax.transAxes, horizontalalignment='left', verticalalignment='top' )
   ax.text( 0.98, 0.98, sims[j], {'color':'black', 'fontsize':16}, bbox={'color':'white', 'pad':0.4},
           transform=ax.transAxes, horizontalalignment='right', verticalalignment='top' )
     
 fig.savefig( outfile )
 plt.close()





print '\nDone!\n'
