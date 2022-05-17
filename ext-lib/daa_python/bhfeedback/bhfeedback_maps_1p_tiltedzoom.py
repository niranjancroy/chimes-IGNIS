#>>> execfile('bhfeedback_maps_1p_tiltedzoom.py')

# python bhfeedback_maps_1p_tiltedzoom.py  h113_HR_sn153/nofback  0  200  > log/map2p_nofback.log 2>&1 &


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
  print 'syntax: python bhfeedback_maps_1p_tiltedzoom.py basedir sname snap_ini snap_end'
  sys.exit()


"""
basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
sname = 'h113_HR_sn277/spawn_s9e1v30f05_w100dt1000m100_TshVsh'
snap_ini = 1000
snap_end = 1001
movie_name = 'h113_HR_sn277_s9_1p_tiltedzoom.mov'
"""


SkipIfDone = 0


simdir = basedir + sname

#outdir = simdir + '/map3p'
outdir = './bhfeedback/' + sname + '/tzgrid'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir


galdata = daa.load_dict_from_hdf5( '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_galdata_fbh_00000.hdf5' )
acirc = np.linspace(0,2*np.pi,100)
xcirc = np.cos(acirc)
ycirc = np.sin(acirc)
ind_snap = np.where( galdata['snapnum'] == int(sname[10:13]))[0]
Rvir = galdata['halo']['Rvir'][ind_snap[0]]



for snapnum in range(snap_ini,snap_end,1):
#for snapnum in [1,100,200,300,400,500,600,700,800,900,1000]:

 outfile = outdir + '/tzgrid_%04d.png' % snapnum
 if SkipIfDone and os.path.exists( outfile ):
   print "\nfile exists:  %s\n" % outfile
   continue

 gridfile = simdir + '/tzgrid/tzgrid_%04d.hdf5' % snapnum
 if not os.path.exists( gridfile ):
   print "\nfile not found:  %s\n" % gridfile
   continue
 grd = daa.load_dict_from_hdf5( gridfile )
 #grd = grd['edgeon']

 print '\n', gridfile, ' --> ', outfile
 fig = plt.figure( figsize=(5.,5.), dpi=150 )
 fig.subplots_adjust(left=0., right=1., bottom=0., top=1., hspace=0., wspace=0.)
 ax = fig.add_subplot(111)

 dx = grd['dx']
 xrg = [ -dx, dx ]
 extent = [xrg[0],xrg[1],xrg[0],xrg[1]]
 ax.set_xlim(xrg); ax.set_ylim(xrg); 
 ax.set_aspect('equal'); ax.axis('off')
 for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_visible(False)


 #cmap_gas = cm.Greys
 #cmap_wind = cm.Spectral

 cmap_gas = cm.Greys
 cmap_wind = cm.jet 
 alpha_wind = 0.25
 

 # --- background gas density
 vgas = [ 5.5, 9.5 ]
 gas_map = grd['gas']['m'][:,:]
 im = ax.imshow( np.log10( gas_map ), vmin=vgas[0], vmax=vgas[1], cmap=cmap_gas, interpolation='bicubic', origin='low', extent=extent )

 # --- AGN wind 
 vwind = [ 3.5, 6 ]
 wind_grid = grd['wind']['m'][:,:]
 wind_map = np.ma.masked_array( wind_grid, wind_grid <= wind_grid.min() )
 im = ax.imshow( np.log10( wind_map ), vmin=vwind[0], vmax=vwind[1], cmap=cmap_wind, interpolation='bicubic', origin='low', extent=extent, alpha=alpha_wind )
 
 """
 ek_grid = grd['wind']['ek'][:,:]
 ek_map = np.ma.masked_array( ek_grid, ek_grid <= ek_grid.min() )
 im = ax.imshow( np.log10( ek_map ), vmin=11, vmax=14, cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent, alpha=0.4 )

 Vr_grid = grd['wind']['Vr'][:,:]
 Vr_map = np.ma.masked_array( np.abs(Vr_grid), wind_grid <= wind_grid.min() )
 im = ax.imshow( np.log10( Vr_map ), vmin=2.5, vmax=4, cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent, alpha=0.4 )

 # --- attenuation
 vy_grid = grd['wind']['vy'][:,:]
 wind_map_above = np.ma.masked_array( wind_grid, vy_grid >= 0. )
 wind_map_below = np.ma.masked_array( wind_grid, vy_grid < 0. )
 AtMask = (vy_grid < 0.) & (np.log10(gas_map) < 9)
 wind_map_attenuated = np.ma.masked_array( wind_grid, AtMask )
 wind_map_rest = np.ma.masked_array( wind_grid, np.logical_not(AtMask) & (wind_grid <= wind_grid.min()) )
 im = ax.imshow( np.log10( wind_map_attenuated ), vmin=vwind[0], vmax=vwind[1], cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent, alpha=0.2 )
 im = ax.imshow( np.log10( wind_map_rest ), vmin=vwind[0], vmax=vwind[1], cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent, alpha=0.4 )
 """
 
 # --- dense ISM
 dense_map = np.ma.masked_array( gas_map, np.log10(gas_map) <= 8.5 )
 im = ax.imshow( np.log10( dense_map ), vmin=vgas[0], vmax=vgas[1], cmap=cmap_gas, interpolation='bicubic', origin='low', extent=extent, alpha=0.15 )
 
 # --- virial radius...
 for rfac in [1.,10.,100.]:
    ax.plot( rfac*xcirc, rfac*ycirc, ':', color='black', linewidth=1)
 #ax.plot( Rvir*xcirc, Rvir*ycirc, ':', color='black', linewidth=1)

 #papertxt = r'$\rm{Angl\'es-Alc\'azar \, et \, al.\, 2017, \, MNRAS, \, 472, \, L109}$'
 #papertxt = r'$\rm{Angl\'es}$-${\rmAlc\'azar \, et \, al.\, 2017, \, MNRAS, \, 472, \, L109}$'
 #ax.text( 0.01, 0.005, papertxt, {'color':'orange', 'fontsize':8}, transform=ax.transAxes, horizontalalignment='left', verticalalignment='bottom' )
     
 fig.savefig( outfile, dpi=150 )
 plt.close()


"""
ffmpeg_str = 'ffmpeg -f image2 -framerate 24 -start_number 0000 -i ' + outdir + '/tzgrid_%04d.png -r 24 -c:v h264 -crf 1 -pix_fmt yuv420p -y ' + outdir + '/' + movie_name
print '\n'+ffmpeg_str+'\n'
os.system(ffmpeg_str)
"""

print '\nDone!\n'


