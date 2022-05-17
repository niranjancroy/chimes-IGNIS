#>>> execfile('bhfeedback_cavitymaps.py')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
import os as os
from shutil import copyfile
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
import daa_bhlib as bhlib
from daa_constants import *

import imp
imp.reload(daa)
imp.reload(bhlib)



sim_list = [ 'nof_s8e1_n128',
             'spawn_res251_s9e01v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f033_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f02_w10dt100m1000_TshVsh',
             #'spawn_res251_s9e1v30f011_w10dt100m1000_TshVsh',
             'spawn_res251_s9e1v30f009_w10dt100m1000_TshVsh'
           ]

tag_list = [ 'nof',
             'e01f05',
             'e1f05',
             'e1f033',
             'e1f02',
             #'e1f011',
             'e1f009'
           ]

basedir = '/simons/scratch/dangles/FIRE/bhfeedback/h113_HR_sn152/'
snap_start = 251

outdir = './bhfeedback/allplots'

skipreading = 0


dx = 1.
times = [ 0.2, 1., 5., 10., 15., 20., 25., 30., 35. ]
SkipIfDone = 1


nsim = len(sim_list)
ntimes = len(times)

if skipreading==0:

 sim = {}
 for sname, stag in zip(sim_list,tag_list):   
    bh_file = basedir + sname + '/cenBH/centralBH.hdf5'
    print '...reading ' + bh_file
    bh = daa.load_dict_from_hdf5( bh_file )
    sim[stag] = bh
 inum = np.where(sim['nof']['snapnum'] == snap_start)[0]
 time_start = sim['nof']['time'][ inum[0] ]
 rotmatrix = bhlib.bhrotmatrix( sim['nof'], Rref=1., inum=inum[0] )

 grids = {}
 for sname, stag in zip(sim_list,tag_list):
  grids[stag] = {}
  simdir = basedir + sname
  t = sim[stag]['time'] - time_start
  for i, thistime in enumerate(times):
    indt = np.argmin( np.abs(t - thistime) )
    snapnum = sim[stag]['snapnum'][indt]
    
    #grid_file = basedir + sname + '/grids/mgrid_%04d.hdf5' % snapnum
    #print '...reading ' + grid_file
    #try:
    #  grd = daa.load_dict_from_hdf5( grid_file )
    #except IOError:
    #  print '        --> not available...'
    #  continue
    #grd['thistime'] = t[indt]
    #grids[stag][str(i)] = grd
    filename = 'tzgrid_ns%04d_dx%.3f' % ( snapnum, dx )
    gridfile = simdir + '/tzgrid' + '/' + filename  + '.hdf5'
    print '...reading ', gridfile
    gridDone = bhlib.make_windgrid( simdir, snapnum, gridfile, dx=dx, rotmatrix=rotmatrix, center='BH', m_wind=3000., DesNgb=48, npix=512, SkipIfDone=SkipIfDone, P=-1 )
    if gridDone:
      grd = daa.load_dict_from_hdf5( gridfile )
    else:
      print '        --> not available...'
      continue
    grids[stag][str(i)] = grd['faceon']


lw_ax = 1.5

vmax = 10.5
vmin = 7

figname = 'bhfeedback_cavitymaps.pdf'
print '\n', figname
fxsz = 10
xyrat = float(nsim) / ntimes
fysz = fxsz * xyrat
ds = 0.01
fig = plt.figure( figsize=(fxsz,fysz) )
gs = GridSpec(nsim,ntimes)
gs.update(left=ds, right=1-ds, bottom=ds, top=1-ds, hspace=0., wspace=0)

for n, stag in enumerate(tag_list):
 for i in range(ntimes):
   ax = plt.subplot(gs[n,i])
   for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_visible(False)
   for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(lw_ax)
   if str(i) not in grids[stag].keys():
     ax.set_xlim(xrg); ax.set_ylim(xrg); ax.set_aspect('equal')
     continue
   """
   thisgrid = grids[stag][str(i)]
   thistime = thisgrid['time']
   grd = thisgrid['faceon']['mgas2']
   print stag, thistime, thisgrid['snapnum']
   dx = grd['dx']
   xrg = [ -dx, dx ]
   yrg = xrg
   extent = [xrg[0],xrg[1],xrg[0],xrg[1]]  
   ax.set_xlim(xrg); ax.set_ylim(xrg); ax.set_aspect('equal')
   im = ax.imshow( np.log10( grd['g'][:,:] ), vmin=vmin, vmax=vmax, cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent )
   """
   grd = grids[stag][str(i)]
   thistime = grd['time'] - time_start
   if np.abs(thistime-times[i]) > 1.:
     continue
   print stag, thistime, grd['snapnum']
   dx = grd['dx']
   xrg = [ -dx, dx ]
   yrg = xrg
   extent = [xrg[0],xrg[1],xrg[0],xrg[1]]
   ax.set_xlim(xrg); ax.set_ylim(xrg); ax.set_aspect('equal')
   gasmap = grd['gas']['m'][:,:] + grd['wind']['m'][:,:]
   im = ax.imshow( np.log10( gasmap ), vmin=vmin, vmax=vmax, cmap=cm.Spectral, interpolation='bicubic', origin='low', extent=extent )

fig.savefig( outdir + '/' + figname, dpi=600 )
plt.close()






print '\n\nDone!\n'













