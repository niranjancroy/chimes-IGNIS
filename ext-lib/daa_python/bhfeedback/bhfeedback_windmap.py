#>>> execfile('bhfeedback_windmap.py')

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




basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
ICdir = 'h113_HR_sn152'
sname = 'spawn_res251_s9e1v30f05_w10dt100m1000_TshVsh'
stag = 'e1f05'

snapnum = int( 252 + 5./0.01 + 5./0.1 )
dx = 10.
SkipIfDone = 1

simdir = basedir + ICdir + '/' + sname

outdir = simdir + '/tzgrid'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir

#figdir = outdir
figdir = './bhfeedback/' + stag + '/tzgrid'
if not os.path.exists(figdir):
  try:
    os.mkdir(figdir)
  except OSError:
    try:
      os.makedirs(figdir)
    except OSError:
      print '...could not create output directory: ' + figdir



bhdata = daa.load_dict_from_hdf5(simdir + '/cenBH/centralBH.hdf5')
time_start = bhdata['time'][0]
t = bhdata['time'] - time_start
ind_snap = np.argmin( np.abs(bhdata['snapnum'][:] - snapnum) )
print 'time =', t[ind_snap], 'Myr', '  snapnum =', snapnum

Rref = 1.   # kpc
ind_Rcm = np.argmin( np.abs(bhdata['radial_bins']-Rref) )
Lgas = bhdata['cum']['L_gas'][0,ind_Rcm,:]
zaxis = np.array([0.,0.,1.])
L_dir = Lgas / np.linalg.norm(Lgas)
rotax = np.cross(L_dir,zaxis)
rotangle = np.arccos(np.dot(L_dir,zaxis))
rotmatrix_z = daa.rotation_matrix( rotax, rotangle )
rotmatrix_tilted = daa.rotation_matrix( [1.,0.,0.], 0.5 )
#rotmatrix = np.matmul( rotmatrix_tilted, rotmatrix_z )
rotmatrix = rotmatrix_z


#for snapnum in [352, 752, 802, 902]:
for snapnum in [902]:
 for dx in [ 0.5, 1., 5., 10., 50., 100. ]:

  filename = 'tzgrid_ns%04d_dx%.3f' % ( snapnum, dx )

  gridfile = outdir + '/' + filename  + '.hdf5'
  gridDone = bhlib.make_windgrid( simdir, snapnum, gridfile, dx=dx, rotmatrix=rotmatrix, m_wind=3000., npix=512, SkipIfDone=SkipIfDone, P=-1 )

  outfile = figdir + '/' + stag + '_' + filename  + '.png'
  plotDone, vmax = bhlib.make_windgrid_plot( gridfile, outfile, vmax=-999 )











print '\n\nDone!\n'













