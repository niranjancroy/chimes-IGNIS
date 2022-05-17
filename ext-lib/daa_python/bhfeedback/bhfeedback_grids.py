# run from python command line
#>>> import sys
#>>> sys.argv = ['bhfeedback_grids.py', 'h113_HR_sn153/kick_s8e1v4p1', 0, 10 ]
#>>> execfile('bhfeedback_grids.py')

# python bhfeedback_grids.py  h113_HR_sn153/kick_s8e1v4p1  snap_ini snap_end > log/grid_kick_s8e1v4p1.log 2>&1 &

import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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
  print 'syntax: python bhfeedback_grids.py sname snap_ini snap_end'
  sys.exit()



SkipIfDone = 1

npix = 512

dx0 = 100.
dx1 = 10.
dx2 = 1. 
dx3 = 0.1 

bfr = 2.


simdir = '/simons/scratch/dangles/FIRE/bhfeedback/old/' + sname

outdir = simdir + '/grids'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir


grid_tmp = { 'snapnum':-1, 'redshift':-1., 
             'mgas0':{ 'dx':dx0, 'Npart':0, 'g':np.zeros([npix,npix]) },
             'mgas1':{ 'dx':dx1, 'Npart':0, 'g':np.zeros([npix,npix]) }, 
             'mgas2':{ 'dx':dx2, 'Npart':0, 'g':np.zeros([npix,npix]) }, 
             'mgas3':{ 'dx':dx3, 'Npart':0, 'g':np.zeros([npix,npix]) },
             'T0':{ 'dx':dx0, 'Npart':0, 'g':np.zeros([npix,npix]) },
             'T1':{ 'dx':dx1, 'Npart':0, 'g':np.zeros([npix,npix]) },
             'T2':{ 'dx':dx2, 'Npart':0, 'g':np.zeros([npix,npix]) },
             'T3':{ 'dx':dx3, 'Npart':0, 'g':np.zeros([npix,npix]) } #,
             #'mstar1':{ 'dx':dx1, 'Npart':np.zeros(dtype=np.int64), 'g':np.zeros([npix,npix]) }, 
             #'mstar2':{ 'dx':dx2, 'Npart':np.zeros(dtype=np.int64), 'g':np.zeros([npix,npix]) }, 
             #'mstar3':{ 'dx':dx3, 'Npart':np.zeros(dtype=np.int64), 'g':np.zeros([npix,npix]) } 
        }


for snapnum in range(snap_ini,snap_end,1):
#for snapnum in [0,1]:

  outfile = outdir + '/grid_%03d.hdf5' % snapnum
  if SkipIfDone and os.path.exists( outfile ):
    print "\nfile exists:  %s\n" % outfile
    continue

  try:
    header = g.readsnap( simdir, snapnum, 0, header_only=1)
  except Exception:
    print simdir, '  snapnum=', snapnum, ' not found!'
    continue

  if header['k'] == -1:
    print simdir, '  snapnum=', snapnum, ' not found!'
    continue

  grid = grid_tmp.copy()

  grid['snapnum'] = snapnum
  grid['redshift'] = header['redshift']

  # --- central black hole ---
  try:
    P = g.readsnap( simdir, snapnum, 5, cosmological=1 )
  except KeyError:
    P = g.readsnap( simdir, snapnum, 5, cosmological=1, skip_bh=1 )

  if P['k'] == -1:
    print 'no snapshot available: ', snapnum
    continue
  if P['id'].size == 1:
    ibh = 0
  else:
    print 'P.size in snapnum ', snapnum, '   ...', P['id'].size
    ibh = np.argmax( P['m'] )

  cen_pos = P['p'][ibh,:]

  if 'mbh' in P.keys():
     print 'ns =', snapnum, '  redshift =', P['redshift'], '  m =', P['m'][ibh]*UnitMass_in_Msun, '  mbh =', P['mbh'][ibh]*UnitMass_in_Msun, \
           '  edd =', P['mdot'][ibh]*UnitMdot_in_Msun_per_yr/daa.MdotEddington(P['mbh'][ibh]*UnitMass_in_Msun), '\n'
  else:
     print 'ns =', snapnum, '  redshift =', P['redshift'], '  m =', P['m'][ibh]*UnitMass_in_Msun

  
  # --- gas ---
  P = g.readsnap( simdir, snapnum, 0, cosmological=1 )
  x = P['p'][:,0] - cen_pos[0]
  y = P['p'][:,1] - cen_pos[1]
  z = P['p'][:,2] - cen_pos[2]
  m = P['m'][:] * UnitMass_in_Msun
  T = g.gas_temperature(P['u'],P['ne'])
  hsml = P['h'][:]
  for im, dx in zip( ['0','1','2','3'], [dx0,dx1,dx2,dx3] ):
     xrg = [ -dx, dx ]
     ind = np.where( (np.abs(x) <= xrg[1]*bfr) & (np.abs(y) <= xrg[1]*bfr) & (np.abs(z) <= xrg[1]*bfr) )[0]
     grid['mgas'+im]['Npart'] = ind.size
     grid['T'+im]['Npart'] = ind.size
     if ind.size <= 1:
       continue
     x_grid = x[ind] 
     y_grid = y[ind] 
     m_grid = m[ind]
     T_grid = T[ind]
     hsml_grid = hsml[ind]
     m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, y_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
     Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
     grid['mgas'+im]['g'][:,:] = Mgas_map.T
     grid['T'+im]['g'][:,:] = Tgas_map.T
  

  # --- stars  
  """
  P = g.readsnap( simdir, snapnum, 4, cosmological=1 )
  x = P['p'][:,0] - cen_pos[0]
  y = P['p'][:,1] - cen_pos[1]
  z = P['p'][:,2] - cen_pos[2]
  m = P['m'][:] * UnitMass_in_Msun
  hsml = daa.get_particle_hsml( x, y, z, DesNgb=32 )
  for im, dx in zip( ['1','2','3'], [dx1,dx2,dx3] ):
     xrg = [ -dx, dx ]
     ind = np.where( (np.abs(x) <= xrg[1]*bfr) & (np.abs(y) <= xrg[1]*bfr) & (np.abs(z) <= xrg[1]*bfr) )[0]
     grid['mstar'+im]['Npart'] = ind.size
     if ind.size <= 1:
       continue
     x_grid = x[ind]
     y_grid = y[ind]
     m_grid = m[ind]
     hsml_grid = hsml[ind]
     m_sMap = daa.make_2Dgrid( x_grid, y_grid, xrg, weight1=m_grid, hsml=hsml_grid, pixels=npix )
     Mstar_map = daa.clip_2Dgrid( m_sMap )
     grid['mstar'+im]['g'][:,:] = Mstar_map.T
  """


  ############################
   ### WRITE DATA TO FILE ###
  ############################

  if os.path.isfile(outfile ):
    os.remove(outfile)
  file = h5py.File(outfile, 'w')
  for grname in grid.keys():
     if type(grid[grname]) is dict:
       group = file.create_group( grname )
       for keyname in grid[grname].keys():
          group.create_dataset(keyname, data=grid[grname][keyname])
     else:
       file.create_dataset(grname, data=grid[grname])
  file.close()



print '\nDone!\n'




