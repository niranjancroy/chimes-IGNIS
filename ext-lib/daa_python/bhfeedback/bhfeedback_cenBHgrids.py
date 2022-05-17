# run from python command line
#>>> import sys
#>>> sys.argv = ['bhfeedback_cenBHgrids.py', 'h113_HR_sn153/kick_s8e1v4p1', 0, 10 ]
#>>> execfile('bhfeedback_cenBHgrids.py')

# python bhfeedback_cenBHgrids.py /simons/scratch/dangles/FIRE/bhfeedback/  h113_HR_sn153/kick_s8e1v4p1  snap_ini snap_end > log/grid_kick_s8e1v4p1.log 2>&1 &

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
from copy import deepcopy
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
  print 'syntax: python bhfeedback_cenBHgrids.py basedir sname snap_ini snap_end'
  sys.exit()


"""
basedir = '/simons/scratch/dangles/FIRE/bhfeedback/old/'
sname = 'h113_HR_sn153/nofback'
snap_ini = 1
snap_end = 2
"""

SkipIfDone = 1

npix = 512

dx0 = 100.
dx1 = 10.
dx2 = 1. 
dx3 = 0.1 

bfr = 2.


simdir = basedir + sname

outdir = simdir + '/grids'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir


grid_vars = { 'mgas0':{ 'dx':dx0, 'Npart':0, 'g':np.zeros([npix,npix]) },
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

grid_tmp = { 'snapnum':-1, 'redshift':-1., 
             'faceon':deepcopy(grid_vars), 'edgeon':deepcopy(grid_vars) }



# --- rotation matrix ---
"""
#bhdata = daa.load_dict_from_hdf5(simdir + '/cenBH/cenBH_0001.hdf5')
filelist, filestr = daa.file_list( simdir + '/cenBH/cenBH', ext='.hdf5' )
bhdata = daa.load_dict_from_hdf5(simdir + '/cenBH/cenBH_%04d.hdf5' % filelist[0])
Rref = 1.   # kpc
ind_Rcm = np.argmin( np.abs(bhdata['cum']['radius']-Rref) )
Lgas = bhdata['cum']['L_gas'][ind_Rcm,:]
"""
bhdata = daa.load_dict_from_hdf5(simdir + '/cenBH/centralBH.hdf5')
Rref = 1.     # kpc
ref_snap = 0  # use first snapshot available to define rotation axis
ind_Rcm = np.argmin( np.abs(bhdata['radial_bins']-Rref) )
Lgas = bhdata['cum']['L_gas'][ref_snap,ind_Rcm,:]
zaxis = np.array([0.,0.,1.])
L_dir = Lgas / np.linalg.norm(Lgas)
rotax = np.cross(L_dir,zaxis)
rotangle = np.arccos(np.dot(L_dir,zaxis))
rotmatrix = daa.rotation_matrix( rotax, rotangle )


for snapnum in range(snap_ini,snap_end,1):
#for snapnum in [0,1]:

  outfile = outdir + '/mgrid_%04d.hdf5' % snapnum
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


  # --- REFERENCE FRAME ---
  bhfile = simdir + '/cenBH/cenBH_%04d.hdf5' % snapnum
  if not os.path.exists( bhfile ):
    print ' --> snapshot ', snapnum, '  not available in bhdata! '
    continue
  bhdata = daa.load_dict_from_hdf5(bhfile)
  cen_pos = bhdata['cum']['COMstar_pos'][ind_Rcm,:] + bhdata['bh']['pos'][:]


  # --- GAS ---
  P = g.readsnap( simdir, snapnum, 0, cosmological=1 )
  r_gas = P['p'][:,:] - cen_pos[:]
  r_gas_rot = np.tensordot( rotmatrix, r_gas[:,:], axes=(1,1) ).T
  m = P['m'][:] * UnitMass_in_Msun
  T = g.gas_temperature(P['u'],P['ne'])
  hsml = P['h'][:]
  for im, dx in zip( ['0','1','2','3'], [dx0,dx1,dx2,dx3] ):
     xrg = [ -dx, dx ]
     ind = np.where( (np.abs(r_gas_rot[:,0]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,1]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,2]) <= xrg[1]*bfr) )[0]
     if ind.size <= 1:
       continue
     x_grid = r_gas_rot[ind,0]
     y_grid = r_gas_rot[ind,1]
     z_grid = r_gas_rot[ind,2]
     m_grid = m[ind]
     T_grid = T[ind]
     hsml_grid = hsml[ind]

     grid['faceon']['mgas'+im]['Npart'] = ind.size
     grid['faceon']['T'+im]['Npart'] = ind.size
     m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, y_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
     Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
     grid['faceon']['mgas'+im]['g'][:,:] = Mgas_map.T
     grid['faceon']['T'+im]['g'][:,:] = Tgas_map.T

     grid['edgeon']['mgas'+im]['Npart'] = ind.size
     grid['edgeon']['T'+im]['Npart'] = ind.size
     m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, z_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
     Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
     grid['edgeon']['mgas'+im]['g'][:,:] = Mgas_map.T
     grid['edgeon']['T'+im]['g'][:,:] = Tgas_map.T


  

  # --- STARS ---  
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

  print '\n... writing', outfile
  if os.path.isfile(outfile ):
    os.remove(outfile)
  file = h5py.File(outfile, 'w')
  for grname in grid.keys():
    if type(grid[grname]) is dict:
      group = file.create_group( grname )
      for keyname in grid[grname].keys():
         if type(grid[grname][keyname]) is dict:
           group2 = group.create_group( keyname )
           for keyname2 in grid[grname][keyname].keys():
              group2.create_dataset(keyname2, data=grid[grname][keyname][keyname2])
         else:
           group.create_dataset(keyname, data=grid[grname][keyname])
    else:
      file.create_dataset(grname, data=grid[grname])
  file.close()



print '\nDone!\n'




