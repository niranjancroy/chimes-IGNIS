#>>> execfile('bhfeedback_cenBHgrids_tiltedzoom.py')

# python bhfeedback_cenBHgrids_tiltedzoom.py /simons/scratch/dangles/FIRE/bhfeedback/  h113_HR_sn153/kick_s8e1v4p1  snap_ini snap_end > log/grid_kick_s8e1v4p1.log 2>&1 &

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
  print 'syntax: python bhfeedback_cenBHgrids_tiltedzoom.py  basedir sname snap_ini snap_end'
  sys.exit()


"""
basedir = '/simons/scratch/dangles/FIRE/bhfeedback/'
sname = 'h113_HR_sn277/spawn_s9e1v30f05_w100dt1000m100_TshVsh'
snap_ini = 1000
snap_end = 1001
"""

SkipIfDone = 0

npix = 512

nsnap_zoom = 1000
dx_min = 0.5
dx_max = 200.
dx_grid = np.logspace( np.log10(dx_min), np.log10(dx_max), nsnap_zoom )

m_wind_cut = 200.
bfr = 2.


simdir = basedir + sname

outdir = simdir + '/tzgrid'
if not os.path.exists(outdir):
  try:
    os.mkdir(outdir)
  except OSError:
    try:
      os.makedirs(outdir)
    except OSError:
      print '...could not create output directory: ' + outdir



grid_tmp = { 'snapnum':-1, 'redshift':-1., 'dx':-1,
             'gas':{ 'Npart':0, 'm':np.zeros([npix,npix]), 'T':np.zeros([npix,npix]), 'sfr':np.zeros([npix,npix]) }, 
             'wind':{ 'Npart':0, 'm':np.zeros([npix,npix]), 'T':np.zeros([npix,npix]), 'Vr':np.zeros([npix,npix]), 'vz':np.zeros([npix,npix]), 'ek':np.zeros([npix,npix]) }
           }



# --- rotation matrix ---
bhdata = daa.load_dict_from_hdf5(simdir + '/cenBH/cenBH_0001.hdf5')
#bhdata = daa.load_dict_from_hdf5(simdir + '/cenBH/cenBH_0032.hdf5')
Rref = 1.   # kpc
ind_Rcm = np.argmin( np.abs(bhdata['cum']['radius']-Rref) )
Lgas = bhdata['cum']['L_gas'][ind_Rcm,:]
zaxis = np.array([0.,0.,1.])
L_dir = Lgas / np.linalg.norm(Lgas)
rotax = np.cross(L_dir,zaxis)
rotangle = np.arccos(np.dot(L_dir,zaxis))
rotmatrix_z = daa.rotation_matrix( rotax, rotangle )
rotmatrix_tilted = daa.rotation_matrix( [1.,0.,0.], 0.5 )
rotmatrix = np.matmul( rotmatrix_tilted, rotmatrix_z )


for snapnum in range(snap_ini,snap_end,1):
#for snapnum in [1,100,200,300,400,500]:#,600,700,800,900,1000]:

  outfile = outdir + '/tzgrid_%04d.hdf5' % snapnum
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

  if snapnum - 1 >= dx_grid.size:
    print '  ... skipping snapnum', snapnum, '  grid size =', dx_grid.size
    continue
  dx = dx_grid[snapnum-1]

  grid = grid_tmp.copy()

  grid['snapnum'] = snapnum
  grid['redshift'] = header['redshift']
  grid['dx'] = dx


  # --- REFERENCE FRAME ---
  bhfile = simdir + '/cenBH/cenBH_%04d.hdf5' % snapnum
  if not os.path.exists( bhfile ):
    print ' --> snapshot ', snapnum, '  not available in bhdata! '
    continue
  bhdata = daa.load_dict_from_hdf5(bhfile)
  cen_pos = bhdata['cum']['COMstar_pos'][ind_Rcm,:] + bhdata['bh']['pos'][:]
  cen_vel = bhdata['cum']['COMstar_vel'][ind_Rcm,:] + bhdata['bh']['vel'][:]


  # --- GAS ---
  P = g.readsnap( simdir, snapnum, 0, cosmological=1 )
  r_gas = P['p'][:,:] - cen_pos[:]
  R_gas = np.sqrt((r_gas*r_gas).sum(axis=1))
  hubble_factor = daa.hubble_z( P['redshift'] )
  v_gas = P['v'][:,:] - cen_vel[:] +  hubble_factor * r_gas * UnitLength_in_cm  / UnitVelocity_in_cm_per_s
  Vr_gas = (v_gas*r_gas).sum(axis=1) / R_gas
  m = P['m'][:] * UnitMass_in_Msun
  T = g.gas_temperature(P['u'],P['ne'])
  sfr = P['sfr'][:]
  hsml = P['h'][:]
  #r_gas_rot_tmp = np.tensordot( rotmatrix, r_gas[:,:], axes=(1,1) ).T
  #r_gas_rot = np.tensordot( rotmatrix_tilted, r_gas_rot_tmp[:,:], axes=(1,1) ).T
  r_gas_rot = np.tensordot( rotmatrix, r_gas[:,:], axes=(1,1) ).T
  v_gas_rot = np.tensordot( rotmatrix, v_gas[:,:], axes=(1,1) ).T

  xrg = [ -dx, dx ]

  # --- normal gas (not wind) ---  
  ind = np.where( ( m > m_wind_cut ) & (np.abs(r_gas_rot[:,0]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,1]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,2]) <= xrg[1]*bfr) )[0]
  if ind.size >= 1:
    x_grid = r_gas_rot[ind,0]
    y_grid = r_gas_rot[ind,1]
    z_grid = r_gas_rot[ind,2]
    m_grid = m[ind]
    T_grid = T[ind]
    sfr_grid = sfr[ind]
    hsml_grid = daa.get_particle_hsml( x_grid, y_grid, z_grid, DesNgb=32 )
    
    grid['gas']['Npart'] = ind.size
    m_sMap, Tm_sMap, sfr_sMap = daa.make_2Dgrid( x_grid, z_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, weight3=sfr_grid, hsml=hsml_grid, pixels=npix )
    Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
    grid['gas']['m'][:,:] = Mgas_map.T
    grid['gas']['T'][:,:] = Tgas_map.T
    grid['gas']['sfr'][:,:] = sfr_sMap.T


  # --- wind particles ---
  ind = np.where( ( m <= m_wind_cut ) & (np.abs(r_gas_rot[:,0]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,1]) <= xrg[1]*bfr) & (np.abs(r_gas_rot[:,2]) <= xrg[1]*bfr) )[0]
  if ind.size >= 1:
    x_grid = r_gas_rot[ind,0]
    y_grid = r_gas_rot[ind,1]
    z_grid = r_gas_rot[ind,2]
    vz_grid = v_gas_rot[ind,2]
    Vr_grid = Vr_gas[ind]
    m_grid = m[ind]
    ek_grid = 0.5 * m_grid * Vr_grid**2.
    T_grid = T[ind]
    hsml_grid = daa.get_particle_hsml( x_grid, y_grid, z_grid, DesNgb=16 )

    grid['wind']['Npart'] = ind.size
    m_sMap, Tm_sMap = daa.make_2Dgrid( x_grid, z_grid, xrg, weight1=m_grid, weight2=T_grid*m_grid, hsml=hsml_grid, pixels=npix )
    Mgas_map, Tgas_map = daa.clip_2Dgrid( m_sMap, map2=Tm_sMap )
    grid['wind']['m'][:,:] = Mgas_map.T
    grid['wind']['T'][:,:] = Tgas_map.T
  
    ek_sMap, Vr_sMap, vz_sMap = daa.make_2Dgrid( x_grid, z_grid, xrg, weight1=ek_grid, weight2=Vr_grid*m_grid, weight3=vz_grid*m_grid, hsml=hsml_grid, pixels=npix )
    ek_map = daa.clip_2Dgrid( ek_sMap )
    Mgas_map, Vr_map  = daa.clip_2Dgrid( m_sMap, map2=Vr_sMap )
    Mgas_map, vz_map = daa.clip_2Dgrid( m_sMap, map2=vz_sMap )
    grid['wind']['ek'][:,:] = ek_map.T
    grid['wind']['Vr'][:,:] = Vr_map.T
    grid['wind']['vz'][:,:] = vz_map.T


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




